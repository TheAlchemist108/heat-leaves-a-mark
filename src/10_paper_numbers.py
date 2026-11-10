"""Write every number and table quoted in the paper from the results files, so the text never
carries a hand-typed estimate.

Usage:
  python src/10_paper_numbers.py <processed_dir> <results_dir> <paper_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def f2(x):
    return f"{x:.2f}".replace("-", "$-$") if pd.notna(x) else "--"


def f3(x):
    return f"{x:.3f}".replace("-", "$-$") if pd.notna(x) else "--"


def pct(x):
    return f"{100 * (np.exp(x) - 1):.0f}" if pd.notna(x) else "--"


def pval(p):
    if pd.isna(p):
        return "--"
    return "$<$0.001" if p < 0.001 else f"{p:.3f}"


def main(proc, res, paper):
    proc, res, paper = Path(proc), Path(res), Path(paper)
    paper.mkdir(parents=True, exist_ok=True)
    m = {}
    info = json.loads((res / "run_info.json").read_text())
    m["NEpisodes"] = info["n_insample_onsets"]
    m["NPulls"] = len(info["pulls"])
    m["NDraws"] = f"{info['n_draws']:,}"
    m["NPoolWeeks"] = info["placebo_pool_weeks"]
    m["AnchorValidityP"] = pval(info["anchor_validity_p"])

    thr = pd.read_csv(proc / "heat_thresholds.csv")
    t = thr[thr.percentile == "p95"].set_index("series").threshold_degC
    m["ThrUK"], m["ThrEngland"], m["ThrScotland"] = f"{t['UK']:.1f}", f"{t['England']:.1f}", f"{t['Scotland']:.1f}"
    m["ThrWales"], m["ThrNI"] = f"{t['Wales']:.1f}", f"{t['Northern Ireland']:.1f}"
    ep = pd.read_csv(proc / "heat_episodes.csv")
    ep = ep[(ep.series == "UK") & (ep.percentile == "p95") & ep.in_sample]
    m["MaxDose"] = f"{ep.dose_degree_days.max():.1f}"
    m["MaxDoseEpisode"] = pd.Timestamp(ep.loc[ep.dose_degree_days.idxmax(), "start"]).strftime("%-d %B %Y")
    m["DoseJulyTwentyTwo"] = f"{ep[ep.start == '2022-07-17'].dose_degree_days.iloc[0]:.1f}"

    cal = pd.read_csv(proc / "era5_calibration.csv")
    nat = cal[cal.series != "UK"]
    m["CalRMin"], m["CalRMax"] = f"{nat.r.min():.3f}", f"{nat.r.max():.3f}"

    sp = pd.read_csv(proc / "trends_sparsity.csv")
    gb = sp[sp.geo == "GB"].set_index("term").share_le1
    m["SparseHeatwave"] = f"{gb['heatwave']:.2f}"
    cop = sp[sp.term.isin(["air conditioning", "portable air conditioner", "heat stroke"])]
    m["SparseCopingMin"], m["SparseCopingMax"] = f"{cop.share_le1.min():.2f}", f"{cop.share_le1.max():.2f}"
    rel = pd.read_csv(proc / "trends_pull_reliability.csv")
    m["PullCorrMin"] = f"{rel[rel.geo == 'GB'].mean_between_pull_corr.min():.2f}"
    wuk = pd.read_csv(proc / "wikipedia_uk_sparsity.csv")
    m["WikiUKPass"] = int(wuk.passes.sum())

    size = pd.read_csv(res / "size_check.csv").set_index("test").rejection_rate_5pct
    m["SizeAOneMin"] = f"{100 * size[[i for i in size.index if i.startswith('A1')]].min():.0f}"
    m["SizeAOneMax"] = f"{100 * size[[i for i in size.index if i.startswith('A1')]].max():.0f}"
    m["SizeATwo"] = f"{100 * size[[i for i in size.index if i.startswith('A2')]].iloc[0]:.0f}"

    a1 = pd.read_csv(res / "analysis1_estimates.csv").set_index("outcome")
    for key, n in [("Sal", "trends_dev2_salience"), ("Cop", "trends_dev2_coping"), ("Con", "trends_dev2_concern"),
                   ("ConPlan", "trends_plan_concern"), ("WikiSal", "wiki_global_salience"),
                   ("WikiCon", "wiki_global_concern"), ("WikiCop", "wiki_global_coping")]:
        r = a1.loc[n]
        m[f"Imp{key}"], m[f"Imp{key}Pct"], m[f"Imp{key}RI"] = f2(r.impact), pct(r.impact), pval(r.impact_ri_p)
        m[f"Imp{key}MDE"] = f2(r.impact_mde)
        m[f"Res{key}"], m[f"Res{key}RI"] = f2(r.residue_vs_baseline), pval(r.residue_vs_baseline_ri_p)
        m[f"Dec{key}"] = f2(r.decay_vs_baseline)
        m[f"Leads{key}"] = pval(r.leads_wald_p)
    for key, n in [("CarIns", "term_car insurance"), ("Recipes", "term_recipes")]:
        m[f"Imp{key}"], m[f"Imp{key}RI"] = f2(a1.loc[n, "impact"]), pval(a1.loc[n, "impact_ri_p"])

    ex = pd.read_csv(res / "exploratory_vs_baseline.csv")
    exr = ex[ex.estimand == "residue_vs_baseline"].set_index("outcome")
    for key, n in [("Sal", "trends_dev2_salience"), ("Cop", "trends_dev2_coping"), ("Con", "trends_dev2_concern")]:
        m[f"Res{key}SECal"] = f2(exr.loc[n, "se_calibrated"])
        m[f"Res{key}Tost"] = pval(exr.loc[n, "p_tost_calibrated"])

    rob = pd.read_csv(res / "robustness.csv")
    py = rob[rob.check.str.startswith("placebo_year") & (rob.estimand == "impact")]
    for key, n in [("Sal", "trends_dev2_salience"), ("Cop", "trends_dev2_coping"), ("Con", "trends_dev2_concern")]:
        v = py[py.outcome == n].estimate
        m[f"PY{key}Min"], m[f"PY{key}Max"] = f2(v.min()), f2(v.max())
        m[f"PY{key}MaxAbs"] = f2(v.abs().max())

    a2 = pd.read_csv(res / "analysis2_estimates.csv")
    for key, fam in [("Cop", "coping"), ("Con", "concern"), ("Sal", "salience")]:
        r = a2[(a2.rule == "dev2") & (a2.family == fam) & ~a2.habituation_model].iloc[0]
        m[f"DR{key}"], m[f"DR{key}RI"] = f2(r.short_run), pval(r.short_run_ri_p)
        h = a2[(a2.rule == "dev2") & (a2.family == fam) & a2.habituation_model].iloc[0]
        m[f"Hab{key}"], m[f"Hab{key}RI"] = f2(h.habituation), pval(h.habituation_ri_p)

    hyp = pd.read_csv(res / "hypotheses.csv")
    m["NPresentCal"] = int((hyp[hyp.rule == "dev2"].result_calibrated.str.startswith("effect present")).sum())
    pre = hyp[(hyp.rule == "dev2") & hyp.result.str.startswith("effect present")]
    m["NPresentPrereg"] = len(pre)
    m["PresentPreregIDs"] = " and ".join(pre.id) if len(pre) else "none"
    m["NTestablePlan"] = int(hyp[hyp.rule == "plan"].testable.sum())

    # p-values with their own relation sign, for use as $p\XEq$ in the text
    for k in [k for k in list(m) if k.endswith("RI") or k.startswith("Leads") or k == "AnchorValidityP"]:
        v = str(m[k])
        m[k + "Eq"] = "<0.001" if v.startswith("$<$") else ("=" + v if v != "--" else "=\\text{--}")
    macro = "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in m.items())
    (paper / "numbers.tex").write_text("% Generated by src/10_paper_numbers.py. Do not edit.\n" + macro + "\n")

    # Table: hypotheses
    rows = []
    for _, r in hyp.iterrows():
        short = {"not testable (no series passes the sparsity rule)": "not testable",
                 "reset (equivalent to zero)": "reset", "effect present (positive)": "present (+)",
                 "effect present (negative)": "present ($-$)", "effect present (opposite sign)": "present (opp.)",
                 "no effect detected": "not detected"}
        rows.append(" & ".join([
            "Plan" if r.rule == "plan" else "Dev.~2", r.id, f2(r.estimate),
            pval(r.p_model), pval(r.p_ri), pval(r.p_calibrated_holm),
            short.get(r.result, r.result), short.get(r.result_calibrated, r.result_calibrated)]) + r" \\")
    tab = [r"\begin{tabular}{llrrrrll}", r"\toprule",
           r"Outcomes & Hyp. & Estimate & $p$ model & $p$ RI & $p$ cal., Holm & Pre-registered rule & Calibrated \\",
           r"\midrule"]
    for i, row in enumerate(rows):
        if i == 7:
            tab.append(r"\midrule")
        tab.append(row)
    tab += [r"\bottomrule", r"\end{tabular}"]
    (paper / "tab_hypotheses.tex").write_text("\n".join(tab) + "\n")

    # Table: robustness of the onset impact
    checks = [("main", "Main specification"), ("drop_concurrent_events", "Drop COP, election, state-funeral weeks"),
              ("threshold_p90", "Threshold: 90th percentile"), ("threshold_p975", "Threshold: 97.5th percentile"),
              ("metoffice_dates", "Met Office--documented dates"),
              ("placebo_year_-2", "Placebo: dates shifted $-2$ years"), ("placebo_year_-1", "Placebo: shifted $-1$ year"),
              ("placebo_year_+1", "Placebo: shifted $+1$ year"), ("placebo_year_+2", "Placebo: shifted $+2$ years")]
    outs = ["trends_dev2_salience", "trends_dev2_coping", "trends_dev2_concern", "term_car insurance"]
    tab = [r"\begin{tabular}{lrrrr}", r"\toprule",
           r"Check & Salience & Coping & Concern & Placebo (car insurance) \\", r"\midrule"]
    for c, lab in checks:
        vals = []
        for o in outs:
            if c == "main":
                vals.append(f2(a1.loc[o, "impact"]))
            else:
                v = rob[(rob.check == c) & (rob.outcome == o) & (rob.estimand == "impact")].estimate
                vals.append(f2(v.iloc[0]) if len(v) else "--")
        tab.append(f"{lab} & " + " & ".join(vals) + r" \\")
        if c == "metoffice_dates":
            tab.append(r"\midrule")
    tab += [r"\bottomrule", r"\end{tabular}"]
    (paper / "tab_robustness.tex").write_text("\n".join(tab) + "\n")
    print(f"wrote {len(m)} numbers and 2 tables to {paper}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
