"""Build outcome panels from the raw Google Trends and Wikipedia files (PREREGISTRATION.md, Section 6).

Google Trends (6.1)
  - values reported as "<1" (and zeros) are coded 0.5
  - log ratio to the anchor: r_j = ln(max(x_j, 0.5)) - ln(max(x_dentist, 0.5)), within each request
  - outcome = mean of r_j over all available pulls; SD across pulls and between-pull correlation reported
  - "climate change" is taken from R1; its R3 copy is used only for the cross-request consistency check
  - sparsity rule: a term-by-geo series is dropped if more than 20% of its weeks have a pull-averaged
    raw value of 1 or less; a geo is dropped if its anchor fails the same rule
  - family indices: mean of the member terms' log ratios (national: terms passing in GB; four-nation
    analysis: terms passing in all four nations)
  - deviation 2 (docs/DEVIATIONS.md): a second set of family indices ("dev2") uses all member terms with
    no sparsity rule; floor_share is the share of member terms at the 0.5 floor in each week
  - alternative anchor ("plumber"): r_j + ln(dentist / plumber) from R3, used only if the anchor
    validity check fails (Section 8.5)
Wikipedia (6.2)
  - global: log of weekly (Sunday-Saturday) total views
  - UK-only: an article is dropped if more than 20% of its days are missing; log of the weekly mean of
    the days present (missing days are below Wikimedia's release threshold)

Usage:
  python src/06_build_outcomes.py <data_dir> <out_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

FAMILIES = {
    "salience": ["heatwave"],
    "coping": ["air conditioning", "portable air conditioner", "heat stroke"],
    "concern": ["climate change", "global warming", "net zero"],
    "transition": ["heat pump"],
    "placebo": ["car insurance", "recipes"],
}
TERM_SOURCE = {  # request each term is taken from
    "heatwave": "R1", "climate change": "R1", "air conditioning": "R1", "car insurance": "R1",
    "global warming": "R2", "net zero": "R2", "portable air conditioner": "R2", "heat stroke": "R2",
    "heat pump": "R3", "recipes": "R3", "plumber": "R3",
}
GEOS = ["GB", "GB-ENG", "GB-SCT", "GB-WLS", "GB-NIR"]
NATION_GEOS = {"GB-ENG": "England", "GB-SCT": "Scotland", "GB-WLS": "Wales", "GB-NIR": "Northern Ireland"}
ANCHOR = "dentist"
SPARSE_SHARE = 0.20
WIKI_FAMILIES = {
    293270: "salience", 319671: "salience",
    7221088: "coping", 24999087: "coping", 12877572: "coping",
    5042951: "concern", 10031054: "concern", 10031822: "concern",
    68316: "transition", 5134: "placebo",
}


def code(fv, v):
    return 0.5 if fv == "<1" else max(float(v), 0.5)


def read_trends(raw_dir):
    """Long table: pull, kind, request, geo, window, date, term, x (coded raw value), r (log ratio)."""
    rows = []
    for pull_dir in sorted(Path(raw_dir).glob("pull*")):
        pull = int(pull_dir.name[4:])
        for f in sorted(pull_dir.glob("*.json")):
            d = json.loads(f.read_text())
            job, tl = d["job"], d["response"]["default"]["timelineData"]
            terms = job["terms"]
            ia = terms.index(ANCHOR)
            for t in tl:
                if t.get("isPartial"):
                    continue
                date = pd.to_datetime(int(t["time"]), unit="s").normalize()
                xs = [code(fv, v) for fv, v in zip(t["formattedValue"], t["value"])]
                for j, term in enumerate(terms):
                    rows.append((pull, job["kind"], job["request"], job["geo"], job["window"], date, term,
                                 xs[j], np.log(xs[j]) - np.log(xs[ia])))
    return pd.DataFrame(rows, columns=["pull", "kind", "request", "geo", "window", "date", "term", "x", "r"])


def pull_summary(df):
    g = df.groupby(["kind", "request", "geo", "date", "term"])
    out = g.agg(x_mean=("x", "mean"), r_mean=("r", "mean"), r_sd=("r", "std"), n_pulls=("pull", "nunique"))
    wide = df.pivot_table(index=["kind", "request", "geo", "date", "term"], columns="pull", values="r")
    wide.columns = [f"r_pull{c}" for c in wide.columns]
    return out.join(wide).reset_index()


def main(data_dir, out_dir):
    data, out = Path(data_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    long = read_trends(data / "raw" / "trends")
    pulls = sorted(long.pull.unique())
    summ = pull_summary(long)
    weekly = summ[summ.kind == "weekly"]
    daily = summ[summ.kind == "daily"]

    # sparsity rule (weekly, per term x geo, using the request each term is taken from)
    sp = []
    for geo in GEOS:
        a = weekly[(weekly.geo == geo) & (weekly.request == "R1") & (weekly.term == ANCHOR)]
        sp.append({"geo": geo, "term": ANCHOR, "request": "R1",
                   "share_le1": round(float((a.x_mean <= 1).mean()), 4)})
        for term, req in TERM_SOURCE.items():
            s = weekly[(weekly.geo == geo) & (weekly.request == req) & (weekly.term == term)]
            sp.append({"geo": geo, "term": term, "request": req,
                       "share_le1": round(float((s.x_mean <= 1).mean()), 4)})
    sp = pd.DataFrame(sp)
    sp["passes"] = sp.share_le1 <= SPARSE_SHARE
    anchor_ok = sp[sp.term == ANCHOR].set_index("geo").passes
    sp["geo_anchor_passes"] = sp.geo.map(anchor_ok)
    sp.to_csv(out / "trends_sparsity.csv", index=False)

    # weekly term series (dentist anchor) and alternative anchor
    terms = []
    for term, req in TERM_SOURCE.items():
        s = weekly[(weekly.request == req) & (weekly.term == term)].copy()
        terms.append(s)
    terms = pd.concat(terms)
    dp = weekly[(weekly.request == "R3") & (weekly.term == "plumber")][["geo", "date", "r_mean"]]
    dp = dp.rename(columns={"r_mean": "ln_plumber_minus_dentist"})
    terms = terms.merge(dp, on=["geo", "date"], how="left")
    terms["r_mean_plumber_anchor"] = terms.r_mean - terms.ln_plumber_minus_dentist
    for p in pulls:
        if f"r_pull{p}" in terms:
            dpp = weekly[(weekly.request == "R3") & (weekly.term == "plumber")][["geo", "date", f"r_pull{p}"]]
            terms = terms.merge(dpp.rename(columns={f"r_pull{p}": f"_pl{p}"}), on=["geo", "date"], how="left")
            terms[f"r_pull{p}_plumber_anchor"] = terms[f"r_pull{p}"] - terms[f"_pl{p}"]
            terms = terms.drop(columns=f"_pl{p}")
    terms = terms.merge(sp[["geo", "term", "passes", "geo_anchor_passes"]], on=["geo", "term"], how="left")
    terms.drop(columns=["kind"]).to_csv(out / "trends_weekly_terms.csv", index=False)

    # family indices
    fam_rows = []
    for col_suffix in [""] + [f"_pull{p}" for p in pulls]:
        for anchor in ["dentist", "plumber"]:
            col = ("r_mean" if col_suffix == "" else f"r{col_suffix}") + ("" if anchor == "dentist" else "_plumber_anchor")
            if col not in terms:
                continue
            for rule in ["plan", "dev2"]:  # dev2: no sparsity rule (docs/DEVIATIONS.md, entry 2)
                for scope in ["national", "four_nation"]:
                    for fam, members in FAMILIES.items():
                        if rule == "dev2":
                            keep = {"GB": members} if scope == "national" else {g: members for g in NATION_GEOS}
                        elif scope == "national":
                            keep = {"GB": [m for m in members if sp.query("geo=='GB' and term==@m").passes.all()]}
                        else:
                            ok = [m for m in members if sp[(sp.term == m) & sp.geo.isin(NATION_GEOS)].passes.all()
                                  and anchor_ok[list(NATION_GEOS)].all()]
                            keep = {g: ok for g in NATION_GEOS}
                        for geo, ms in keep.items():
                            if not ms:
                                continue
                            sub = terms[(terms.geo == geo) & terms.term.isin(ms)]
                            s = sub.groupby("date")[col].mean()
                            fl = sub.assign(fl=sub.x_mean <= 0.5).groupby("date").fl.mean()
                            for date, v in s.items():
                                fam_rows.append((rule, scope, geo, fam, "|".join(ms), anchor,
                                                 "mean" if col_suffix == "" else col_suffix[1:], date, v, fl[date]))
    fam = pd.DataFrame(fam_rows, columns=["rule", "scope", "geo", "family", "terms", "anchor", "pull", "date",
                                          "value", "floor_share"])
    fam.to_csv(out / "trends_weekly_families.csv", index=False)

    # anchor validity series and cross-request check
    av = weekly[(weekly.geo == "GB") & (weekly.request == "R3") & (weekly.term == "plumber")]
    av = av[["date", "r_mean"]].assign(value=lambda d: -d.r_mean)[["date", "value"]]
    av.to_csv(out / "anchor_validity_GB.csv", index=False)
    cc = weekly[(weekly.term == "climate change")].pivot_table(index=["geo", "date"], columns="request", values="r_mean")
    pd.DataFrame({"geo": [g for g in GEOS],
                  "corr_R1_R3_climate_change": [round(float(cc.loc[g].corr().loc["R1", "R3"]), 4) for g in GEOS]}
                 ).to_csv(out / "trends_cross_request_check.csv", index=False)

    # between-pull agreement
    if len(pulls) > 1:
        rel = []
        for (geo, term), s in terms.groupby(["geo", "term"]):
            cols = [f"r_pull{p}" for p in pulls]
            c = s[cols].corr().values[np.triu_indices(len(pulls), 1)]
            rel.append({"geo": geo, "term": term, "mean_between_pull_corr": round(float(np.nanmean(c)), 4),
                        "mean_sd_across_pulls": round(float(s.r_sd.mean()), 4)})
        pd.DataFrame(rel).to_csv(out / "trends_pull_reliability.csv", index=False)

    # daily national series (descriptive decay only)
    dd = []
    for term, req in TERM_SOURCE.items():
        if req == "R3":
            continue
        dd.append(daily[(daily.request == req) & (daily.term == term)])
    pd.concat(dd).drop(columns=["kind"]).to_csv(out / "trends_daily_GB.csv", index=False)

    # Wikipedia, global
    gl = json.loads((data / "raw" / "wikipedia" / "global_rest_pageviews.json").read_text())
    wrows = []
    for pid, a in gl["articles"].items():
        for it in a["response"]["items"]:
            wrows.append((int(pid), a["title"], pd.to_datetime(it["timestamp"][:8]), it["views"]))
    wg = pd.DataFrame(wrows, columns=["page_id", "title", "date", "views"])
    wg["week"] = wg.date - pd.to_timedelta((wg.date.dt.dayofweek + 1) % 7, unit="D")
    wgw = wg.groupby(["page_id", "title", "week"]).agg(views=("views", "sum"), n_days=("views", "size")).reset_index()
    wgw = wgw[wgw.n_days == 7]
    wgw["value"] = np.log(wgw.views)
    wgw["family"] = wgw.page_id.map(WIKI_FAMILIES)
    wgw["source"] = "global"

    # Wikipedia, UK-only (differentially private)
    uk = pd.read_csv(data / "raw" / "wikipedia" / "uk_dp_rows.tsv", sep="\t", parse_dates=["date"])
    man = pd.read_csv(data / "raw" / "wikipedia" / "uk_dp_manifest.csv", parse_dates=["date"])
    published = man[man.status == "ok"].date
    uk = uk.groupby(["page_id", "date"]).views.sum().reset_index()
    cov = uk.groupby("page_id").date.nunique() / len(published)
    ukcov = pd.DataFrame({"page_id": list(WIKI_FAMILIES)})
    ukcov["share_days_present"] = ukcov.page_id.map(cov).fillna(0).round(4)
    ukcov["passes"] = ukcov.share_days_present >= 1 - SPARSE_SHARE
    ukcov["family"] = ukcov.page_id.map(WIKI_FAMILIES)
    ukcov.to_csv(out / "wikipedia_uk_sparsity.csv", index=False)
    keep = ukcov[ukcov.passes].page_id
    uk = uk[uk.page_id.isin(keep)].copy()
    uk["week"] = uk.date - pd.to_timedelta((uk.date.dt.dayofweek + 1) % 7, unit="D")
    ukw = uk.groupby(["page_id", "week"]).agg(mean_views=("views", "mean"), n_days=("views", "size")).reset_index()
    ukw["value"] = np.log(ukw.mean_views)
    ukw["family"] = ukw.page_id.map(WIKI_FAMILIES)
    ukw["source"] = "uk_dp"
    titles = wg.drop_duplicates("page_id").set_index("page_id").title
    ukw["title"] = ukw.page_id.map(titles)
    pd.concat([wgw[["source", "page_id", "title", "family", "week", "value"]],
               ukw[["source", "page_id", "title", "family", "week", "value"]]]
              ).to_csv(out / "wikipedia_weekly.csv", index=False)
    print(f"pulls used: {pulls}")
    print(sp[~sp.passes].to_string(index=False) if (~sp.passes).any() else "all Trends series pass sparsity")
    print(ukcov.to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:3])
