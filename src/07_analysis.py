"""Pre-registered analyses (PREREGISTRATION.md, Sections 7-8; deviation 2 in docs/DEVIATIONS.md).

Analysis 1: national event study (UK, weekly), onset-aligned (1a) and end-aligned (1b), distributed-lag
  form with binned endpoints (k = -4 collects all earlier leads, k = 12 all later lags), reference k = -1,
  week-of-year and calendar-year fixed effects, Newey-West standard errors (8 lags).
Analysis 2: four-nation dose-response, lags k = 0..12 of weekly heat dose (per 10 degree-days),
  nation-by-year, nation-by-week-of-year and week fixed effects, Driscoll-Kraay standard errors (8 lags),
  randomisation inference by permuting whole nation dose series (24 permutations).
Confirmatory tests H1-H6 (Section 7.3), Holm-adjusted as one family, run twice:
  "plan" (sparsity rule as written) and "dev2" (deviation 2: no sparsity rule, floor safeguard).
Placebo and robustness checks: Section 8.

Usage:
  python src/07_analysis.py <processed_dir> <config_dir> <results_dir> [n_draws]
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml
from scipy import optimize, stats

SEED = 20261105
SAMPLE_START, SAMPLE_END = pd.Timestamp("2021-11-07"), pd.Timestamp("2026-10-25")  # first/last week starts
KS = [-4, -3, -2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]  # k = -1 omitted
NW_LAGS, DK_LAGS = 8, 8
EQ_BOUND = 0.10
NATIONS = ["England", "Scotland", "Wales", "Northern Ireland"]
NATION_GEO = {"England": "GB-ENG", "Scotland": "GB-SCT", "Wales": "GB-WLS", "Northern Ireland": "GB-NIR"}
WEEKS = pd.date_range(SAMPLE_START, SAMPLE_END, freq="7D")


def week_of(d):
    d = pd.DatetimeIndex(pd.to_datetime(d))
    return d - pd.to_timedelta((d.dayofweek + 1) % 7, unit="D")


def woy(weeks):
    return np.minimum((pd.DatetimeIndex(weeks).dayofyear - 1) // 7 + 1, 52)


# ---------------------------------------------------------------- Analysis 1 design
def event_matrix(event_weeks, weeks=WEEKS):
    """Columns b_k: number of events k weeks before week w; endpoints binned (Schmidheiny-Siegloch)."""
    ev = np.array(sorted(pd.DatetimeIndex(event_weeks)), dtype="datetime64[ns]")
    w = np.array(weeks, dtype="datetime64[ns]")
    rel = ((w[:, None] - ev[None, :]) / np.timedelta64(7, "D")).round().astype(int) if len(ev) else np.zeros((len(w), 0), int)
    cols = {}
    for k in KS:
        if k == -4:
            cols[k] = (rel <= -4).sum(1)
        elif k == 12:
            cols[k] = (rel >= 12).sum(1)
        else:
            cols[k] = (rel == k).sum(1)
    return pd.DataFrame(cols, index=weeks).rename(columns=lambda k: f"b{k}")


def fe_matrix(weeks=WEEKS):
    weeks = pd.DatetimeIndex(weeks)
    f = pd.DataFrame({"woy": woy(weeks), "year": weeks.year}, index=weeks)
    return pd.get_dummies(f.astype(str), drop_first=True, dtype=float).assign(const=1.0)


def fit1(y, X, F):
    d = pd.concat([y.rename("y"), X, F], axis=1).dropna()
    m = sm.OLS(d["y"], d.drop(columns="y")).fit(cov_type="HAC", cov_kwds={"maxlags": NW_LAGS})
    return m


def lincomb(m, names):
    w = np.zeros(len(m.params))
    for n in names:
        w[list(m.params.index).index(n)] = 1 / len(names)
    est = float(w @ m.params.values)
    se = float(np.sqrt(w @ m.cov_params().values @ w))
    return est, se


def fast_coef(y, X, F, names_sets):
    d = pd.concat([y.rename("y"), X, F], axis=1).dropna()
    A = d.drop(columns="y").values
    b = np.linalg.lstsq(A, d["y"].values, rcond=None)[0]
    idx = {c: i for i, c in enumerate(d.drop(columns="y").columns)}
    return [float(np.mean([b[idx[n]] for n in ns])) for ns in names_sets]


IMPACT, DECAY, RESIDUE = ["b0", "b1"], [f"b{k}" for k in range(4, 9)], [f"b{k}" for k in range(9, 13)]


def tost(est, se, bound=EQ_BOUND):
    p_lo = 1 - stats.norm.cdf((est + bound) / se)
    p_hi = stats.norm.cdf((est - bound) / se)
    return float(max(p_lo, p_hi))


def pz(est, se):
    return float(2 * (1 - stats.norm.cdf(abs(est / se))))


def half_life(betas, ks=np.arange(13)):
    f = lambda k, A, lam: A * np.exp(-lam * k)
    try:
        (A, lam), _ = optimize.curve_fit(f, ks, betas, p0=[betas[0] if betas[0] != 0 else 0.1, 0.5],
                                         bounds=([-np.inf, 1e-3], [np.inf, 5.0]), maxfev=20000)
    except Exception:
        return np.nan, np.nan, True
    at_bound = lam <= 1.001e-3 or lam >= 4.999
    return float(np.log(2) / lam), float(A), bool(at_bound)


def block_boot_resid(m, n, block=8, rng=None):
    e = m.resid.values
    T = len(e)
    starts = rng.integers(0, T - block + 1, size=int(np.ceil(T / block)))
    return np.concatenate([e[s:s + block] for s in starts])[:T]


# ---------------------------------------------------------------- Analysis 2 design
def weekly_dose(heat_daily, pct="p95"):
    d = heat_daily[heat_daily.percentile == pct].copy()
    d["week"] = week_of(d.date)
    g = d.groupby(["series", "week"]).agg(dose=("excess_degC", "sum"), n=("excess_degC", "size")).reset_index()
    return g[g.n == 7].pivot(index="week", columns="series", values="dose")


def residualise(M, F):
    q, _ = np.linalg.qr(F)
    return M - q @ (q.T @ M)


def panel_design(Y, dose, habit=False):
    """Y: weeks x nations outcome; dose: weeks x nations (degree-days). Returns stacked arrays."""
    rows = []
    full = dose.reindex(pd.date_range(dose.index.min(), SAMPLE_END, freq="7D")).fillna(0.0) / 10.0
    for n in NATIONS:
        for w in WEEKS:
            if pd.isna(Y.at[w, n]) if w in Y.index else True:
                continue
            lags = [full.at[w - pd.Timedelta(weeks=k), n] if (w - pd.Timedelta(weeks=k)) in full.index else np.nan
                    for k in range(13)]
            post = float(w.year >= 2023)
            extra = [lags[0] * post, lags[1] * post] if habit else []
            rows.append([n, w, Y.at[w, n]] + lags + extra)
    cols = ["nation", "week", "y"] + [f"d{k}" for k in range(13)] + (["d0_post", "d1_post"] if habit else [])
    d = pd.DataFrame(rows, columns=cols).dropna()
    F = pd.get_dummies(pd.DataFrame({
        "ny": d.nation + "_" + d.week.dt.year.astype(str),
        "nw": d.nation + "_" + pd.Series(woy(d.week), index=d.index).astype(str),
        "w": d.week.astype(str)}), drop_first=False, dtype=float).values
    return d, F


def fit2(d, F, xcols):
    X = residualise(d[xcols].values, F)
    y = residualise(d[["y"]].values, F)[:, 0]
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    # Driscoll-Kraay: HAC on cross-sectional sums of scores
    t_index = d.week.values
    weeks = np.unique(t_index)
    H = np.vstack([(X[t_index == t] * e[t_index == t, None]).sum(0) for t in weeks])
    S = H.T @ H
    for l in range(1, DK_LAGS + 1):
        G = H[l:].T @ H[:-l]
        S += (1 - l / (DK_LAGS + 1)) * (G + G.T)
    V = XtX_inv @ S @ XtX_inv
    return pd.Series(b, index=xcols), pd.DataFrame(V, index=xcols, columns=xcols)


def comb2(b, V, names, mean=False):
    w = pd.Series(0.0, index=b.index)
    for n in names:
        w[n] = 1 / len(names) if mean else 1.0
    return float(w @ b), float(np.sqrt(w @ V @ w))


def fast2(d, F, xcols):
    X = residualise(d[xcols].values, F)
    y = residualise(d[["y"]].values, F)[:, 0]
    return pd.Series(np.linalg.lstsq(X, y, rcond=None)[0], index=xcols)


# ---------------------------------------------------------------- helpers
def episodes(ep, series="UK", pct="p95"):
    e = ep[(ep.series == series) & (ep.percentile == pct)].copy()
    e["onset_week"] = pd.to_datetime(e.onset_week)
    e["end_week"] = pd.to_datetime(e.end_week)
    return e


def heat_weeks(heat_daily, series="UK", pct="p95"):
    d = heat_daily[(heat_daily.series == series) & (heat_daily.percentile == pct) & heat_daily.heatwave_day]
    return set(week_of(d.date))


def placebo_pool(hw_weeks):
    pool = []
    for w in WEEKS:
        if w.month in range(4, 10) and all(abs((w - h).days) > 28 for h in hw_weeks):
            pool.append(w)
    return pool


def holm(pvals):
    p = np.array(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    mult = np.empty(m)
    for i, j in enumerate(order):
        run = max(run, min(1.0, (m - i) * p[j]))
        adj[j] = run
        mult[j] = m - i
    return adj, mult


# ---------------------------------------------------------------- main
def main(proc, cfg, res, n_draws=2000):
    proc, cfg, res = Path(proc), Path(cfg), Path(res)
    res.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    heat_daily = pd.read_csv(proc / "heat_daily.csv", parse_dates=["date"])
    ep_all = pd.read_csv(proc / "heat_episodes.csv")
    fam = pd.read_csv(proc / "trends_weekly_families.csv", parse_dates=["date"])
    terms = pd.read_csv(proc / "trends_weekly_terms.csv", parse_dates=["date"])
    wiki = pd.read_csv(proc / "wikipedia_weekly.csv", parse_dates=["week"])
    anchor_av = pd.read_csv(proc / "anchor_validity_GB.csv", parse_dates=["date"]).set_index("date").value
    conc = yaml.safe_load(open(cfg / "concurrent_events.yaml"))
    mo = yaml.safe_load(open(cfg / "events_metoffice.yaml"))

    ep = episodes(ep_all)
    onsets, ends = list(ep.onset_week), list(ep.end_week)
    insample = ep[(ep.onset_week >= SAMPLE_START) & (ep.onset_week <= SAMPLE_END)]
    hw_weeks = heat_weeks(heat_daily)
    pool = placebo_pool(hw_weeks)
    F = fe_matrix()
    X1a, X1b = event_matrix(onsets), event_matrix(ends)

    # ---- anchor validity (Section 8.5) decides the primary anchor
    m_av = fit1(anchor_av.reindex(WEEKS), X1a, F)
    av_est, av_se = lincomb(m_av, IMPACT)
    anchor = "plumber" if pz(av_est, av_se) < 0.05 else "dentist"

    def fam_series(rule, scope, family, geo="GB", pull="mean", anc=None):
        s = fam[(fam.rule == rule) & (fam.scope == scope) & (fam.family == family) & (fam.geo == geo)
                & (fam.anchor == (anc or anchor)) & (fam.pull == pull)]
        return s.set_index("date").value.reindex(WEEKS) if len(s) else None

    def fam_floor(rule, family):
        s = fam[(fam.rule == rule) & (fam.scope == "national") & (fam.family == family) & (fam.geo == "GB")
                & (fam.anchor == anchor) & (fam.pull == "mean")]
        return s.set_index("date").floor_share.reindex(WEEKS) if len(s) else None

    def term_series(term, geo="GB"):
        col = "r_mean" if anchor == "dentist" else "r_mean_plumber_anchor"
        s = terms[(terms.term == term) & (terms.geo == geo)].set_index("date")[col]
        return s.reindex(WEEKS)

    # ---- outcomes for Analysis 1
    outcomes = {}
    for rule in ["plan", "dev2"]:
        for f_ in ["salience", "coping", "concern", "transition", "placebo"]:
            s = fam_series(rule, "national", f_)
            if s is not None:
                outcomes[f"trends_{rule}_{f_}"] = s
        c, k = fam_series(rule, "national", "coping"), fam_series(rule, "national", "concern")
        if c is not None and k is not None:
            outcomes[f"trends_{rule}_coping_minus_concern"] = c - k
    for t in ["heatwave", "air conditioning", "portable air conditioner", "heat stroke", "climate change",
              "global warming", "net zero", "heat pump", "car insurance", "recipes"]:
        outcomes[f"term_{t}"] = term_series(t)
    for src in ["global", "uk_dp"]:
        w = wiki[wiki.source == src]
        for f_, g in w.groupby("family"):
            outcomes[f"wiki_{src}_{f_}"] = g.groupby("week").value.mean().reindex(WEEKS)
        for (pid, title), g in w.groupby(["page_id", "title"]):
            outcomes[f"wikiart_{src}_{title}"] = g.set_index("week").value.reindex(WEEKS)
    outcomes["anchor_dentist_minus_plumber"] = anchor_av.reindex(WEEKS)
    step = pd.DataFrame({"s1": (WEEKS >= pd.Timestamp("2023-02-05")).astype(float),
                         "s2": (WEEKS >= pd.Timestamp("2023-05-21")).astype(float)}, index=WEEKS)

    def F_for(name):
        return pd.concat([F, step], axis=1) if name.startswith("wiki") and "uk_dp" in name else F

    # ---- Analysis 1: main estimates
    a1, coefs = [], []
    for name, y in outcomes.items():
        if y is None or y.notna().sum() < 100:
            continue
        Fn = F_for(name)
        m_a, m_b = fit1(y, X1a, Fn), fit1(y, X1b, Fn)
        imp, imp_se = lincomb(m_a, IMPACT)
        dec, dec_se = lincomb(m_b, DECAY)
        resd, res_se = lincomb(m_b, RESIDUE)
        wald = m_a.wald_test("b-4 = 0, b-3 = 0, b-2 = 0", scalar=True)
        betas = np.array([m_b.params[f"b{k}"] for k in range(13)])
        hl, A, at_b = half_life(betas)
        a1.append({"outcome": name, "n_weeks": int(m_a.nobs), "impact": imp, "impact_se": imp_se,
                   "impact_p": pz(imp, imp_se), "decay": dec, "decay_se": dec_se, "decay_p": pz(dec, dec_se),
                   "decay_tost_p": tost(dec, dec_se), "residue": resd, "residue_se": res_se,
                   "residue_p": pz(resd, res_se), "residue_tost_p": tost(resd, res_se),
                   "leads_wald_p": float(wald.pvalue), "half_life_weeks": hl, "half_life_at_bound": at_b})
        for lab, m in [("1a", m_a), ("1b", m_b)]:
            for k in KS:
                coefs.append({"outcome": name, "model": lab, "k": k, "beta": m.params[f"b{k}"],
                              "se": m.bse[f"b{k}"]})
            coefs.append({"outcome": name, "model": lab, "k": -1, "beta": 0.0, "se": 0.0})
    a1 = pd.DataFrame(a1).set_index("outcome")
    pd.DataFrame(coefs).to_csv(res / "analysis1_event_coefficients.csv", index=False)

    # ---- placebo dates RI + MDE + half-life bootstrap for key outcomes
    key = [n for n in a1.index if n.startswith("trends_") or n.startswith("term_")
           or n.startswith("wiki_") or n == "anchor_dentist_minus_plumber"]
    n_on = len(insample)
    draws = [sorted(rng.choice(len(pool), size=n_on, replace=False)) for _ in range(n_draws)]
    ri = {n: {"impact": [], "decay": [], "residue": []} for n in key}
    band_out = [n for n in key if n.startswith("trends_") or n == "term_heatwave"]
    band = {n: [] for n in band_out}
    allk = [[f"b{k}"] for k in KS]
    for dr in draws:
        fake = [pool[i] for i in dr]
        Xf = event_matrix(fake)
        for n in key:
            i_, d_, r_ = fast_coef(outcomes[n], Xf, F_for(n), [IMPACT, DECAY, RESIDUE])
            ri[n]["impact"].append(i_); ri[n]["decay"].append(d_); ri[n]["residue"].append(r_)
            if n in band:
                band[n].append(fast_coef(outcomes[n], Xf, F_for(n), allk))
    bq = []
    for n, arr in band.items():
        arr = np.array(arr)
        for j, k in enumerate(KS):
            bq.append({"outcome": n, "k": k, "q025": np.percentile(arr[:, j], 2.5),
                       "q975": np.percentile(arr[:, j], 97.5), "sd": arr[:, j].std(ddof=1)})
    pd.DataFrame(bq).to_csv(res / "placebo_beta_bands.csv", index=False)
    for n in key:
        for est in ["impact", "decay", "residue"]:
            arr = np.array(ri[n][est])
            a1.loc[n, f"{est}_ri_p"] = float((1 + np.sum(np.abs(arr) >= abs(a1.loc[n, est]))) / (1 + len(arr)))
            a1.loc[n, f"{est}_mde"] = float(2.8 * arr.std(ddof=1))

    hl_rows = []
    for n in [k for k in key if k.startswith("trends_") or k in ("term_heatwave",)]:
        y = outcomes[n]
        m_b = fit1(y, X1b, F)
        d = pd.concat([y.rename("y"), X1b, F], axis=1).dropna()
        A_ = d.drop(columns="y")
        fitted = m_b.fittedvalues.values
        hls = []
        for _ in range(n_draws):
            ystar = fitted + block_boot_resid(m_b, len(fitted), 8, rng)
            b = np.linalg.lstsq(A_.values, ystar, rcond=None)[0]
            bb = np.array([b[list(A_.columns).index(f"b{k}")] for k in range(13)])
            hls.append(half_life(bb)[0])
        hls = np.array(hls)
        lo, hi = np.nanpercentile(hls, [2.5, 97.5])
        hl, _, at_b = a1.loc[n, "half_life_weeks"], None, a1.loc[n, "half_life_at_bound"]
        hl_rows.append({"outcome": n, "half_life_weeks": hl, "ci_low": lo, "ci_high": hi,
                        "identifiable": (not at_b) and hi <= 12})
    pd.DataFrame(hl_rows).to_csv(res / "analysis1_half_life.csv", index=False)

    # ---- floor shares (deviation 2 safeguard)
    floor_rows = []
    for rule in ["plan", "dev2"]:
        for f_ in ["salience", "coping", "concern"]:
            fl = fam_floor(rule, f_)
            if fl is None:
                continue
            for lab, ks in [("decay", range(4, 9)), ("residue", range(9, 13))]:
                wk = [e + pd.Timedelta(weeks=k) for e in insample.end_week for k in ks]
                wk = [w for w in wk if w in fl.index]
                floor_rows.append({"rule": rule, "family": f_, "window": lab,
                                   "floor_share": float(fl.reindex(wk).mean()) if wk else np.nan})
    floor = pd.DataFrame(floor_rows)
    floor.to_csv(res / "floor_shares.csv", index=False)

    a1.reset_index().to_csv(res / "analysis1_estimates.csv", index=False)

    # ---- Analysis 2
    dose = weekly_dose(heat_daily)
    a2_rows, perms = [], list(itertools.permutations(range(4)))
    a2_store = {}
    for rule in ["plan", "dev2"]:
        for f_ in ["coping", "concern", "salience", "placebo", "transition"]:
            Y = pd.DataFrame({n: fam_series(rule, "four_nation", f_, NATION_GEO[n]) for n in NATIONS}) \
                if fam_series(rule, "four_nation", f_, "GB-ENG") is not None else None
            if Y is None:
                continue
            for habit in [False, True]:
                d, Fp = panel_design(Y, dose, habit)
                xc = [f"d{k}" for k in range(13)] + (["d0_post", "d1_post"] if habit else [])
                b, V = fit2(d, Fp, xc)
                row = {"rule": rule, "family": f_, "habituation_model": habit, "n_obs": len(d)}
                sr, sr_se = comb2(b, V, ["d0", "d1"])
                rs, rs_se = comb2(b, V, [f"d{k}" for k in range(9, 13)], mean=True)
                row.update({"short_run": sr, "short_run_se": sr_se, "short_run_p": pz(sr, sr_se),
                            "residue": rs, "residue_se": rs_se, "residue_p": pz(rs, rs_se)})
                if habit:
                    hb, hb_se = comb2(b, V, ["d0_post", "d1_post"])
                    row.update({"habituation": hb, "habituation_se": hb_se, "habituation_p": pz(hb, hb_se)})
                pd.DataFrame({"rule": rule, "family": f_, "habituation_model": habit, "term": b.index,
                              "beta": b.values, "se": np.sqrt(np.diag(V.values))}).to_csv(
                    res / f"analysis2_lags_{rule}_{f_}_{'habit' if habit else 'main'}.csv", index=False)
                # permutation RI: permute whole nation dose series across nations
                pr_sr, pr_hb, pr_lags = [], [], []
                for p in perms:
                    dz = dose.copy()
                    dz[NATIONS] = dose[[NATIONS[i] for i in p]].values
                    dp, Fq = panel_design(Y, dz, habit)
                    bp = fast2(dp, Fq, xc)
                    pr_lags.append(bp[[f"d{k}" for k in range(13)]].values)
                    pr_sr.append(bp["d0"] + bp["d1"])
                    if habit:
                        pr_hb.append(bp["d0_post"] + bp["d1_post"])
                pl = np.array(pr_lags)
                pd.DataFrame({"k": range(13), "perm_min": pl.min(0), "perm_max": pl.max(0)}).to_csv(
                    res / f"analysis2_permbands_{rule}_{f_}_{'habit' if habit else 'main'}.csv", index=False)
                row["short_run_ri_p"] = float(np.mean(np.abs(pr_sr) >= abs(sr) - 1e-12))
                row["short_run_mde"] = float(2.8 * np.std(pr_sr, ddof=1))
                if habit:
                    row["habituation_ri_p"] = float(np.mean(np.abs(pr_hb) >= abs(row["habituation"]) - 1e-12))
                    row["habituation_mde"] = float(2.8 * np.std(pr_hb, ddof=1))
                # year-shift placebo
                for shift in [-2, -1, 1, 2]:
                    dz = dose.copy()
                    dz.index = dz.index + pd.Timedelta(weeks=52 * shift)
                    dp, Fq = panel_design(Y, dz.reindex(dose.index), habit)
                    if len(dp) > 200:
                        bp = fast2(dp, Fq, xc)
                        row[f"short_run_yearshift_{shift:+d}"] = float(bp["d0"] + bp["d1"])
                a2_rows.append(row)
                a2_store[(rule, f_, habit)] = row
    a2 = pd.DataFrame(a2_rows)
    a2.to_csv(res / "analysis2_estimates.csv", index=False)

    # ---- Section 8 robustness (Analysis 1, confirmatory national outcomes)
    rob_out = [n for n in a1.index if n.startswith("trends_") or n in ("term_heatwave", "term_car insurance",
                                                                          "term_recipes")]
    rob = []
    # placebo years
    for shift in [-2, -1, 1, 2]:
        sh = []
        for o in insample.onset_week:
            s = week_of([o + pd.DateOffset(years=shift)])[0]
            if SAMPLE_START <= s <= SAMPLE_END and all(abs((s - h).days) > 28 for h in hw_weeks):
                sh.append(s)
        if not sh:
            continue
        Xs = event_matrix(sh)
        for n in rob_out:
            m = fit1(outcomes[n], Xs, F)
            e_, s_ = lincomb(m, IMPACT)
            rob.append({"check": f"placebo_year_{shift:+d}", "outcome": n, "estimand": "impact",
                        "estimate": e_, "se": s_, "p": pz(e_, s_), "n_events": len(sh)})
    # concurrent events dropped
    drop = set()
    for ev in conc["cop"] + conc["other"]:
        for w in WEEKS:
            if w <= pd.Timestamp(ev["window_end"]) and w + pd.Timedelta(days=6) >= pd.Timestamp(ev["window_start"]):
                drop.add(w)
    keep = [w for w in WEEKS if w not in drop]
    for n in rob_out:
        y = outcomes[n].reindex(keep)
        ma, mb = fit1(y, X1a.reindex(keep), F.reindex(keep)), fit1(y, X1b.reindex(keep), F.reindex(keep))
        for lab, m, ks in [("impact", ma, IMPACT), ("decay", mb, DECAY), ("residue", mb, RESIDUE)]:
            e_, s_ = lincomb(m, ks)
            rob.append({"check": "drop_concurrent_events", "outcome": n, "estimand": lab,
                        "estimate": e_, "se": s_, "p": pz(e_, s_), "n_events": len(drop)})
    # thresholds and Met Office dating
    alt = {}
    for pct in ["p90", "p975"]:
        e2 = episodes(ep_all, pct=pct)
        alt[f"threshold_{pct}"] = (list(e2.onset_week), list(e2.end_week))
    mo_on = [week_of([e["start"]])[0] for e in mo["events"]]
    mo_end = [week_of([e["end"]])[0] for e in mo["events"]]
    alt["metoffice_dates"] = (mo_on, mo_end)
    for lab, (on, en) in alt.items():
        Xa, Xb = event_matrix(on), event_matrix(en)
        for n in rob_out:
            ma, mb = fit1(outcomes[n], Xa, F), fit1(outcomes[n], Xb, F)
            for est, m, ks in [("impact", ma, IMPACT), ("decay", mb, DECAY), ("residue", mb, RESIDUE)]:
                e_, s_ = lincomb(m, ks)
                rob.append({"check": lab, "outcome": n, "estimand": est, "estimate": e_, "se": s_,
                            "p": pz(e_, s_), "n_events": len(on)})
    # pull noise: single pulls
    pulls = sorted(p for p in fam.pull.unique() if p != "mean")
    for p in pulls:
        for rule in ["plan", "dev2"]:
            for f_ in ["salience", "coping", "concern"]:
                s = fam_series(rule, "national", f_, pull=p)
                if s is None:
                    continue
                ma, mb = fit1(s, X1a, F), fit1(s, X1b, F)
                for est, m, ks in [("impact", ma, IMPACT), ("residue", mb, RESIDUE)]:
                    e_, s_ = lincomb(m, ks)
                    rob.append({"check": f"single_{p}", "outcome": f"trends_{rule}_{f_}", "estimand": est,
                                "estimate": e_, "se": s_, "p": pz(e_, s_), "n_events": len(onsets)})
    rob = pd.DataFrame(rob)
    rob.to_csv(res / "robustness.csv", index=False)

    # ---- Section 7.4: per-episode impacts (descriptive)
    epi_rows = []
    for rule in ["plan", "dev2"]:
        for f_ in ["salience", "coping", "concern"]:
            y = fam_series(rule, "national", f_)
            if y is None:
                continue
            X = X1a.drop(columns=["b0", "b1"]).copy()
            for i, r in insample.reset_index(drop=True).iterrows():
                X[f"e{i}_0"] = (WEEKS == r.onset_week).astype(float)
                X[f"e{i}_1"] = (WEEKS == r.onset_week + pd.Timedelta(weeks=1)).astype(float)
            m = fit1(y, X, F)
            for i, r in insample.reset_index(drop=True).iterrows():
                e_, s_ = lincomb(m, [f"e{i}_0", f"e{i}_1"])
                epi_rows.append({"rule": rule, "family": f_, "episode_start": r.start, "peak_tmax": r.peak_tmax_degC,
                                 "dose": r.dose_degree_days, "impact": e_, "se": s_})
    epi = pd.DataFrame(epi_rows)
    epi.to_csv(res / "per_episode_impacts.csv", index=False)

    # ---- confirmatory hypotheses
    def a1g(n, col):
        return a1.loc[n, col] if n in a1.index else np.nan

    def a2g(rule, f_, habit, col):
        r = a2_store.get((rule, f_, habit))
        return r.get(col, np.nan) if r else np.nan

    fl = floor.set_index(["rule", "family", "window"]).floor_share if len(floor) else pd.Series(dtype=float)
    hyp_all = []
    for rule in ["plan", "dev2"]:
        H = [
            ("H1", "Heat raises salience", "trends heatwave (UK)", "impact", "positive",
             a1g(f"trends_{rule}_salience", "impact"), a1g(f"trends_{rule}_salience", "impact_se"),
             a1g(f"trends_{rule}_salience", "impact_p"), a1g(f"trends_{rule}_salience", "impact_ri_p"),
             a1g(f"trends_{rule}_salience", "impact_mde"), None),
            ("H2", "Coping attention decays within about 2 months", "coping index (UK)", "decay", "equivalent",
             a1g(f"trends_{rule}_coping", "decay"), a1g(f"trends_{rule}_coping", "decay_se"),
             a1g(f"trends_{rule}_coping", "decay_tost_p"), a1g(f"trends_{rule}_coping", "decay_ri_p"),
             a1g(f"trends_{rule}_coping", "decay_mde"), fl.get((rule, "coping", "decay"), np.nan)),
            ("H3", "Concern attention leaves no residue", "concern index (UK)", "residue", "equivalent",
             a1g(f"trends_{rule}_concern", "residue"), a1g(f"trends_{rule}_concern", "residue_se"),
             a1g(f"trends_{rule}_concern", "residue_tost_p"), a1g(f"trends_{rule}_concern", "residue_ri_p"),
             a1g(f"trends_{rule}_concern", "residue_mde"), fl.get((rule, "concern", "residue"), np.nan)),
            ("H4", "Coping responds more than concern", "coping minus concern (UK)", "impact", "positive",
             a1g(f"trends_{rule}_coping_minus_concern", "impact"), a1g(f"trends_{rule}_coping_minus_concern", "impact_se"),
             a1g(f"trends_{rule}_coping_minus_concern", "impact_p"),
             a1g(f"trends_{rule}_coping_minus_concern", "impact_ri_p"),
             a1g(f"trends_{rule}_coping_minus_concern", "impact_mde"), None),
            ("H5a", "Hotter nations respond more (coping)", "coping index (4 nations)", "short-run", "positive",
             a2g(rule, "coping", False, "short_run"), a2g(rule, "coping", False, "short_run_se"),
             a2g(rule, "coping", False, "short_run_p"), a2g(rule, "coping", False, "short_run_ri_p"),
             a2g(rule, "coping", False, "short_run_mde"), None),
            ("H5b", "Hotter nations respond more (concern)", "concern index (4 nations)", "short-run", "two-sided",
             a2g(rule, "concern", False, "short_run"), a2g(rule, "concern", False, "short_run_se"),
             a2g(rule, "concern", False, "short_run_p"), a2g(rule, "concern", False, "short_run_ri_p"),
             a2g(rule, "concern", False, "short_run_mde"), None),
            ("H6", "Habituation after 2022", "coping index (4 nations)", "interaction", "two-sided",
             a2g(rule, "coping", True, "habituation"), a2g(rule, "coping", True, "habituation_se"),
             a2g(rule, "coping", True, "habituation_p"), a2g(rule, "coping", True, "habituation_ri_p"),
             a2g(rule, "coping", True, "habituation_mde"), None),
        ]
        h = pd.DataFrame(H, columns=["id", "hypothesis", "outcome", "estimand", "prediction", "estimate", "se",
                                     "p_model", "p_ri", "mde", "floor_share"])
        h.insert(0, "rule", rule)
        h["testable"] = h.estimate.notna()
        t = h[h.testable].copy()
        adj, mult = holm(t.p_model.values)
        h.loc[t.index, "p_holm"] = adj
        h.loc[t.index, "holm_multiplier"] = mult
        labels = []
        for _, r in h.iterrows():
            if not r.testable:
                labels.append("not testable (no series passes the sparsity rule)")
                continue
            if r.prediction == "equivalent":
                p_two = min(1.0, pz(r.estimate, r.se) * r.holm_multiplier)
                if r.floor_share is not None and not pd.isna(r.floor_share) and r.floor_share > 0.5 and rule == "dev2":
                    labels.append("floor-limited")
                elif r.p_holm < 0.05:
                    labels.append("reset (equivalent to zero)")
                elif p_two < 0.05 and r.p_ri < 0.05:
                    labels.append("residue" if r.estimate > 0 else "undershoot")
                else:
                    labels.append("inconclusive")
            else:
                present = r.p_holm < 0.05 and r.p_ri < 0.05
                if not present:
                    labels.append("no effect detected")
                elif r.prediction == "positive":
                    labels.append("effect present (positive)" if r.estimate > 0 else "effect present (opposite sign)")
                else:
                    labels.append("effect present (" + ("positive" if r.estimate > 0 else "negative") + ")")
        h["result"] = labels
        # calibrated version (docs/DEVIATIONS.md, entry 3)
        h["se_calibrated"] = np.maximum(h.se, h.mde / 2.8)
        h["p_calibrated"] = [tost(r.estimate, r.se_calibrated) if r.prediction == "equivalent" else r.p_ri
                             for _, r in h.iterrows()]
        if len(t):
            adj_c, mult_c = holm(h.loc[t.index, "p_calibrated"].values)
            h.loc[t.index, "p_calibrated_holm"] = adj_c
            h.loc[t.index, "holm_multiplier_calibrated"] = mult_c
        labels_c = []
        for _, r in h.iterrows():
            if not r.testable:
                labels_c.append("not testable (no series passes the sparsity rule)")
            elif r.prediction == "equivalent":
                if rule == "dev2" and not pd.isna(r.floor_share) and r.floor_share > 0.5:
                    labels_c.append("floor-limited")
                elif r.p_calibrated_holm < 0.05:
                    labels_c.append("reset (equivalent to zero)")
                elif min(1.0, r.p_ri * r.holm_multiplier_calibrated) < 0.05:
                    labels_c.append("residue" if r.estimate > 0 else "undershoot")
                else:
                    labels_c.append("inconclusive")
            elif r.p_calibrated_holm < 0.05:
                if r.prediction == "positive":
                    labels_c.append("effect present (positive)" if r.estimate > 0 else "effect present (opposite sign)")
                else:
                    labels_c.append("effect present (" + ("positive" if r.estimate > 0 else "negative") + ")")
            else:
                labels_c.append("no effect detected")
        h["result_calibrated"] = labels_c
        hyp_all.append(h)
    hyp = pd.concat(hyp_all)
    hyp.to_csv(res / "hypotheses.csv", index=False)

    meta = {"anchor_used": anchor, "anchor_validity_impact": av_est, "anchor_validity_p": pz(av_est, av_se),
            "n_insample_onsets": n_on, "insample_onsets": [str(x.date()) for x in insample.onset_week],
            "placebo_pool_weeks": len(pool), "n_draws": n_draws, "seed": SEED,
            "pulls": pulls, "dropped_concurrent_weeks": len(drop)}
    (res / "run_info.json").write_text(json.dumps(meta, indent=1))
    pd.set_option("display.width", 220)
    print(json.dumps(meta, indent=1))
    print(hyp[["rule", "id", "estimate", "se", "p_model", "p_holm", "p_ri", "mde", "floor_share", "result",
               "p_calibrated_holm", "result_calibrated"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 2000)
