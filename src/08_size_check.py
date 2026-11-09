"""Size check of the pre-registered model-based tests on simulated white noise (docs/DEVIATIONS.md, entry 3).

Uses the real weeks, fixed effects and heatwave timing, with outcomes replaced by independent N(0, 1)
noise. A correctly sized test rejects 5% of the time at the 5% level.

Usage:
  python src/08_size_check.py <processed_dir> <results_dir> [n_sims]
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("analysis", Path(__file__).with_name("07_analysis.py"))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def main(proc, res, n=500):
    proc, res = Path(proc), Path(res)
    res.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.SEED)
    ep = a.episodes(pd.read_csv(proc / "heat_episodes.csv"))
    F, X1a, X1b = a.fe_matrix(), a.event_matrix(ep.onset_week), a.event_matrix(ep.end_week)
    dose = a.weekly_dose(pd.read_csv(proc / "heat_daily.csv", parse_dates=["date"]))
    z = {"A1 impact (Newey-West)": [], "A1 decay (Newey-West)": [], "A1 residue (Newey-West)": [],
         "A2 short-run (Driscoll-Kraay)": []}
    for _ in range(n):
        y = pd.Series(rng.normal(size=len(a.WEEKS)), index=a.WEEKS)
        ma, mb = a.fit1(y, X1a, F), a.fit1(y, X1b, F)
        for lab, m, ks in [("impact", ma, a.IMPACT), ("decay", mb, a.DECAY), ("residue", mb, a.RESIDUE)]:
            e, s = a.lincomb(m, ks)
            z[f"A1 {lab} (Newey-West)"].append(e / s)
        Y = pd.DataFrame(rng.normal(size=(len(a.WEEKS), 4)), index=a.WEEKS, columns=a.NATIONS)
        d, Fp = a.panel_design(Y, dose, False)
        b, V = a.fit2(d, Fp, [f"d{k}" for k in range(13)])
        e, s = a.comb2(b, V, ["d0", "d1"])
        z["A2 short-run (Driscoll-Kraay)"].append(e / s)
    out = pd.DataFrame([{"test": k, "n_sims": n, "rejection_rate_5pct": float(np.mean(np.abs(v) > 1.96)),
                         "sd_of_z": float(np.std(v))} for k, v in z.items()])
    out.to_csv(res / "size_check.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 500)
