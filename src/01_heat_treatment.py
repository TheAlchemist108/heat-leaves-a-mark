"""Build the heat treatment from HadUK-Grid daily maximum temperature (PREREGISTRATION.md, Section 4).

For each series (England, Scotland, Wales, Northern Ireland, and a population-weighted UK):
  - threshold T* = percentile of daily Tmax over June-August 1991-2020 (95th primary; 90th and 97.5th robustness)
  - heatwave day = a day in a run of >= 3 consecutive days with Tmax >= T*
  - spells separated by fewer than 14 days (i.e. a gap of 13 or fewer non-heatwave days) merge into one episode
  - onset week / end week = Sunday-Saturday week containing the first / last heatwave day of an episode
  - weekly heat dose = sum over the week's days of max(0, Tmax - T*), in degree-days

Usage:
  python src/01_heat_treatment.py <haduk_country_tasmax.nc> <population_weights.yaml> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

NATIONS = {"England": "England", "Scotland": "Scotland", "Wales": "Wales", "Northern Ireland": "Northern Ireland"}
PERCENTILES = {"p95": 95.0, "p90": 90.0, "p975": 97.5}
BASE_START, BASE_END = "1991-01-01", "2020-12-31"
MIN_RUN = 3
MERGE_GAP = 13  # max non-heatwave days between spells that still merge (fewer than 14)
SAMPLE_START = "2021-11-07"  # first Sunday of the primary sample (Section 3)


def load_series(nc_path, weights):
    ds = xr.open_dataset(nc_path)
    names = [g.decode().strip() for g in ds.geo_region.values]
    df = ds.tasmax.to_pandas()
    df.columns = names
    df.index = pd.to_datetime(df.index).normalize()
    df = df[list(NATIONS)]
    w = pd.Series({n: weights[n] for n in NATIONS}, dtype=float)
    df["UK"] = (df[list(NATIONS)] * w).sum(axis=1, min_count=4) / w.sum()
    return df


def week_start(dates):
    # Sunday that starts the Sunday-Saturday week containing each date
    dates = pd.DatetimeIndex(dates)
    return dates - pd.to_timedelta((dates.dayofweek + 1) % 7, unit="D")


def heatwave_days(tmax, thr):
    hot = (tmax >= thr).to_numpy()
    out = np.zeros(len(hot), dtype=bool)
    i = 0
    while i < len(hot):
        if hot[i]:
            j = i
            while j < len(hot) and hot[j]:
                j += 1
            if j - i >= MIN_RUN:
                out[i:j] = True
            i = j
        else:
            i += 1
    return pd.Series(out, index=tmax.index)


def episodes(hw_days):
    days = hw_days[hw_days].index
    if len(days) == 0:
        return []
    eps, start, prev = [], days[0], days[0]
    for d in days[1:]:
        if (d - prev).days - 1 <= MERGE_GAP:
            prev = d
        else:
            eps.append((start, prev))
            start, prev = d, d
    eps.append((start, prev))
    return eps


def main(nc_path, weights_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    weights = yaml.safe_load(open(weights_path))["population"]
    tmax = load_series(nc_path, weights)

    if tmax.isnull().any().any():
        missing = tmax.isnull().sum()
        print("WARNING: missing temperature values (not imputed):\n", missing[missing > 0])

    base = tmax.loc[BASE_START:BASE_END]
    base = base[base.index.month.isin([6, 7, 8])]

    thr_rows, ep_rows, daily_rows, weekly = [], [], [], []
    for series in tmax.columns:
        s = tmax[series]
        for label, p in PERCENTILES.items():
            thr = float(np.nanpercentile(base[series], p))
            thr_rows.append({"series": series, "percentile": label, "threshold_degC": round(thr, 3)})
            hw = heatwave_days(s, thr)
            for k, (a, b) in enumerate(episodes(hw)):
                if b < pd.Timestamp(SAMPLE_START) - pd.Timedelta(days=120):
                    continue
                ep_days = s.loc[a:b]
                ep_rows.append({
                    "series": series, "percentile": label,
                    "in_sample": a >= pd.Timestamp(SAMPLE_START),
                    "start": a.date(), "end": b.date(),
                    "length_days": (b - a).days + 1,
                    "heatwave_days": int(hw.loc[a:b].sum()),
                    "onset_week": week_start([a])[0].date(),
                    "end_week": week_start([b])[0].date(),
                    "peak_tmax_degC": round(float(ep_days.max()), 2),
                    "dose_degree_days": round(float((ep_days - thr).clip(lower=0).sum()), 2),
                })
            excess = (s - thr).clip(lower=0)
            d = pd.DataFrame({"date": s.index, "series": series, "percentile": label,
                              "tmax_degC": s.values, "threshold_degC": thr,
                              "heatwave_day": hw.values, "excess_degC": excess.values})
            daily_rows.append(d[d.date >= pd.Timestamp(SAMPLE_START) - pd.Timedelta(days=120)])
            wk = pd.DataFrame({"week": week_start(s.index), "dose": excess.values, "hw": hw.values})
            wk = wk.groupby("week").agg(dose_degree_days=("dose", "sum"), heatwave_days=("hw", "sum"),
                                        n_days=("dose", "size")).reset_index()
            wk["series"], wk["percentile"] = series, label
            weekly.append(wk[wk.week >= pd.Timestamp(SAMPLE_START)])

    pd.DataFrame(thr_rows).to_csv(out / "heat_thresholds.csv", index=False)
    pd.DataFrame(ep_rows).to_csv(out / "heat_episodes.csv", index=False)
    pd.concat(daily_rows).to_csv(out / "heat_daily.csv", index=False)
    weekly = pd.concat(weekly)
    weekly = weekly[weekly.n_days == 7]  # drop the incomplete final week
    weekly["dose_degree_days"] = weekly["dose_degree_days"].round(3)
    weekly.to_csv(out / "heat_weekly.csv", index=False)
    print(f"Wrote thresholds, {len(ep_rows)} episode rows and weekly dose to {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
