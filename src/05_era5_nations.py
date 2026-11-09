"""ERA5 nation averages and calibration to HadUK-Grid (PREREGISTRATION.md, Section 4.1).

1. Average ERA5 daily maximum 2 m temperature over each nation's land area, using the
   Natural Earth 1:10m admin-0 map units (England, Scotland, Wales, Northern Ireland).
   Each grid cell is weighted by the fraction of it inside the nation and by cos(latitude).
2. Calibration, for each nation, on April-October 2016-2025 daily values, as written in the plan:
   ERA5 = a + b * HadUK. Calibrated 2026 values are (ERA5 - a) / b, i.e. on the HadUK scale.
   The direct regression HadUK = c + d * ERA5 is also reported as a check.
3. Rule: if the daily correlation is below 0.95 for a nation, 2026 is excluded for that nation.
   The UK series is the population-weighted mean of the four calibrated nation series; its
   correlation with the HadUK-based UK series is checked against the same 0.95 rule.

Usage:
  python src/05_era5_nations.py <era5_dir> <naturalearth_map_units.zip> <haduk_country_tasmax.nc>
                                <population_weights.yaml> <out_dir>
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import regionmask
import xarray as xr
import yaml

NATIONS = {"ENG": "England", "SCT": "Scotland", "WLS": "Wales", "NIR": "Northern Ireland"}
CAL_YEARS = range(2016, 2026)
MIN_R = 0.95


def era5_nations(era5_dir, ne_zip):
    files = sorted(p for p in Path(era5_dir).glob("era5_tmax_uk_*.nc"))
    t = xr.concat([xr.open_dataset(f)["t2m"].load() for f in files], dim="valid_time").sortby("valid_time")
    t = (t - 273.15).rename({"valid_time": "time"})
    g = gpd.read_file(f"zip://{ne_zip}")
    g = g[g["GU_A3"].isin(NATIONS)].reset_index(drop=True)
    regions = regionmask.from_geopandas(g, names="GU_A3", abbrevs="GU_A3", overlap=False)
    frac = regions.mask_3D_frac_approx(t.longitude, t.latitude)
    w = frac * np.cos(np.deg2rad(t.latitude))
    out = {}
    for code, name in NATIONS.items():
        wi = w.sel(region=frac.abbrevs == code).squeeze("region", drop=True)
        out[name] = (t * wi).sum(("latitude", "longitude")) / wi.sum()
    df = pd.DataFrame({k: v.to_series() for k, v in out.items()})
    df.index = pd.to_datetime(df.index).normalize()
    return df, files


def haduk(nc_path):
    ds = xr.open_dataset(nc_path)
    names = [x.decode().strip() for x in ds.geo_region.values]
    df = ds.tasmax.to_pandas()
    df.columns = names
    df.index = pd.to_datetime(df.index).normalize()
    return df[list(NATIONS.values())]


def main(era5_dir, ne_zip, haduk_nc, weights_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    e5, files = era5_nations(era5_dir, ne_zip)
    hd = haduk(haduk_nc)
    wts = yaml.safe_load(open(weights_path))["population"]
    w = pd.Series({n: wts[n] for n in NATIONS.values()}, dtype=float)

    rows, cal = [], pd.DataFrame(index=e5.index)
    cal_days = e5.index[(e5.index.year.isin(CAL_YEARS)) & (e5.index.month.isin(range(4, 11)))]
    for name in NATIONS.values():
        x = hd.loc[cal_days, name]
        y = e5.loc[cal_days, name]
        ok = x.notna() & y.notna()
        b, a = np.polyfit(x[ok], y[ok], 1)       # ERA5 = a + b * HadUK (as in the plan)
        d, c = np.polyfit(y[ok], x[ok], 1)       # HadUK = c + d * ERA5 (check)
        r = float(np.corrcoef(x[ok], y[ok])[0, 1])
        cal[name] = (e5[name] - a) / b
        resid = cal.loc[cal_days, name][ok] - x[ok]
        rows.append({"series": name, "n_days": int(ok.sum()), "r": round(r, 4),
                     "passes_rule": r >= MIN_R, "a": round(a, 4), "b": round(b, 4),
                     "check_c": round(c, 4), "check_d": round(d, 4),
                     "rmse_after_calibration_degC": round(float(np.sqrt((resid ** 2).mean())), 3)})

    uk_cal = (cal[list(NATIONS.values())] * w).sum(axis=1, min_count=4) / w.sum()
    uk_hd = (hd[list(NATIONS.values())] * w).sum(axis=1, min_count=4) / w.sum()
    r_uk = float(np.corrcoef(uk_hd.loc[cal_days], uk_cal.loc[cal_days])[0, 1])
    rows.append({"series": "UK", "n_days": len(cal_days), "r": round(r_uk, 4), "passes_rule": r_uk >= MIN_R})
    pd.DataFrame(rows).to_csv(out / "era5_calibration.csv", index=False)

    ext = cal.loc["2026-01-01":"2026-12-31"].copy()
    for row in rows[:-1]:
        if not row["passes_rule"]:
            ext[row["series"]] = np.nan
    ext.index.name = "date"
    ext.round(3).to_csv(out / "tmax_2026_era5_calibrated.csv")
    e5.round(3).to_csv(out / "era5_nation_daily_uncalibrated.csv")
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"ERA5 files used: {[f.name for f in files]}")


if __name__ == "__main__":
    main(*sys.argv[1:6])
