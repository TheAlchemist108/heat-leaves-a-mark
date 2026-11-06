"""Download ERA5 daily maximum 2 m temperature over the UK (PREREGISTRATION.md, Section 4.1).

Dataset: "ERA5 post-processed daily statistics on single levels from 1940 to present"
(derived-era5-single-levels-daily-statistics), Copernicus Climate Data Store.
  - 2026, January-October: the treatment for 2026 (HadUK-Grid does not yet cover 2026)
  - 2016-2025, April-October: calibration against HadUK-Grid

Requests are submitted without waiting, their IDs saved, and results collected later, so the
script can be run in short sessions:
  python src/03_download_era5.py submit  <out_dir>
  python src/03_download_era5.py collect <out_dir>

The API key is read from the file named by the CDSAPI_RC environment variable (or ~/.cdsapirc);
it is never printed or saved elsewhere.
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from ecmwf.datastores import Client  # installed with cdsapi >= 0.7.7

DATASET = "derived-era5-single-levels-daily-statistics"
AREA = [61, -9, 49.5, 2]  # N, W, S, E: all four UK nations
DAYS = [f"{d:02d}" for d in range(1, 32)]


def requests():
    for year in range(2016, 2026):
        yield f"era5_tmax_uk_{year}_04-10", year, [f"{m:02d}" for m in range(4, 11)]
    yield "era5_tmax_uk_2026_01-10", 2026, [f"{m:02d}" for m in range(1, 11)]


def body(year, months):
    return {"product_type": "reanalysis", "variable": ["2m_temperature"], "year": str(year),
            "month": months, "day": DAYS, "daily_statistic": "daily_maximum",
            "time_zone": "utc+00:00", "frequency": "1_hourly", "area": AREA}


def client():
    rc = Path(os.environ.get("CDSAPI_RC", Path.home() / ".cdsapirc"))
    cfg = dict(line.split(":", 1) for line in rc.read_text(encoding="utf-8-sig").splitlines() if ":" in line)
    return Client(url=cfg["url"].strip(), key=cfg["key"].strip(), progress=False)


def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def submit(out):
    state_file = out / "requests.json"
    state = json.loads(state_file.read_text()) if state_file.exists() else {}
    c = client()
    for name, year, months in requests():
        if name in state or (out / f"{name}.nc").exists():
            continue
        req = body(year, months)
        remote = c.submit(DATASET, req)
        state[name] = {"request_id": remote.request_id, "request": req, "submitted_utc": stamp()}
        state_file.write_text(json.dumps(state, indent=1))
        print(f"submitted {name}: {state[name]['request_id']}", flush=True)


def collect(out):
    state = json.loads((out / "requests.json").read_text())
    c = client()
    for name, info in state.items():
        path = out / f"{name}.nc"
        if path.exists():
            continue
        remote = c.get_remote(info["request_id"])
        status = remote.status
        if status in ("successful", "completed"):
            remote.download(str(path))
            info["downloaded_utc"] = stamp()
            info["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            print(f"downloaded {name}", flush=True)
        else:
            print(f"{name}: {status}", flush=True)
    (out / "requests.json").write_text(json.dumps(state, indent=1))
    done = sum((out / f"{n}.nc").exists() for n in state)
    print(f"{done}/{len(state)} files downloaded")


if __name__ == "__main__":
    out = Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    {"submit": submit, "collect": collect}[sys.argv[1]](out)
