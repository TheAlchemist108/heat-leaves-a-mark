"""Download one complete Google Trends pull (PREREGISTRATION.md, Section 6.1).

Each request is saved as the raw JSON returned by Google, unedited, together with the request
parameters and download time. Re-running the same pull resumes and skips files already saved,
so the pull can be run in several short sessions. Google rate-limits heavily; on a limit the
script waits and retries.

Usage:
  python src/02_download_trends.py <pull_number> <out_dir> [time_budget_seconds]
  python src/02_download_trends.py manifest <pull_number> <out_dir>
"""
import os
import hashlib
import json
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from pytrends.request import TrendReq

REQUESTS = {
    "R1": ["dentist", "heatwave", "climate change", "air conditioning", "car insurance"],
    "R2": ["dentist", "global warming", "net zero", "portable air conditioner", "heat stroke"],
    "R3": ["dentist", "heat pump", "recipes", "plumber", "climate change"],
}
GEOS = ["GB", "GB-ENG", "GB-SCT", "GB-WLS", "GB-NIR"]
WEEKLY_WINDOW = "2021-11-07 2026-10-31"
DAILY_WINDOWS = {
    "2022": "2022-04-01 2022-12-15",
    "2023": "2023-04-01 2023-12-15",
    "2024": "2024-04-01 2024-12-15",
    "2025": "2025-04-01 2025-12-15",
    "2026": "2026-04-01 2026-11-03",
}


def jobs():
    for r, terms in REQUESTS.items():
        for geo in GEOS:
            yield {"kind": "weekly", "request": r, "terms": terms, "geo": geo,
                   "window": "all", "timeframe": WEEKLY_WINDOW}
    for r in ("R1", "R2"):
        for label, tf in DAILY_WINDOWS.items():
            yield {"kind": "daily", "request": r, "terms": REQUESTS[r], "geo": "GB",
                   "window": label, "timeframe": tf}


def fetch(pt, job):
    pt.build_payload(job["terms"], cat=0, timeframe=job["timeframe"], geo=job["geo"], gprop="")
    widget = pt.interest_over_time_widget
    params = {"req": json.dumps(widget["request"]), "token": widget["token"], "tz": pt.tz}
    return pt._get_data(url=TrendReq.INTEREST_OVER_TIME_URL, method=TrendReq.GET_METHOD,
                        trim_chars=5, params=params)


def fname(job):
    return f"{job['kind']}_{job['request']}_{job['geo']}_{job['window']}.json"


def main(pull, out_dir, budget=None):
    start = time.monotonic()
    left = lambda: float("inf") if budget is None else budget - (time.monotonic() - start)
    out = Path(out_dir) / f"pull{pull}"
    out.mkdir(parents=True, exist_ok=True)
    pt = TrendReq(hl="en-GB", tz=0, timeout=(10, 30))
    todo = list(jobs())
    for i, job in enumerate(todo, 1):
        path = out / fname(job)
        if path.exists():
            continue
        data = None
        while left() > 25:
            try:
                data = fetch(pt, job)
                break
            except Exception as e:  # rate limits (429) and transient errors
                wait = min(45 + random.uniform(0, 15), left() - 25)
                print(f"[{i}/{len(todo)}] {path.name}: {type(e).__name__}; waiting {max(wait, 0):.0f}s", flush=True)
                if wait > 0:
                    time.sleep(wait)
                pt = TrendReq(hl="en-GB", tz=0, timeout=(10, 30))
        if data is None:
            break
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        record = {"job": job, "downloaded_utc": stamp, "response": data}
        tmp = path.with_suffix(".tmp")
        tmp.write_bytes(json.dumps(record, ensure_ascii=False, indent=1).encode("utf-8"))
        os.replace(tmp, path)  # atomic: a file either exists complete or not at all
        print(f"[{i}/{len(todo)}] {path.name}: saved", flush=True)
        pause = random.uniform(10, 20)
        if left() - pause < 25:
            break
        time.sleep(pause)
    done = sum((out / fname(j)).exists() for j in todo)
    print(f"Pull {pull}: {done}/{len(todo)} requests saved in {out}", flush=True)


def write_manifest(pull, out_dir):
    out = Path(out_dir) / f"pull{pull}"
    rows = ["file,kind,request,geo,window,timeframe,terms,downloaded_utc,sha256"]
    for job in jobs():
        path = out / fname(job)
        if not path.exists():
            continue
        raw = path.read_bytes()
        stamp = json.loads(raw)["downloaded_utc"]
        rows.append(f"{path.name},{job['kind']},{job['request']},{job['geo']},{job['window']},"
                    f"\"{job['timeframe']}\",\"{'|'.join(job['terms'])}\",{stamp},"
                    f"{hashlib.sha256(raw).hexdigest()}")
    (out / "manifest.csv").write_text("\n".join(rows) + "\n")
    print(f"Manifest written with {len(rows) - 1} files")


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        write_manifest(int(sys.argv[2]), sys.argv[3])
    else:
        main(int(sys.argv[1]), sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else None)
