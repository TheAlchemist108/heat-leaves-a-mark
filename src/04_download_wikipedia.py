"""Download Wikipedia pageviews for the pre-registered articles (PREREGISTRATION.md, Section 6.2).

(a) UK-only, primary: Wikimedia differentially private daily pageviews by country, project and page.
    One file per day for all countries is streamed; only rows with country code GB, project
    en.wikipedia and one of the article page IDs are kept. The SHA-256 of every daily source file
    is logged so the filtering can be checked against the published files.
      historical release: 2017-02-09 to 2023-02-05
      current release:    2023-02-06 onwards
(b) Global, secondary: Wikimedia REST API per-article daily pageviews (all-access, user).
    Page titles are resolved from page IDs, and page moves during the sample are recorded.

Usage:
  python src/04_download_wikipedia.py uk     <out_dir> [workers]
  python src/04_download_wikipedia.py global <out_dir>
"""
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import requests

UA = {"User-Agent": "heat-leaves-a-mark research project (https://github.com/TheAlchemist108/heat-leaves-a-mark)"}
ARTICLES = {  # page ID: family (Section 6.2)
    293270: "salience", 319671: "salience",
    7221088: "coping", 24999087: "coping", 12877572: "coping",
    5042951: "concern", 10031054: "concern", 10031822: "concern",
    68316: "transition", 5134: "placebo",
}
START, END = date(2021, 11, 7), date(2026, 10, 31)
CURRENT_FROM = date(2023, 2, 6)
BASE = "https://analytics.wikimedia.org/published/datasets/"


def stamp():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def day_url(d):
    folder = "country_project_page" if d >= CURRENT_FROM else "country_project_page_historical"
    return f"{BASE}{folder}/{d.isoformat()}.tsv"


def fetch_day(d, session):
    url = day_url(d)
    for attempt in range(5):
        try:
            r = session.get(url, headers=UA, stream=True, timeout=(15, 120))
            if r.status_code == 404:
                return d, url, "missing", "", []
            r.raise_for_status()
            h, buf, rows = hashlib.sha256(), b"", []
            for chunk in r.iter_content(1 << 20):
                h.update(chunk)
                buf += chunk
                lines = buf.split(b"\n")
                buf = lines.pop()
                for line in lines:
                    if b"\tGB\ten.wikipedia\t" not in line:
                        continue
                    f = line.decode("utf-8").split("\t")
                    if int(f[3]) in ARTICLES:
                        rows.append(f)
            if buf and b"\tGB\ten.wikipedia\t" in buf:
                f = buf.decode("utf-8").split("\t")
                if int(f[3]) in ARTICLES:
                    rows.append(f)
            return d, url, "ok", h.hexdigest(), rows
        except Exception as e:
            time.sleep(10 * (attempt + 1))
    return d, url, "failed", "", []


def uk(out, workers=6):
    out.mkdir(parents=True, exist_ok=True)
    rows_path, man_path = out / "uk_dp_rows.tsv", out / "uk_dp_manifest.csv"
    done = set()
    if man_path.exists():
        with man_path.open() as f:
            done = {r["date"] for r in csv.DictReader(f) if r["status"] in ("ok", "missing")}
    else:
        man_path.write_text("date,source_url,status,source_sha256,rows_kept,downloaded_utc\n")
        rows_path.write_text("date\tcountry\tcountry_code\tproject\tpage_id\tpage_title\titem_id\tviews\n")
    days = [START + timedelta(n) for n in range((END - START).days + 1)]
    todo = [d for d in days if d.isoformat() not in done]
    print(f"{len(todo)} days to fetch", flush=True)
    session = requests.Session()
    with ThreadPoolExecutor(workers) as ex, rows_path.open("a") as rf, man_path.open("a") as mf:
        futures = [ex.submit(fetch_day, d, session) for d in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            d, url, status, sha, rows = fut.result()
            for f in rows:
                rf.write(d.isoformat() + "\t" + "\t".join(x.strip() for x in f) + "\n")
            mf.write(f"{d.isoformat()},{url},{status},{sha},{len(rows)},{stamp()}\n")
            rf.flush(); mf.flush()
            if i % 50 == 0:
                print(f"{i}/{len(todo)} days", flush=True)
    print("done", flush=True)


def api(params):
    params = {**params, "format": "json", "formatversion": 2}
    r = requests.get("https://en.wikipedia.org/w/api.php", params=params, headers=UA, timeout=60)
    r.raise_for_status()
    return r.json()


def global_views(out):
    out.mkdir(parents=True, exist_ok=True)
    ids = "|".join(str(i) for i in ARTICLES)
    pages = api({"action": "query", "pageids": ids})["query"]["pages"]
    titles = {p["pageid"]: p["title"] for p in pages}
    moves = {}
    for pid in titles:
        # page moves leave a revision whose comment starts "moved page"; scan the sample period
        found, cont = [], {}
        while True:
            q = api({"action": "query", "prop": "revisions", "pageids": pid, "rvprop": "timestamp|comment",
                     "rvlimit": "max", "rvstart": START.isoformat() + "T00:00:00Z",
                     "rvend": END.isoformat() + "T23:59:59Z", "rvdir": "newer", **cont})
            for rev in q["query"]["pages"][0].get("revisions", []):
                if "moved page" in rev.get("comment", "").lower():
                    found.append(rev)
            if "continue" not in q:
                break
            cont = q["continue"]
        moves[pid] = found
    raw = {"titles": titles, "moves_in_sample": moves, "retrieved_utc": stamp(), "articles": {}}
    for pid, title in titles.items():
        t = title.replace(" ", "_")
        url = ("https://wikimedia.org/api/rest_v1/metrics/pageviews/per-article/en.wikipedia/"
               f"all-access/user/{requests.utils.quote(t, safe='')}/daily/"
               f"{START.strftime('%Y%m%d')}00/{END.strftime('%Y%m%d')}00")
        r = requests.get(url, headers=UA, timeout=60)
        r.raise_for_status()
        raw["articles"][str(pid)] = {"title": title, "url": url, "response": r.json()}
        time.sleep(1)
    path = out / "global_rest_pageviews.json"
    data = json.dumps(raw, ensure_ascii=False, indent=1).encode("utf-8")
    path.write_bytes(data)
    print(f"saved {path} sha256={hashlib.sha256(data).hexdigest()}")
    for pid, title in titles.items():
        print(pid, ARTICLES[pid], title, "moves in sample:", len(moves[pid]))


if __name__ == "__main__":
    out = Path(sys.argv[2])
    if sys.argv[1] == "uk":
        uk(out, int(sys.argv[3]) if len(sys.argv) > 3 else 6)
    else:
        global_views(out)
