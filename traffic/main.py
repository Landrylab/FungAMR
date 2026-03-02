#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import json
import os
from pathlib import Path
from datetime import datetime, timezone
import sys
import time
import urllib.error
import urllib.request

# ========= Env variables =========
TOKEN = os.getenv("TRAFFIC_ACTION_TOKEN")
OWNER = os.getenv("TRAFFIC_ACTION_OWNER")
REPO  = os.getenv("TRAFFIC_ACTION_REPO")

if not TOKEN or not OWNER or not REPO:
    print("ERROR: Missing one or more env vars: "
          "TRAFFIC_ACTION_TOKEN, TRAFFIC_ACTION_OWNER, TRAFFIC_ACTION_REPO",
          file=sys.stderr)
    sys.exit(1)

# ========= API constants =========
API_VERSION = "2022-11-28"  # pin a stable REST API version
BASE_URL = f"https://api.github.com/repos/{OWNER}/{REPO}/traffic"
HEADERS = {
    "Accept": "application/vnd.github+json",
    "Authorization": f"Bearer {TOKEN}",
    "X-GitHub-Api-Version": API_VERSION,
    "User-Agent": f"github-traffic-archiver ({OWNER}/{REPO})"
}

# ========= Output structure =========
OUT_DIR = Path("traffic")
RAW_DIR = OUT_DIR / "raw"
OUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

# keep one folder per UTC date for raw payloads
today_utc = datetime.now(timezone.utc).date().isoformat()
RAW_TODAY = RAW_DIR / today_utc
RAW_TODAY.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = OUT_DIR / "summary.csv"
CSV_FIELDS = ["date", "views", "unique_views", "clones", "unique_cloners"]

# ========= Simple HTTP GET with retries =========
def http_get(url: str, retries: int = 3, backoff: float = 2.0):
    req = urllib.request.Request(url, headers=HEADERS, method="GET")
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                code = resp.getcode()
                data = resp.read().decode("utf-8")
                if code == 200:
                    return json.loads(data)
                else:
                    # include body for diagnostics
                    raise RuntimeError(f"HTTP {code}: {data}")
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            last_err = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last_err}")

# ========= Endpoints (daily buckets for views/clones) =========
def fetch_views_daily():
    # per=day ensures daily buckets (last 14 days)
    return http_get(f"{BASE_URL}/views?per=day")  # {count, uniques, views:[{timestamp,count,uniques},...]}  # noqa

def fetch_clones_daily():
    return http_get(f"{BASE_URL}/clones?per=day")  # {count, uniques, clones:[{timestamp,count,uniques},...]}  # noqa

def fetch_popular_paths():
    return http_get(f"{BASE_URL}/popular/paths")

def fetch_popular_referrers():
    return http_get(f"{BASE_URL}/popular/referrers")

# ========= CSV I/O =========
def read_summary():
    rows = {}
    if SUMMARY_CSV.exists():
        with SUMMARY_CSV.open(newline="", encoding="utf-8") as f:
            r = csv.DictReader(f)
            for row in r:
                rows[row["date"]] = {
                    "views": int(row.get("views", 0) or 0),
                    "unique_views": int(row.get("unique_views", 0) or 0),
                    "clones": int(row.get("clones", 0) or 0),
                    "unique_cloners": int(row.get("unique_cloners", 0) or 0),
                }
    return rows

def write_summary(rows_dict):
    dates = sorted(rows_dict.keys())
    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        for d in dates:
            rec = rows_dict[d]
            w.writerow({
                "date": d,
                "views": rec.get("views", 0),
                "unique_views": rec.get("unique_views", 0),
                "clones": rec.get("clones", 0),
                "unique_cloners": rec.get("unique_cloners", 0),
            })

def to_date(iso_ts: str) -> str:
    # GitHub returns UTC timestamps aligned to midnight for daily buckets
    # e.g., "2024-06-01T00:00:00Z" -> "2024-06-01"
    return iso_ts.split("T", 1)[0]

# ========= Main =========
def main():
    try:
        views = fetch_views_daily()   # dict with keys: count, uniques, views:[...]
        clones = fetch_clones_daily() # dict with keys: count, uniques, clones:[...]
        paths  = fetch_popular_paths()
        refs   = fetch_popular_referrers()
    except RuntimeError as e:
        print("ERROR calling GitHub REST API:", e, file=sys.stderr)
        print("Tip: ensure the token has access and the caller has push/write permission to the repo.", file=sys.stderr)
        sys.exit(2)

    # save raw payloads for today (debug/trace)
    (RAW_TODAY / "views.json").write_text(json.dumps(views, indent=2))
    (RAW_TODAY / "clones.json").write_text(json.dumps(clones, indent=2))
    (RAW_TODAY / "popular_paths.json").write_text(json.dumps(paths, indent=2))
    (RAW_TODAY / "popular_referrers.json").write_text(json.dumps(refs, indent=2))

    # merge into summary.csv by date
    combined = read_summary()

    for item in views.get("views", []):
        d = to_date(item["timestamp"])
        rec = combined.get(d, {"views": 0, "unique_views": 0, "clones": 0, "unique_cloners": 0})
        rec["views"] = int(item.get("count", 0))
        rec["unique_views"] = int(item.get("uniques", 0))
        combined[d] = rec

    for item in clones.get("clones", []):
        d = to_date(item["timestamp"])
        rec = combined.get(d, {"views": 0, "unique_views": 0, "clones": 0, "unique_cloners": 0})
        rec["clones"] = int(item.get("count", 0))
        rec["unique_cloners"] = int(item.get("uniques", 0))
        combined[d] = rec

    write_summary(combined)

    # console summary
    print(f"Saved raw payloads under {RAW_TODAY}")
    print(f"Wrote {SUMMARY_CSV} with {len(combined)} days tracked.")
    for d in sorted(combined.keys())[-5:]:
        r = combined[d]
        print(f"{d}: views={r['views']} (unique={r['unique_views']}), "
              f"clones={r['clones']} (unique={r['unique_cloners']})")

if __name__ == "__main__":
    main()
