"""Season playback builder: score every day of a fall season from archived weather + cached iNat observations.

  python scripts/build_season.py --year 2024 --year 2025 --year current

Outputs site/seasons/<year>/
  timeline.json      per-day regional summary + species stats, plus the season's observations (for pins)
  <date>.bin.gz      uint8 [4 species, rows, cols] scores for that day (same layout as site/data/scores.bin.gz)
and site/seasons/index.json listing the seasons available.

Weather comes from Open-Meteo's archive (ERA5 reanalysis) rather than the forecast model used by the live map, so
playback is for timing and broad pattern, not for grading individual cells.
"""
import argparse
import datetime as dt
import gzip
import json
import os
import sys
import time
from zoneinfo import ZoneInfo

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C                      # noqa: E402
import engine as EN                     # noqa: E402
import model as M                       # noqa: E402

ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"
INAT = "https://api.inaturalist.org/v1/observations"
SEASONS_DIR = os.path.join(C.ROOT, "site", "seasons")
WX_CACHE_DIR = os.path.join(C.DATA_STATIC, "..", "seasons")     # data/seasons/weather_<year>.npz (committed)
# Archive requests are weighted by points x days against Open-Meteo's free hourly quota, so playback uses a
# 0.2 deg grid (195 points, 2 requests per season; ERA5 is ~25 km anyway) and caches past seasons on disk.
AR_STEP = 0.2
AR_LATS = np.round(np.arange(C.LAT_MIN, 47.0 + 1e-9, AR_STEP), 4)      # 13
AR_LONS = np.round(np.arange(C.LON_MIN, C.LON_MAX + 1e-9, AR_STEP), 4)  # 15
AR_POINTS = [(la, lo) for la in AR_LATS for lo in AR_LONS]
PLAYBACK_START = (8, 15)     # first day shown
WX_LEAD_DAYS = 31            # history needed before the first day (30-day rain)


def fetch_archive(start, end):
    pts = AR_POINTS
    print(f"[archive] {len(pts)} points {start}..{end}")
    results = []
    for i in range(0, len(pts), 100):
        chunk = pts[i:i + 100]
        d = C.http_json(ARCHIVE, dict(
            latitude=",".join(f"{p[0]:.2f}" for p in chunk), longitude=",".join(f"{p[1]:.2f}" for p in chunk),
            start_date=start.isoformat(), end_date=end.isoformat(),
            daily="precipitation_sum,temperature_2m_max,temperature_2m_min", timezone="America/Los_Angeles"), timeout=180)
        results += [d] if isinstance(d, dict) else d
        time.sleep(8.0)          # archive requests are weighted heavily against the per-minute quota
    return EN.parse_weather(results, len(pts), AR_LATS, AR_LONS)


def season_weather(year, current, wx_start, last):
    """Past seasons: load data/seasons/weather_<year>.npz if present, else fetch once and save it (commit it)."""
    path = os.path.normpath(os.path.join(WX_CACHE_DIR, f"weather_{year}.npz"))
    if not current and os.path.exists(path):
        z = np.load(path, allow_pickle=False)
        dates = [dt.date.fromisoformat(d) for d in z["dates"]]
        print(f"[archive] {year}: cached {path} ({dates[0]}..{dates[-1]})")
        return dates, z["P"].astype(np.float32), z["TX"].astype(np.float32), z["TN"].astype(np.float32), z["E"].astype(np.float32)
    dates, P, TX, TN, E = fetch_archive(wx_start, last)
    if not current:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        np.savez_compressed(path, dates=np.array([d.isoformat() for d in dates]), P=P.astype(np.float16),
                            TX=TX.astype(np.float16), TN=TN.astype(np.float16), E=E.astype(np.float16))
        print(f"[archive] saved {path} ({os.path.getsize(path) / 1024:.0f} KB) - commit it")
    return dates, P, TX, TN, E


def load_observations(year, current):
    """All cached observations (+ a fresh pull for the current season)."""
    hist = json.load(open(os.path.join(C.DATA_STATIC, "inat_historical.json")))
    obs = {k: {o[0]: o for o in v} for k, v in hist["species"].items()}     # id -> [id, date, lat, lon, obsc, grade]
    if current:
        d1 = dt.date(year, 8, 1).isoformat()
        for key, sp in M.SPECIES.items():
            id_above = 0
            try:
                while True:
                    d = C.http_json(INAT, dict(
                        taxon_id=",".join(map(str, sp["inat_taxa"])), nelat=C.LAT_MAX, nelng=C.LON_MAX,
                        swlat=C.LAT_MIN, swlng=C.LON_MIN, quality_grade="research,needs_id", d1=d1,
                        per_page=200, order_by="id", order="asc", id_above=id_above), tries=3)
                    res = d["results"]
                    for r in res:
                        g = r.get("geojson")
                        if not g or not r.get("observed_on"):
                            continue
                        obs[key][r["id"]] = [r["id"], r["observed_on"], round(g["coordinates"][1], 4),
                                             round(g["coordinates"][0], 4),
                                             1 if (r.get("obscured") or r.get("geoprivacy") == "obscured"
                                                   or r.get("taxon_geoprivacy") == "obscured") else 0,
                                             r["quality_grade"][0]]
                    if len(res) < 200:
                        break
                    id_above = res[-1]["id"]
                    time.sleep(1.0)
            except Exception as e:                      # noqa: BLE001
                print(f"WARNING: iNat pull for current season failed ({e}); using cache only", file=sys.stderr)
            time.sleep(1.0)
    return {k: sorted(v.values(), key=lambda o: o[1]) for k, v in obs.items()}


def build(year, current, static, out_root):
    today = dt.datetime.now(ZoneInfo("America/Los_Angeles")).date()
    first = dt.date(year, *PLAYBACK_START)
    last = dt.date(year, 12, 31) if not current else today - dt.timedelta(days=1)
    if last < first:
        print(f"[season {year}] nothing to build yet (starts {first})")
        return None
    wx_start = first - dt.timedelta(days=WX_LEAD_DAYS)
    dates, P, TX, TN, E = season_weather(year, current, wx_start, last)
    wx = EN.CellWeather(static, dates, P, TX, TN, E, AR_LATS, AR_LONS)
    obs = load_observations(year, current)
    obs_ord = {k: np.array([dt.date.fromisoformat(o[1]).toordinal() for o in v]) for k, v in obs.items()}
    season_obs = []
    for k, v in obs.items():
        for o in v:
            d = dt.date.fromisoformat(o[1])
            if first - dt.timedelta(days=21) <= d <= last:
                season_obs.append([o[0], o[1], k, o[2], o[3], o[4], o[5]])
    season_obs.sort(key=lambda o: o[1])

    out_dir = os.path.join(out_root, str(year))
    os.makedirs(out_dir, exist_ok=True)
    land = ~static.ocean
    days, total_bytes = [], 0
    day = first
    while day <= last:
        if day not in dates:
            break
        cond = wx.conditions(day)
        recent_n = EN.count_recent(obs_ord, day)
        res = EN.score_day(static, cond, day, recent_n)
        data = gzip.compress(np.ascontiguousarray(res["scores"]).tobytes(), 6)
        with open(os.path.join(out_dir, f"{day.isoformat()}.bin.gz"), "wb") as f:
            f.write(data)
        total_bytes += len(data)
        it = dates.index(day)
        summ = EN.species_summary(static, res)
        rec = {
            "date": day.isoformat(),
            "rain_day": round(float(wx.P[it - 1][land].mean()), 1),        # previous complete day
            "rain14": round(float(cond["rain14"][land].mean()), 1),
            "rain30": round(float(cond["rain30"][land].mean()), 1),
            "temp7": round(float(cond["tmean7"][land].mean()), 1),
            "trig_pct": round(float((~np.isnan(cond["ds_first"]))[land].mean() * 100), 1),
            "frost_pct": round(float((cond["frost_days"] >= M.FROST_DAYS)[land].mean() * 100), 1),
            "species": {k: {"obs": int(np.sum(obs_ord[k] == day.toordinal())), "recent21": recent_n[k],
                            "max": summ[k]["score_max"], "mean": summ[k]["score_mean_habitat"],
                            "ge50": summ[k]["cells_ge50"]} for k in M.SPECIES_ORDER},
        }
        days.append(rec)
        day += dt.timedelta(days=1)
    trig_date = next((d["date"] for d in days if d["trig_pct"] >= 50), None)
    timeline = {
        "year": year, "current": current, "start": days[0]["date"], "end": days[-1]["date"],
        "generated_at": dt.datetime.now(ZoneInfo("America/Los_Angeles")).isoformat(timespec="seconds"),
        "first_trigger_date": trig_date,
        "species": M.SPECIES_ORDER,
        "days": days,
        "obs_fields": ["id", "date", "species", "lat", "lon", "obscured", "grade"],
        "obs": season_obs,
        "note": "Weather from Open-Meteo ERA5 archive; observations as currently on iNaturalist (incl. later uploads).",
    }
    json.dump(timeline, open(os.path.join(out_dir, "timeline.json"), "w"), separators=(",", ":"))
    peak = max(days, key=lambda d: sum(s["max"] for s in d["species"].values()))
    print(f"[season {year}] {len(days)} days {days[0]['date']}..{days[-1]['date']}, {len(season_obs)} obs, "
          f"{total_bytes / 1e6:.1f} MB of score files; first trigger {trig_date}; peak day {peak['date']} "
          + " ".join(f"{k[:4]}={peak['species'][k]['max']}" for k in M.SPECIES_ORDER))
    return {"year": year, "current": current, "start": days[0]["date"], "end": days[-1]["date"],
            "days": len(days), "first_trigger_date": trig_date}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", action="append", required=True, help="season year (Aug-Dec) or 'current'; repeatable")
    ap.add_argument("--out", default=SEASONS_DIR)
    a = ap.parse_args()
    today = dt.datetime.now(ZoneInfo("America/Los_Angeles")).date()
    cur_year = EN.season_year(today)
    static = EN.Static()
    index_path = os.path.join(a.out, "index.json")
    index = json.load(open(index_path)) if os.path.exists(index_path) else {"seasons": []}
    failed = []
    for y in a.year:
        year = cur_year if y == "current" else int(y)
        try:
            entry = build(year, year == cur_year, static, a.out)
        except Exception as e:                      # noqa: BLE001
            print(f"ERROR: season {year} failed: {e}", file=sys.stderr)
            failed.append(year)
            continue
        index["seasons"] = [s for s in index["seasons"] if s["year"] != year] + ([entry] if entry else [])
        index["seasons"].sort(key=lambda s: s["year"])
        index["generated_at"] = dt.datetime.now(ZoneInfo("America/Los_Angeles")).isoformat(timespec="seconds")
        os.makedirs(a.out, exist_ok=True)
        json.dump(index, open(index_path, "w"), indent=1)
    print(f"[done] {index_path}: {[s['year'] for s in index['seasons']]}" + (f"; FAILED {failed}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
