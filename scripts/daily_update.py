"""Daily pipeline (GitHub Actions): weather + recent iNat -> score grids -> site/data/.

Outputs (all in site/data/):
  meta.json            grid def, timestamps, species list, legend, factor byte layout, regional summary
  scores.bin.gz        uint8 [4 species, rows, cols] row-major, 0-100 relative suitability
  factors.bin.gz       uint8 [rows, cols, FACTOR_BYTES] per-cell explanation bytes (see FACTOR_LAYOUT)
  recent_obs.geojson   iNat observations from the last 30 days
  public_lands.geojson copied from data/static
Exit status is non-zero if weather cannot be fetched.  A failed iNat fetch reuses the previous
recent_obs.geojson (or writes an empty one) with a warning.
"""
import datetime as dt
import gzip
import json
import os
import shutil
import sys
import time
from zoneinfo import ZoneInfo

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C                      # noqa: E402
import engine as EN                     # noqa: E402
import model as M                       # noqa: E402

TZ = ZoneInfo("America/Los_Angeles")
OPEN_METEO = "https://api.open-meteo.com/v1/forecast"
INAT = "https://api.inaturalist.org/v1/observations"
NS = EN.NS

# per-cell explanation bytes
FACTOR_LAYOUT = [
    ["rain14_mm", "mm, capped 255"],
    ["rain30_mm_x2", "mm / 2, capped 255"],
    ["fc3_mm", "forecast precip next 3 days, mm, capped 255"],
    ["temp7_enc", "7-day mean air temp at cell elevation: C = v/4 - 20"],
    ["days_since_first_trigger", "255 = none yet this season"],
    ["days_since_latest_trigger", "255 = none"],
    ["elev_x16", "metres / 16"],
    ["eco_class", "index into meta.eco_classes"],
    ["flags", "bit0 frost-killed, bit1 ocean, bit2 no-ecoregion"],
    ["tree_pct", "tree cover %"],
] + [[f"{k}_dyn", "weather x phenology 0-100"] for k in M.SPECIES_ORDER] \
  + [[f"{k}_hab", "habitat 0-100"] for k in M.SPECIES_ORDER] \
  + [[f"{k}_prior", "iNat prior 0-100 (50 = no history)"] for k in M.SPECIES_ORDER]
FACTOR_BYTES = len(FACTOR_LAYOUT)


# ---------------------------------------------------------------------------
# weather
# ---------------------------------------------------------------------------
def fetch_weather(today):
    season_start = dt.date(EN.season_year(today), *M.SEASON_START)
    past_days = int(np.clip((today - season_start).days + 5, 31, 92))
    pts = EN.WX_POINTS
    print(f"[weather] {len(pts)} points, past_days={past_days}, forecast_days=3")
    results = []
    for i in range(0, len(pts), 100):
        chunk = pts[i:i + 100]
        d = C.http_json(OPEN_METEO, dict(
            latitude=",".join(f"{p[0]:.2f}" for p in chunk), longitude=",".join(f"{p[1]:.2f}" for p in chunk),
            daily="precipitation_sum,temperature_2m_max,temperature_2m_min",
            past_days=past_days, forecast_days=3, timezone="America/Los_Angeles"), timeout=120)
        results += [d] if isinstance(d, dict) else d
        time.sleep(0.5)
    return EN.parse_weather(results, len(pts))


# ---------------------------------------------------------------------------
# iNaturalist recent
# ---------------------------------------------------------------------------
def fetch_recent_inat(today, days=30):
    d1 = (today - dt.timedelta(days=days)).isoformat()
    feats = []
    for key, sp in M.SPECIES.items():
        id_above = 0
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
                obscured = bool(r.get("obscured") or r.get("geoprivacy") == "obscured"
                                or r.get("taxon_geoprivacy") == "obscured")
                feats.append({"type": "Feature",
                              "geometry": {"type": "Point", "coordinates": [round(g["coordinates"][0], 4),
                                                                            round(g["coordinates"][1], 4)]},
                              "properties": {"species": key, "date": r["observed_on"], "obscured": obscured,
                                             "grade": r["quality_grade"], "taxon": r["taxon"]["name"],
                                             "url": r["uri"], "id": r["id"]}})
            if len(res) < 200:
                break
            id_above = res[-1]["id"]
            time.sleep(1.0)
        time.sleep(1.0)
    return {"type": "FeatureCollection", "features": feats,
            "window_days": days, "fetched": today.isoformat()}


# ---------------------------------------------------------------------------
def main():
    now = dt.datetime.now(TZ)
    today = now.date()
    doy = today.timetuple().tm_yday
    os.makedirs(C.SITE_DATA, exist_ok=True)

    static = EN.Static()
    grid = json.load(open(os.path.join(C.DATA_STATIC, "grid.json")))
    elev, ocean, eco_class = static.elev, static.ocean, static.eco_class

    # ---- weather (fatal on failure) --------------------------------------
    dates, P, TX, TN, E = fetch_weather(today)
    if today not in dates:
        raise C.FetchError(f"weather: today {today} missing from date axis {dates[0]}..{dates[-1]}")
    it = dates.index(today)
    print(f"[weather] {dates[0]} .. {dates[-1]}; {it} history days, today={today}")
    wx = EN.CellWeather(static, dates, P, TX, TN, E)
    cond = wx.conditions(today)
    rain14, rain30, fc3, Tmean7 = cond["rain14"], cond["rain30"], cond["fc3"], cond["tmean7"]
    frost_days, ds_first, ds_latest = cond["frost_days"], cond["ds_first"], cond["ds_latest"]

    # ---- iNat recent (non-fatal) ------------------------------------------
    recent_path = os.path.join(C.SITE_DATA, "recent_obs.geojson")
    try:
        recent = fetch_recent_inat(today)
        json.dump(recent, open(recent_path, "w"), separators=(",", ":"))
        print(f"[inat] {len(recent['features'])} observations in last 30 days")
    except Exception as e:                      # noqa: BLE001
        print(f"WARNING: iNat recent fetch failed ({e}); reusing cached recent_obs.geojson", file=sys.stderr)
        if os.path.exists(recent_path):
            recent = json.load(open(recent_path))
        else:
            recent = {"type": "FeatureCollection", "features": [], "stale": True}
            json.dump(recent, open(recent_path, "w"))
    cutoff = (today - dt.timedelta(days=M.RECENT_DAYS)).isoformat()
    recent_n = {k: sum(1 for f in recent["features"]
                       if f["properties"]["species"] == k and f["properties"]["date"] >= cutoff)
                for k in M.SPECIES_ORDER}

    # ---- scores -------------------------------------------------------------
    res = EN.score_day(static, cond, today, recent_n)
    scores = res["scores"]
    dyn_b = np.round(res["dyn"] * 100).astype(np.uint8)
    hab_b = np.round(static.habitat * 100).astype(np.uint8)
    pri_b = np.round(static.prior * 100).astype(np.uint8)
    summary = EN.species_summary(static, res)
    for key in M.SPECIES_ORDER:
        summary[key]["recent_obs_21d"] = recent_n[key]
        sm = summary[key]
        print(f"[score] {key:18s} max {sm['score_max']:3d}  mean(habitat cells) {sm['score_mean_habitat']:5.1f}  "
              f">=15: {sm['cells_ge15']:6d}  >=50: {sm['cells_ge50']:5d}  recent obs {recent_n[key]}")

    # ---- factors -------------------------------------------------------------
    F = np.zeros((C.NROWS, C.NCOLS, FACTOR_BYTES), np.uint8)
    F[..., 0] = np.clip(np.round(rain14), 0, 255)
    F[..., 1] = np.clip(np.round(rain30 / 2), 0, 255)
    F[..., 2] = np.clip(np.round(fc3), 0, 255)
    F[..., 3] = np.clip(np.round((Tmean7 + 20) * 4), 0, 255)
    F[..., 4] = np.where(np.isnan(ds_first), 255, np.clip(np.nan_to_num(ds_first), 0, 254))
    F[..., 5] = np.where(np.isnan(ds_latest), 255, np.clip(np.nan_to_num(ds_latest), 0, 254))
    F[..., 6] = np.clip(np.round(elev / 16), 0, 255)
    F[..., 7] = eco_class
    F[..., 8] = (frost_days >= M.FROST_DAYS).astype(np.uint8) | (ocean.astype(np.uint8) << 1) | ((eco_class == 0).astype(np.uint8) << 2)
    F[..., 9] = static.tree_pct
    for si in range(NS):
        F[..., 10 + si] = dyn_b[si]
        F[..., 10 + NS + si] = hab_b[si]
        F[..., 10 + 2 * NS + si] = pri_b[si]

    # ---- write ----------------------------------------------------------------
    def wgz(name, arr):
        data = gzip.compress(np.ascontiguousarray(arr).tobytes(), 9)
        with open(os.path.join(C.SITE_DATA, name), "wb") as f:
            f.write(data)
        print(f"  wrote {name}: {len(data) / 1024:.0f} KB gz ({arr.nbytes / 1024:.0f} KB raw)")

    wgz("scores.bin.gz", scores)
    wgz("factors.bin.gz", F)
    shutil.copy(os.path.join(C.DATA_STATIC, "public_lands.geojson"), os.path.join(C.SITE_DATA, "public_lands.geojson"))
    land_all = ~ocean
    meta = {
        "generated_at": now.isoformat(timespec="seconds"),
        "date": today.isoformat(), "doy": doy,
        "weather": {"history_start": dates[0].isoformat(), "history_end": dates[it - 1].isoformat(),
                    "forecast_end": dates[-1].isoformat(), "source": "Open-Meteo (0.1 deg grid, lapse-adjusted)",
                    "regional_mean": {"rain14_mm": round(float(rain14[land_all].mean()), 1),
                                      "rain30_mm": round(float(rain30[land_all].mean()), 1),
                                      "fc3_mm": round(float(fc3[land_all].mean()), 1),
                                      "temp7_c": round(float(Tmean7[land_all].mean()), 1)},
                    "cells_triggered_pct": round(float((~np.isnan(ds_first))[land_all].mean() * 100), 1)},
        "grid": grid["grid"],
        "species": [{"key": k, "index": i, "name": M.SPECIES[k]["name"], "latin": M.SPECIES[k]["latin"],
                     "color": M.SPECIES[k]["color"], "inat_taxa": M.SPECIES[k]["inat_taxa"],
                     "summary": summary[k],
                     "params": {"elev_zones": [[z["label"], z["elev"]] for z in M.SPECIES[k]["zones"]],
                                "rain14_mm": M.SPECIES[k]["rain14"], "temp_c": M.SPECIES[k]["temp"],
                                "peak_doy": M.SPECIES[k]["pheno"]["peak_doy"],
                                "trigger": M.SPECIES[k]["trigger"]}}
                    for i, k in enumerate(M.SPECIES_ORDER)],
        "score": {"range": [0, 100], "min_display": M.SCORE_MIN_DISPLAY,
                  "note": "Relative suitability index, not a probability. 100 = best conditions the model can express."},
        "eco_classes": grid["eco_classes"], "eco_labels": grid["eco_labels"],
        "factors": {"bytes_per_cell": FACTOR_BYTES, "layout": FACTOR_LAYOUT},
        "files": {"scores": "scores.bin.gz", "factors": "factors.bin.gz",
                  "recent_obs": "recent_obs.geojson", "public_lands": "public_lands.geojson",
                  "seasons_index": "../seasons/index.json"},
        "inat_historical": {"fetched": grid["inat_fetched"], "counts": grid["inat_counts"]},
        "model": {"gamma": M.GAMMA, "trigger_mm_48h": M.TRIGGER_MM_48H, "season_start": list(M.SEASON_START),
                  "frost": [M.FROST_TEMP_C, M.FROST_DAYS, M.FROST_WINDOW], "kde_sigma_km": M.KDE_SIGMA_KM},
    }
    json.dump(meta, open(os.path.join(C.SITE_DATA, "meta.json"), "w"), indent=1)
    total = sum(os.path.getsize(os.path.join(C.SITE_DATA, f)) for f in os.listdir(C.SITE_DATA))
    print(f"[done] site/data total {total / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
