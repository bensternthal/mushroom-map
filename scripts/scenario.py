"""What-if runner for tuning: score the static layers under synthetic conditions (no network).

  python scripts/scenario.py --rain14 60 --rain30 120 --temp 14 --doy 288 --first 12 --latest 3 --recent 5
  python scripts/scenario.py --frost 4          # sustained hard frost everywhere

--temp is the lowland (sea-level) 7-day mean air temperature; cells are lapse-adjusted from it.
Prints the score distribution per species; --png writes data/static/preview/scenario_<species>.png.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C                      # noqa: E402
import model as M                       # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rain14", type=float, default=60)
    ap.add_argument("--rain30", type=float, default=120)
    ap.add_argument("--temp", type=float, default=14, help="lowland 7-day mean air temp, C")
    ap.add_argument("--doy", type=int, default=288, help="day of year (288 = Oct 15)")
    ap.add_argument("--first", type=float, default=12, help="days since first season trigger (-1 = none)")
    ap.add_argument("--latest", type=float, default=3, help="days since latest big rain (-1 = none)")
    ap.add_argument("--recent", type=int, default=0, help="recent iNat sightings (21 d) for every species")
    ap.add_argument("--frost", type=int, default=0, help="frost days in the last 14 (>=3 kills)")
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()

    L = np.load(os.path.join(C.DATA_STATIC, "layers.npz"))
    kde = np.load(os.path.join(C.DATA_STATIC, "inat_kde.npz"))["kde"].astype(np.float32) / 255
    elev = np.maximum(L["elev"].astype(np.float32), 0)
    tf, eco, ck, ocean = L["tree_pct"] / 100.0, L["eco_class"], L["coast_km"].astype(np.float32), L["ocean"]
    bf = L["built_pct"] / 100.0 if "built_pct" in L else np.zeros_like(tf)
    full = lambda v: np.full(elev.shape, np.nan if v < 0 else v, np.float32)   # noqa: E731
    t7 = a.temp - M.LAPSE_C_PER_KM * elev / 1000
    print(f"rain14={a.rain14} rain30={a.rain30} lowlandT={a.temp} doy={a.doy} first={a.first} latest={a.latest} recent={a.recent} frost={a.frost}")
    for si, k in enumerate(M.SPECIES_ORDER):
        sp = M.SPECIES[k]
        hab = M.habitat(sp, elev, tf, eco, ck, bf)
        dyn = M.dynamic_score(sp, M.rain_score(sp, full(a.rain14), full(a.rain30)),
                              M.temp_score(sp, t7, np.full(elev.shape, a.frost)),
                              M.phenology(sp, a.doy, elev, a.recent),
                              M.trigger_score(sp, full(a.first), full(a.latest)))
        s = M.total_score(hab, M.prior_from_density(kde[si]), dyn)
        s[ocean] = 0
        land = ~ocean & (hab > 0.05)
        print(f"  {k:18s} max {s.max():3d}  p50 {np.percentile(s[land], 50):4.0f}  p90 {np.percentile(s[land], 90):4.0f}"
              f"  >=15 {(s >= 15).sum():6d}  >=50 {(s >= 50).sum():6d}")
        if a.png:
            from build_static import preview
            preview(f"scenario_{k}", s, "YlOrRd", 0, 100)


if __name__ == "__main__":
    main()
