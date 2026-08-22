# PNW Mushroom Forecast Map

A static, daily-updated deck.gl map of northwest Oregon / southwest Washington (lat 44.5–47.0, lon −124.2…−121.4,
"coast to Cascades") on a ~1 km grid, colour-coded by the **relative likelihood** of finding four fall species:

| Species | Key | Habitat encoded |
|---|---|---|
| Matsutake (*Tricholoma murrillianum*, incl. *T. magnivelare* records) | `matsutake` | coastal shore-pine/dune zone (≤8 km from the ocean, <250 m) **or** Cascade hemlock / true fir 600–1600 m |
| White chanterelle (*Cantharellus subalbidus*) | `white_chanterelle` | Doug-fir / western hemlock 0–1000 m |
| Pacific golden chanterelle (*Cantharellus formosus*) | `golden_chanterelle` | Doug-fir / western hemlock 0–1200 m, Coast Range + W Cascades |
| Fall king bolete (*Boletus edulis* group) | `king_bolete` | coastal Sitka spruce / shore pine (≤10 km from ocean) **or** Cascade true fir / hemlock 800–1500 m (weighted 0.7) |

The number shown is a **0–100 relative suitability index, not a probability** and not a guarantee. It is meant to
answer "where and when is it worth looking?", and to be tuned by someone who actually picks.

Live site: GitHub Pages, rebuilt every morning by GitHub Actions (see *Setup*).

## How a cell is scored

Everything tunable lives in one dict at the top of [`scripts/model.py`](scripts/model.py) (`SPECIES`) plus a few
globals right below it. Per cell and species:

```
score = 100 · (habitat · prior · dynamic) ^ GAMMA            GAMMA = 0.5

habitat  = forest_score · max over habitat zones ( ecoregion_weight · elevation_band · coast_band )
           forest_score   ramp on ESA WorldCover tree-cover fraction (0.15 → 0.5)
           ecoregion      EPA Level IV ecoregion mapped to an assemblage class (Sitka spruce zone, western hemlock
                          zone, silver fir zone, Willamette Valley floor ≈ 0, …)
           elevation_band trapezoid in metres, per zone
           coast_band     trapezoid on distance-to-ocean (km), only for coastal zones

prior    = 0.5 + 0.5 · normalised iNaturalist density (all years, research + needs-ID, Gaussian KDE σ = 8 km
           because fungi coordinates on iNat are obscured to ~0.2°). Unobserved-but-suitable habitat still scores.

dynamic  = rain^w · temp^w · phenology^w · trigger^w        (weights per species)
           rain       0.7·ramp(14-day precip, 25→80 mm) + 0.3·ramp(30-day precip, 50→160 mm); floor 0.1
           temp       trapezoid on 7-day mean air temp, lapse-adjusted (−6.5 °C/km) from the weather-grid
                      elevation to the cell elevation; ideal ≈ 7–15 °C, matsutake a bit colder; too-warm floor 0.1;
                      0 after ≥3 days below −2 °C in the last 14 (hard frost)
           phenology  Gaussian in day-of-year (peak / sigma per species; matsutake peaks earlier at elevation);
                      recent regional iNat sightings (last 21 d) raise its floor ("the flush has started")
           trigger    response to the season-opening rain event (≥15 mm in 48 h after Aug 15): rises over ~10 days,
                      plateaus, decays to a baseline. Chanterelles/matsutake key off the *first* event, kings off
                      the *latest* big rain (they flush fast and briefly).
```

The popup on each cell shows the pieces: score, habitat, history prior, "now" (dynamic), rain 14/30 d, 3-day
forecast, soil-temp proxy, days since the trigger, elevation, ecoregion class, tree cover.

## Repo layout

```
scripts/build_static.py   one-time static layers -> data/static/  (committed)
scripts/daily_update.py   daily: weather + recent iNat -> site/data/  (built by Actions, deployed, never committed)
scripts/model.py          ALL species parameters + scoring functions
scripts/common.py         grid definition, paths, polite HTTP
data/static/              elevation, tree cover, ecoregion class, coast distance, iNat KDE + raw cache, PAD-US
data/static/preview/      PNG eyeball checks of every static layer
site/index.html           single-file frontend (deck.gl + MapLibre, CDN, no build step)
site/data/                generated output (a seed copy is committed so the site works before the first Action run)
.github/workflows/daily.yml
```

### Data sources (all free, no keys)

* Weather — [Open-Meteo](https://open-meteo.com) forecast API, 0.1° grid (754 points), `past_days` up to 92 so the
  season-trigger lookback always reaches Aug 15, plus a 3-day forecast. Interpolated bilinearly to cells; temperature
  lapse-adjusted to cell elevation.
* Observations — [iNaturalist API](https://api.inaturalist.org/v1/docs/). Historical pull cached in
  `data/static/inat_historical.json`; daily run fetches only the last 30 days. Obscured coordinates are flagged.
* Elevation — [AWS Terrain Tiles](https://registry.opendata.aws/terrain-tiles/) (terrarium PNG, z10), also used for
  the hillshade in the browser.
* Tree cover — [ESA WorldCover 2021](https://esa-worldcover.org) 10 m, read from the COG overview (~75 m) with Pillow.
* Ecoregions — [EPA Level IV](https://www.epa.gov/eco-research/ecoregions) for OR and WA.
* Public land — [USGS PAD-US 4.1](https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-overview)
  fee layer: USFS, BLM, and state forest / state resource-management lands, simplified to ~150 m.
* Basemap — [OpenFreeMap](https://openfreemap.org) Positron style.

**Public land ≠ permission to pick.** Personal-use harvest rules vary: Mt Hood, Gifford Pinchot and Siuslaw National
Forests issue free-use permits with limits; WA DNR lands need a permit and a Discover Pass; Oregon state forests and
BLM districts have their own rules; national parks and most state parks prohibit harvest. The overlay shows where
picking is *plausible*, not legal advice. Respect private land and closures.

## Setup

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/build_static.py      # one time, ~3-5 min; downloads ~170 MB into cache/ (gitignored)
.venv/bin/python scripts/daily_update.py      # ~30 s; writes site/data/
python3 -m http.server 8765 --directory site  # open http://localhost:8765
```

Deploy:

1. Create a GitHub repo and push (`gh repo create mushroom-map --public --source=. --push`).
2. Settings → Pages → *Build and deployment* → Source: **GitHub Actions**.
3. Actions → *Daily forecast* → *Run workflow* once. The site is live at `https://<user>.github.io/mushroom-map/`
   and rebuilds at 06:30 Pacific daily. Daily data is deployed as a Pages artifact, never committed, so the repo
   stays small.

If the weather fetch fails the workflow exits non-zero (the previous deployment stays live). If the iNaturalist
recent-sightings fetch fails, the run continues with the committed seed `recent_obs.geojson` and a warning.

## Tuning

Edit `SPECIES` in `scripts/model.py` — elevation bands, ecoregion weights, rain thresholds, temperature band,
phenology peak/sigma, trigger curve, factor weights, colour. Habitat is composed at run time from the committed static
layers, so a tuning change only needs `daily_update.py` (no static rebuild). Re-run `build_static.py` only when the
*inputs* change (`--refresh-inat` to re-pull the historical observations, `--refresh-padus` for public land), or to add
a species (add its iNat taxon ids to `SPECIES`, then `--refresh-inat`).

To see how the model behaves under other conditions before the rain arrives (or to check a tuning change):

```bash
.venv/bin/python scripts/scenario.py --rain14 60 --rain30 120 --temp 14 --doy 288 --first 12 --latest 3 --recent 5 --png
```

prints the score distribution per species for synthetic conditions and writes `data/static/preview/scenario_*.png`.

## Frontend

`site/index.html` loads `meta.json`, gunzips `scores.bin.gz` (uint8 `[4, rows, cols]`) and `factors.bin.gz`
(uint8 `[rows, cols, 22]`, layout documented in `meta.json`) in the browser, paints a canvas and renders it as a
deck.gl `BitmapLayer` with nearest-neighbour filtering so cells stay crisp when zoomed. Species chips toggle any
combination (cell colour = highest-scoring species, opacity ∝ score); a slider hides cells below a score (default 15).
Public land (GeoJsonLayer), recent sightings (ScatterplotLayer, dashed 10 km halo for obscured records) and the
hillshade are toggleable. Daily payload ≈ 0.8 MB.
