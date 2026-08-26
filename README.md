# PNW Mushroom Forecast Map

A static, daily-updated deck.gl map of northwest Oregon / southwest Washington (lat 44.5–47.0, lon −124.2…−121.4,
"coast to Cascades") on a ~1 km grid, colour-coded by the **relative likelihood** of finding the region's common edible (and two Psilocybe) species:

| Species | Key | Habitat encoded |
|---|---|---|
| Matsutake (*Tricholoma murrillianum*, incl. *T. magnivelare* records) | `matsutake` | coastal shore-pine/dune zone (≤8 km from the ocean, <250 m) **or** Cascade hemlock / true fir 600–1600 m |
| White chanterelle (*Cantharellus subalbidus*) | `white_chanterelle` | Doug-fir / western hemlock 0–1000 m |
| Pacific golden chanterelle (*Cantharellus formosus*) | `golden_chanterelle` | Doug-fir / western hemlock 0–1200 m, Coast Range + W Cascades |
| Fall king bolete (*Boletus edulis* group) | `king_bolete` | coastal Sitka spruce / shore pine (≤10 km from ocean) **or** Cascade true fir / hemlock 800–1500 m (weighted 0.7) |
| *Psilocybe azurescens* | `psilocybe_azurescens` | coastal dune grass / drift-wood debris, ≤6 km from the ocean, <120 m; no forest requirement; fruits through frost (Oct–Dec) |
| *Psilocybe cyanescens* | `psilocybe_cyanescens` | wood-chip / mulch landscaping (WorldCover built-up fraction) in lowland towns, secondarily lowland alder/cottonwood debris; <600 m; fruits through frost (Oct–Jan) |
| Morels, grouped (*M. elata, esculenta, deliciosa, americana, snyderi*) | `morchella` | **spring** (Apr–Jun, later upslope): conifer/mixed forest <1400 m, secondarily valley riparian / gallery forest and old orchards; fall rain-trigger weight 0; burn-year flushes not modelled |

The number shown is a **0–100 relative suitability index, not a probability** and not a guarantee. The two *Psilocybe*
layers are natural-history information (habitat model + public iNaturalist records, which iNat obscures for these
taxa); they are controlled substances whose legal status differs between Oregon and Washington and changes over time —
the map is not collection guidance. It is meant to
answer "where and when is it worth looking?", and to be tuned by someone who actually picks.

Live site: Cloudflare Workers (`mushroom-map.<subdomain>.workers.dev`), rebuilt every morning by GitHub Actions (see *Setup*).

## How a cell is scored

Everything tunable lives in one dict at the top of [`scripts/model.py`](scripts/model.py) (`SPECIES`) plus a few
globals right below it. Per cell and species:

```
score = 100 · (habitat · prior · dynamic) ^ GAMMA            GAMMA = 0.5

habitat  = forest_score · max over habitat zones ( ecoregion_weight · elevation_band · coast_band )
           forest_score   ramp on ESA WorldCover tree-cover fraction (0.15 → 0.5); none for non-forest species
           builtup        optional trapezoid on WorldCover built-up fraction (landscaping / wood-chip species)
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
scripts/build_season.py   season playback: archived weather -> site/seasons/<year>/  (built by Actions)
scripts/engine.py         scoring engine shared by the two above
scripts/model.py          ALL species parameters + scoring functions
scripts/common.py         grid definition, paths, polite HTTP
data/static/              elevation, tree cover, ecoregion class, coast distance, iNat KDE + raw cache, PAD-US
data/static/preview/      PNG eyeball checks of every static layer
site/index.html           single-file frontend (deck.gl + MapLibre, CDN, no build step)
site/data/                generated output (a seed copy is committed so the site works before the first Action run)
.github/workflows/daily.yml   build on Actions, deploy to Cloudflare
wrangler.toml             Cloudflare Worker (static assets) config
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

Deploy (private repo → Cloudflare Workers static assets, free):

GitHub Pages is not available on a private repo with a free GitHub plan, so the workflow builds on GitHub Actions
and uploads `site/` straight to **Cloudflare** as a Worker with static assets (`wrangler.toml`; direct upload, so
Cloudflare never needs access to the repo). One-time setup:

1. Cloudflare dashboard → *Workers & Pages* → *Create* → **Upload assets** (or *Start with Hello World*) → name
   `mushroom-map`. (The dashboard creates a Worker, URL `https://mushroom-map.<subdomain>.workers.dev`; the
   workflow overwrites it on every run.)
2. Cloudflare → *My Profile* → *API Tokens* → *Create Token* → template **Edit Cloudflare Workers**. Copy the token
   and your *Account ID* (account Overview page, right-hand sidebar).
3. Add both as repo secrets:
   ```bash
   gh secret set CLOUDFLARE_API_TOKEN -R bensternthal/mushroom-map
   gh secret set CLOUDFLARE_ACCOUNT_ID -R bensternthal/mushroom-map
   ```
4. Actions → *Daily forecast* → *Run workflow*. The site rebuilds at 06:30 Pacific daily. Daily data is uploaded as
   a deployment, never committed, so the repo stays small. Each run also keeps `site/` as a 7-day workflow artifact.
   Request metrics: Workers & Pages → mushroom-map → *Metrics*.

(To use GitHub Pages instead — public repo or a Pro plan — replace the deploy step with
`actions/upload-pages-artifact` + `actions/deploy-pages` and set Pages → Source to *GitHub Actions*.)

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

## Season playback

`scripts/build_season.py --year 2024 --year 2025 --year current` scores every day of a season (Mar 1 → Dec 31 —
spring morels through the fall flush — or up to yesterday for the current one) from Open-Meteo's ERA5 **archive** and the cached iNaturalist records, and
writes `site/seasons/<year>/<date>.bin.gz` (one score grid per day, ~50–100 KB) plus `timeline.json` (regional rain /
temperature / trigger coverage / sighting counts per day, and the season's observations for the pins). The workflow
rebuilds these every run (nothing is committed); the frontend's *Season playback* control loads days lazily as you
scrub or play, with a timeline strip (rain bars, max-score line, sightings stacked by species, trigger date dashed,
frost shaded). Playback uses reanalysis weather rather than the forecast model and iNat records as they exist today
(including later uploads), so it is for checking timing and broad pattern — and for tuning — not for grading single
cells. To add a season, add another `--year` in `daily.yml`.

## Frontend

`site/index.html` loads `meta.json`, gunzips `scores.bin.gz` (uint8 `[4, rows, cols]`) and `factors.bin.gz`
(uint8 `[rows, cols, 22]`, layout documented in `meta.json`) in the browser, paints a canvas and renders it as a
deck.gl `BitmapLayer` with nearest-neighbour filtering so cells stay crisp when zoomed. Species chips toggle any
combination (cell colour = highest-scoring species, opacity ∝ score); a slider hides cells below a score (default 15).
Public land (GeoJsonLayer), recent sightings (ScatterplotLayer, dashed 10 km halo for obscured records) and the
hillshade are toggleable. Daily payload ≈ 0.8 MB.
