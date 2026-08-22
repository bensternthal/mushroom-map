# Project Brief: PNW Mushroom Forecast Map

Paste this file into a new Claude Code session as the opening context. It captures all decisions already made — build exactly this unless something proves technically wrong, in which case explain and propose the closest alternative.

## Goal

A static website (GitHub Pages) showing a daily-updated, interactive deck.gl map of northwest Oregon and southwest Washington with a ~1 km grid, color-coded by relative likelihood of finding each of four mushroom species. Data refreshes daily via GitHub Actions. Owner: Ben (ben@sternthal.org), an experienced mushroom hunter.

## Region

Bounding box: **lat 44.5 → 47.0, lon -124.2 → -121.4** ("coast to Cascades": Clatsop/Tillamook coast, Coast Range, Willamette Valley foothills, Mt Hood NF, Gifford Pinchot NF, SW Washington coast). 1 km grid ≈ 55,000 cells. Grid definition: fixed row/col lattice in EPSG:4326, cell size 0.009° lat × 0.0125° lon (≈1 km at 45.7°N). Cell identity = (row, col) index so daily outputs stay aligned.

## Species and scoring factors

Output for each species is a 0–100 **relative suitability index** (be explicit in the UI that it's not a true probability). Compute per cell as a weighted product/sum of factor scores:

### Shared dynamic factors (from weather)
- **14-day and 30-day precipitation accumulation** (mm). Fall fruiting needs a soaking; score ramps up above ~25 mm/14d, saturates ~80 mm.
- **Days since first "season-opening" rain event** (≥15 mm in 48h after Aug 15). Most fall species flush 7–21 days after this trigger — model a lag curve peaking ~10–14 days post-trigger.
- **Soil temperature proxy**: 7-day mean air temp, lapse-rate adjusted to cell elevation (-6.5 °C/km from the weather grid point's elevation). Ideal band ~7–15 °C; penalize hard above 18 °C, and score 0 after sustained hard frost (min temp < -2 °C for 3+ days) except where noted.
- **Day-of-year phenology curve** per species (Gaussian-ish window).

### Per-species static habitat factors
1. **Matsutake** (*Tricholoma murrillianum* — the PNW matsutake; iNat has largely reclassified *T. magnivelare* to this in the West, query both). Habitat: coastal shore-pine/sandy-soil zones (dunes, Clatsop/Long Beach), and Cascade western hemlock / true fir / mountain hemlock at **600–1600 m**. Prefers poor, well-drained soils. Phenology: mid-Sept → mid-Nov, earlier at elevation. Slightly colder optimum than chanterelles.
2. **White chanterelle** (*Cantharellus subalbidus*). Doug-fir/western hemlock, **0–1000 m**, often slightly drier/better-drained sites than formosus. Sept–Nov.
3. **Pacific golden chanterelle** (*Cantharellus formosus*). Doug-fir/western hemlock throughout Coast Range and west Cascades, **0–1200 m**. Sept → Dec (until sustained frost); most rain-responsive of the four; long season, tolerates ongoing wet.
4. **Fall king bolete** (*Boletus edulis* group; on the coast often *B. edulis* var. *grandedulis*; query the *edulis* complex). Strongest on the **coastal Sitka spruce zone** (spruce/shore pine near dunes), secondarily true fir/hemlock in Cascades 800–1500 m. Flushes fast, 1–3 weeks after big rain; shorter window than chanterelles. Sept–Nov.

Habitat suitability per cell = combination of: forest cover, elevation band, distance-to-coast (spruce zone ≈ within ~8 km of coastline), and ecoregion class (see data sources). Encode per-species habitat 0–1 as a static layer.

### Observation prior
- **iNaturalist historical density**: all research-grade + needs-ID observations of each taxon in the bbox (all years), kernel-density smoothed (σ ≈ 5–10 km — coordinates for fungi are often obscured ~0.2°, so do NOT trust them at 1 km; the KDE bandwidth must reflect that). Use as a multiplicative prior (e.g., 0.5 + 0.5·normalized density) so unobserved-but-suitable habitat still scores.
- **Recent activity signal**: observations of each species in the last 21 days anywhere in the region bump the regional phenology term (evidence the flush has started).

## Data sources (all free, no API keys)

- **Weather — Open-Meteo** (`api.open-meteo.com/v1/forecast` with `past_days`, and `archive-api.open-meteo.com` for longer history). Supports batched multi-point requests (comma-separated latitude/longitude lists). Fetch on a **0.1° coarse grid (~750 points)** — native model resolution is coarser than 1 km anyway — then interpolate to cells with elevation lapse-rate adjustment for temperature and nearest/bilinear for precip. Daily variables: `precipitation_sum, temperature_2m_max, temperature_2m_min`. Be polite: chunk requests, retry with backoff.
- **iNaturalist API** (`api.inaturalist.org/v1/observations`). Resolve taxon IDs first via `/v1/taxa?q=`. Paginate `per_page=200`, respect ~1 req/s, cache the historical pull as a committed file and only fetch a rolling recent window daily. Include `geoprivacy`/`taxon_geoprivacy` fields; note obscured coords in popups.
- **Elevation — AWS Terrain Tiles** (`s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png`, no auth). Zoom 10–11, decode `(R*256 + G + B/256) - 32768`. Build the 1 km elevation grid **once**, commit as a small compressed .npz.
- **Habitat/tree proxy**: EPA **Level IV ecoregions** (downloadable shapefile/GeoJSON) mapped to assemblages (Sitka spruce zone, coastal dunes, western hemlock zone, silver fir / mountain hemlock zone, Willamette Valley — score valley ag land ~0) + a forest mask from **ESA WorldCover 10 m** (public S3/Zenodo) or, simpler, tree-cover from terrain of OSM landuse if WorldCover proves painful. This is a one-time static layer; commit it. Don't over-engineer: ecoregion + elevation + forest mask is enough.
- **Public land — PAD-US** (USGS ArcGIS REST export for the bbox, one-time). Keep federal (USFS, BLM) and state forest polygons, simplify aggressively (shapely `simplify`), commit as GeoJSON. Note in README: personal-use harvest permits vary (Mt Hood NF, Gifford Pinchot NF, Siuslaw NF free-use permits; WA DNR needs a permit + Discover Pass) — the overlay shows where picking is *plausible*, not legal advice.

## Pipeline architecture

Python (`scripts/`), minimal deps (requests, numpy, shapely, pillow). Two stages:

1. **`build_static.py`** — run once (locally, committed outputs): elevation grid, habitat layers per species, iNat historical KDE, public-land GeoJSON, grid metadata. Outputs to `data/static/`.
2. **`daily_update.py`** — run by Actions: fetch 30-day weather + 3-day forecast for coarse grid, fetch recent iNat obs, compute the four score grids, write compact outputs to `site/data/`:
   - `meta.json` (grid def, generated timestamp, species list, score legend)
   - `scores.bin` — Uint8, shape [4 species × rows × cols], row-major (≈220 KB — tiny)
   - `recent_obs.geojson` — last 30 days iNat pins with species, date, obscured flag
   - copy static `public_lands.geojson`

Also emit a per-cell **explanation** without bloating data: the frontend recomputes nothing; instead include compact per-cell factor bytes in a second binary `factors.bin` (rain14 score, temp score, habitat score, prior score as Uint8 each) so the popup can show "why this score."

## Frontend (`site/index.html`, single file, CDN imports, no build step)

- **deck.gl + MapLibre GL JS**. Terrain-flavored basemap: free vector style (OpenFreeMap/Carto) + MapLibre **hillshade** from AWS terrarium raster-DEM tiles.
- Species **toggle chips** (any combination; when multiple are on, render the max or side-by-side opacity — pick max with per-species hue). Color ramp per species (e.g., matsutake orange, golden chanterelle gold, white chanterelle teal, king bolete brown-red), opacity ∝ score, cells below score 15 fully transparent.
- Render the grid as a deck.gl **GridCellLayer/BitmapLayer** (BitmapLayer over the bbox from a canvas-generated PNG is the fastest for 55k cells; GridCellLayer is fine too — pick whichever performs).
- **Click popup** per cell: score per active species, elevation, habitat class, 14/30-day rain, soil-temp proxy, historical-observation density, "days since season trigger."
- Toggleable **public land** overlay (GeoJsonLayer, subtle fill + outline, hover shows unit name).
- Toggleable **recent iNat sightings** pins (IconLayer/ScatterplotLayer), popup links to the observation; obscured ones drawn with a dashed halo + "location obscured ~20 km" note.
- Header: last-updated timestamp, small legend, disclaimer (relative index, not a guarantee; respect permits/private land).
- Must work zoomed from region-wide down to individual cells; keep total daily payload < 2 MB.

## GitHub Actions (`.github/workflows/daily.yml`)

- `schedule: cron "30 13 * * *"` (≈6:30 am Pacific) + `workflow_dispatch`.
- Steps: checkout → setup-python → pip install → `python scripts/daily_update.py` → upload `site/` with `actions/upload-pages-artifact` → `actions/deploy-pages`. **Deploy as a Pages artifact rather than committing daily data** to avoid repo bloat. Static data files are committed once by `build_static.py`.
- Fail loudly (non-zero exit) if any fetch fully fails; tolerate partial iNat failure by reusing last cached recent_obs.

## Repo layout

```
mushroom-map/
  README.md            # what it is, how scoring works, setup: create repo, push,
                       # Settings→Pages→GitHub Actions, run workflow manually once
  scripts/build_static.py
  scripts/daily_update.py
  scripts/model.py     # scoring functions, per-species params in a dict at top (tunable)
  data/static/         # committed one-time layers
  site/index.html
  site/data/           # generated daily (gitignored except a seed copy)
  .github/workflows/daily.yml
  requirements.txt
```

## Build order & verification

1. Verify all endpoints with real requests before writing the full pipeline.
2. `build_static.py`, run it, sanity-check layers (print stats, render a quick matplotlib PNG of each layer and eyeball: coast spruce zone hugging the coast, elevation bands correct, Willamette Valley scored low).
3. `model.py` + `daily_update.py`, run for real, check score distributions per species (shouldn't be all-0 or all-100; in late August scores should be low-to-moderate — season is just starting).
4. Frontend against the real generated data; test toggles, popup, zoom.
5. Workflow; then `gh repo create`, push, enable Pages, trigger a manual run, confirm the live site.

## Tuning expectations

Ben will iterate ("kings overpredicted on the coast," "add spring porcini," "widen the matsutake elevation band"). Keep all species parameters (elevation bands, rain thresholds, phenology windows, weights) in one obvious dict in `model.py` with comments, so tuning is a one-line change.
