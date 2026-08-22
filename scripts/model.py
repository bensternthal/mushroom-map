"""Scoring model for the PNW mushroom forecast.

EVERYTHING TUNABLE LIVES IN `SPECIES` (per species) AND THE GLOBALS BELOW IT.
Scores are a 0-100 *relative suitability index*, not a probability.

Per cell:  score = 100 * (habitat * prior * dynamic) ** GAMMA
  habitat  = forest_score * max over habitat zones (ecoregion weight * elevation band * coast band)
  prior    = 0.5 + 0.5 * normalised iNaturalist historical density (KDE, sigma ~8 km)
  dynamic  = rain^w_rain * temp^w_temp * phenology^w_pheno * trigger^w_trigger
"""
import numpy as np

# ---------------------------------------------------------------------------
# Ecoregion assemblage classes (EPA Level IV -> forest assemblage)
# ---------------------------------------------------------------------------
ECO_CLASSES = [
    "none",                 # 0 water / outside data
    "coastal_lowland",      # 1 Sitka spruce zone, dunes, shore pine (L4 1a)
    "coastal_upland",       # 2 spruce / western hemlock coastal hills (1b)
    "coast_range",          # 3 Coast Range western hemlock zone (1c-1h)
    "willamette_valley",    # 4 valley floor ag/urban (3a-3c)
    "valley_foothills",     # 5 valley margin Doug-fir / oak (3d)
    "puget_lowland",        # 6 Puget/Cowlitz lowlands (2a-2i)
    "west_cascades_low",    # 7 western hemlock zone, W Cascades lowlands (4a)
    "west_cascades_montane",# 8 Pacific silver fir zone (4b)
    "cascade_crest",        # 9 mountain hemlock / crest forest (4c, 4e)
    "subalpine",            # 10 subalpine/alpine (4d, 77c)
    "north_cascades",       # 11 North Cascades forests (77a, 77b)
    "east_cascades",        # 12 east-side pine/fir (9x)
    "other",                # 13 anything else (Columbia Plateau etc.)
]
ECO_LABELS = {
    "none": "water / no data",
    "coastal_lowland": "Coastal lowlands (Sitka spruce / dunes)",
    "coastal_upland": "Coastal uplands (spruce-hemlock)",
    "coast_range": "Coast Range (western hemlock zone)",
    "willamette_valley": "Willamette Valley floor",
    "valley_foothills": "Valley foothills (Doug-fir / oak)",
    "puget_lowland": "Puget / Cowlitz lowlands",
    "west_cascades_low": "W Cascades lowlands (western hemlock)",
    "west_cascades_montane": "W Cascades montane (silver fir)",
    "cascade_crest": "Cascade crest (mountain hemlock)",
    "subalpine": "Subalpine / alpine",
    "north_cascades": "North Cascades forest",
    "east_cascades": "East Cascades (pine / fir)",
    "other": "Other (non-forest region)",
}
ECO_L4_TO_CLASS = {
    "1a": "coastal_lowland", "1b": "coastal_upland",
    "1c": "coast_range", "1d": "coast_range", "1e": "coast_range", "1f": "coast_range",
    "1g": "coast_range", "1h": "coast_range", "1i": "coast_range",
    "3a": "willamette_valley", "3b": "willamette_valley", "3c": "willamette_valley",
    "3d": "valley_foothills",
    "4a": "west_cascades_low", "4b": "west_cascades_montane",
    "4c": "cascade_crest", "4e": "cascade_crest", "4f": "west_cascades_low",
    "4d": "subalpine", "77c": "subalpine",
    "77a": "north_cascades", "77b": "north_cascades", "77i": "subalpine",
}
for _c in "abcdefghi":
    ECO_L4_TO_CLASS.setdefault("2" + _c, "puget_lowland")
for _c in "abcdefghijklmnopqrstuvwxyz":
    ECO_L4_TO_CLASS.setdefault("9" + _c, "east_cascades")


def eco_class_id(l4code):
    return ECO_CLASSES.index(ECO_L4_TO_CLASS.get(l4code, "other"))


# ---------------------------------------------------------------------------
# Species parameters  (order here = species index in scores.bin)
# ---------------------------------------------------------------------------
# Conventions
#   elev       : trapezoid (a, b, c, d) in metres -> 0 below a, 1 between b..c, 0 above d
#   coast_km   : trapezoid on distance-to-ocean in km (omit = anywhere)
#   eco        : ecoregion class -> weight 0..1 (missing = 0)
#   forest     : (lo, hi) ramp on tree-cover fraction of the cell
#   rain14/30  : (lo, hi) mm ramp; rain30_weight blends the 30-day term in
#   temp       : trapezoid (a, b, c, d) on 7-day mean air temp (C), lapse-adjusted to cell
#   pheno      : day-of-year Gaussian; peak shifts `elev_shift_days_per_100m` for every
#                100 m above `elev_ref`;  recent iNat sightings raise the floor
#   trigger    : response to "season-opening" rain (>= TRIGGER_MM_48H mm in 48 h after Aug 15)
#                mode 'first' = first such event of the season, 'latest' = most recent one
#                pre = value before any trigger, rise/plateau/tail in days, post = value after tail
SPECIES = {
    "matsutake": dict(
        name="Matsutake", latin="Tricholoma murrillianum (incl. T. magnivelare)",
        inat_taxa=[521711, 62483], color=[240, 120, 30],
        zones=[
            dict(label="coastal shore pine / dunes",
                 eco={"coastal_lowland": 1.0, "coastal_upland": 0.4},
                 elev=(-1, 0, 100, 250), coast_km=(-1, 0, 3, 8)),
            dict(label="Cascade hemlock / true fir 600-1600 m",
                 eco={"west_cascades_low": 0.5, "west_cascades_montane": 1.0, "cascade_crest": 1.0,
                      "north_cascades": 1.0, "subalpine": 0.3, "east_cascades": 0.7,
                      "coast_range": 0.2, "other": 0.1},
                 elev=(450, 600, 1600, 1900)),
        ],
        forest=(0.15, 0.5),
        rain14=(25, 80), rain30=(50, 160), rain30_weight=0.3,
        temp=(1, 5, 13, 18), frost_kill=True,
        pheno=dict(peak_doy=283, sigma_days=20, elev_shift_days_per_100m=-1.2, elev_ref=600),
        trigger=dict(mode="first", pre=0.25, rise=10, plateau=16, tail=35, post=0.6),
        weights=dict(rain=1.0, temp=1.0, pheno=1.0, trigger=1.0),
    ),
    "white_chanterelle": dict(
        name="White chanterelle", latin="Cantharellus subalbidus",
        inat_taxa=[54132], color=[30, 170, 160],
        zones=[
            dict(label="Doug-fir / western hemlock 0-1000 m",
                 eco={"coastal_lowland": 0.5, "coastal_upland": 0.9, "coast_range": 1.0,
                      "valley_foothills": 0.6, "willamette_valley": 0.05, "puget_lowland": 0.5,
                      "west_cascades_low": 1.0, "west_cascades_montane": 0.6, "cascade_crest": 0.2,
                      "north_cascades": 0.6, "subalpine": 0.05, "east_cascades": 0.2, "other": 0.1},
                 elev=(-1, 0, 1000, 1300)),
        ],
        forest=(0.15, 0.5),
        rain14=(25, 80), rain30=(50, 160), rain30_weight=0.3,
        temp=(2, 7, 15, 19), frost_kill=True,
        pheno=dict(peak_doy=283, sigma_days=22, elev_shift_days_per_100m=0.0, elev_ref=0),
        trigger=dict(mode="first", pre=0.25, rise=10, plateau=18, tail=40, post=0.6),
        weights=dict(rain=1.0, temp=1.0, pheno=1.0, trigger=1.0),
    ),
    "golden_chanterelle": dict(
        name="Pacific golden chanterelle", latin="Cantharellus formosus",
        inat_taxa=[120443], color=[245, 190, 20],
        zones=[
            dict(label="Doug-fir / western hemlock 0-1200 m",
                 eco={"coastal_lowland": 0.8, "coastal_upland": 1.0, "coast_range": 1.0,
                      "valley_foothills": 0.5, "willamette_valley": 0.05, "puget_lowland": 0.5,
                      "west_cascades_low": 1.0, "west_cascades_montane": 0.8, "cascade_crest": 0.3,
                      "north_cascades": 0.8, "subalpine": 0.05, "east_cascades": 0.2, "other": 0.1},
                 elev=(-1, 0, 1200, 1500)),
        ],
        forest=(0.15, 0.5),
        rain14=(20, 70), rain30=(40, 140), rain30_weight=0.35,   # most rain-responsive
        temp=(2, 7, 15, 19), frost_kill=True,
        pheno=dict(peak_doy=293, sigma_days=30, elev_shift_days_per_100m=0.0, elev_ref=0),
        trigger=dict(mode="first", pre=0.3, rise=8, plateau=21, tail=50, post=0.7),  # long season
        weights=dict(rain=1.2, temp=1.0, pheno=1.0, trigger=0.8),
    ),
    "king_bolete": dict(
        name="Fall king bolete", latin="Boletus edulis group (incl. var. grandedulis)",
        inat_taxa=[48701], color=[170, 60, 40],
        zones=[
            dict(label="coastal Sitka spruce / shore pine",
                 eco={"coastal_lowland": 1.0, "coastal_upland": 0.7, "coast_range": 0.3},
                 elev=(-1, 0, 200, 400), coast_km=(-1, 0, 5, 10)),
            dict(label="Cascade true fir / hemlock 800-1500 m", weight=0.7,
                 eco={"west_cascades_low": 0.3, "west_cascades_montane": 0.8, "cascade_crest": 0.8,
                      "north_cascades": 0.8, "subalpine": 0.2, "east_cascades": 0.5, "other": 0.1},
                 elev=(650, 800, 1500, 1800)),
        ],
        forest=(0.15, 0.5),
        rain14=(25, 80), rain30=(50, 160), rain30_weight=0.2,
        temp=(2, 7, 15, 19), frost_kill=True,
        pheno=dict(peak_doy=278, sigma_days=18, elev_shift_days_per_100m=0.0, elev_ref=0),
        trigger=dict(mode="latest", pre=0.15, rise=7, plateau=14, tail=28, post=0.3),  # flushes fast after big rain
        weights=dict(rain=1.0, temp=1.0, pheno=1.0, trigger=1.3),
    ),
}
SPECIES_ORDER = list(SPECIES)

# ---- Global tunables ---------------------------------------------------------
GAMMA = 0.5                  # score = 100 * (habitat*prior*dynamic)**GAMMA; <1 spreads low values out
RAIN_FLOOR = 0.10            # bone-dry ground still scores this much (keeps habitat structure visible)
TEMP_HOT_FLOOR = 0.10        # too-warm soil penalty floor (heat delays fruiting, does not kill the mycelium)
PHENO_FLOOR = 0.05           # off-season floor on the phenology curve
LAPSE_C_PER_KM = 6.5         # air temp lapse rate used to move weather-grid temps to cell elevation
TRIGGER_MM_48H = 15.0        # season-opening rain event threshold
SEASON_START = (8, 15)       # (month, day) trigger events only count from here
FROST_TEMP_C = -2.0          # hard-frost threshold on lapse-adjusted daily minimum
FROST_DAYS = 3               # this many frost days in FROST_WINDOW -> season over (frost_kill species)
FROST_WINDOW = 14
KDE_SIGMA_KM = 8.0           # iNat coordinates for fungi are obscured ~0.2 deg; smooth accordingly
PRIOR_FLOOR = 0.5            # prior = floor + (1-floor) * normalised density
RECENT_DAYS = 21             # window for the "flush has started" sighting signal
RECENT_FLOOR_BASE, RECENT_FLOOR_PER_OBS, RECENT_FLOOR_CAP = 0.45, 0.08, 0.9
SCORE_MIN_DISPLAY = 15       # frontend hides cells below this


# ---------------------------------------------------------------------------
# Factor functions (all vectorised over numpy arrays)
# ---------------------------------------------------------------------------
def ramp(x, lo, hi):
    return np.clip((np.asarray(x, dtype=np.float32) - lo) / float(hi - lo), 0.0, 1.0)


def trapezoid(x, a, b, c, d):
    x = np.asarray(x, dtype=np.float32)
    up = np.ones_like(x) if b <= a else np.clip((x - a) / float(b - a), 0.0, 1.0)
    down = np.ones_like(x) if d <= c else np.clip((d - x) / float(d - c), 0.0, 1.0)
    return np.minimum(up, down)


def habitat(sp, elev, tree_frac, eco_class, coast_km):
    """Static 0-1 habitat suitability for one species."""
    forest = ramp(tree_frac, *sp["forest"])
    best = np.zeros(np.shape(elev), dtype=np.float32)
    for z in sp["zones"]:
        w = np.array([z["eco"].get(c, 0.0) for c in ECO_CLASSES], dtype=np.float32)[eco_class]
        v = w * trapezoid(elev, *z["elev"]) * float(z.get("weight", 1.0))
        if "coast_km" in z:
            v = v * trapezoid(coast_km, *z["coast_km"])
        best = np.maximum(best, v)
    return forest * best


def prior_from_density(kde_norm):
    return PRIOR_FLOOR + (1.0 - PRIOR_FLOOR) * np.asarray(kde_norm, dtype=np.float32)


def rain_score(sp, rain14, rain30):
    w = sp["rain30_weight"]
    r = (1 - w) * ramp(rain14, *sp["rain14"]) + w * ramp(rain30, *sp["rain30"])
    return np.maximum(r, RAIN_FLOOR)


def temp_score(sp, tmean7, frost_days):
    s = trapezoid(tmean7, *sp["temp"])
    hot = np.asarray(tmean7) > sp["temp"][2]
    s = np.where(hot, np.maximum(s, TEMP_HOT_FLOOR), s)     # too warm: floor; too cold: real zero
    if sp.get("frost_kill", True):
        s = np.where(np.asarray(frost_days) >= FROST_DAYS, 0.0, s)
    return s


def phenology(sp, doy, elev, recent_n):
    p = sp["pheno"]
    doy = float(doy)
    if doy < 150:            # Jan-May belongs to the previous fall season's tail
        doy += 365
    shift = p["elev_shift_days_per_100m"] * np.maximum(np.asarray(elev, dtype=np.float32) - p["elev_ref"], 0) / 100.0
    center = p["peak_doy"] + shift
    g = np.maximum(np.exp(-0.5 * ((doy - center) / p["sigma_days"]) ** 2), PHENO_FLOOR).astype(np.float32)
    if recent_n > 0:
        floor = min(RECENT_FLOOR_CAP, RECENT_FLOOR_BASE + RECENT_FLOOR_PER_OBS * recent_n)
        g = np.maximum(g, floor)
    return g


def trigger_score(sp, days_since_first, days_since_latest):
    t = sp["trigger"]
    d = np.asarray(days_since_first if t["mode"] == "first" else days_since_latest, dtype=np.float32)
    none = ~np.isfinite(d)
    d = np.where(none, 0, d)
    pre, post = t["pre"], t["post"]
    rise, plateau, tail = t["rise"], t["plateau"], t["tail"]
    s = np.where(d < rise, pre + (1 - pre) * d / rise, 1.0)
    s = np.where(d > plateau, 1.0 - (1 - post) * np.clip((d - plateau) / max(tail - plateau, 1), 0, 1), s)
    return np.where(none, pre, s).astype(np.float32)


def dynamic_score(sp, rain, temp, pheno, trig):
    w = sp["weights"]
    return (rain ** w["rain"]) * (temp ** w["temp"]) * (pheno ** w["pheno"]) * (trig ** w["trigger"])


def total_score(hab, prior, dyn):
    v = np.clip(hab * prior * dyn, 0, 1) ** GAMMA
    return np.round(100.0 * v).astype(np.uint8)
