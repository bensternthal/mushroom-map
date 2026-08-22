"""Shared grid definition, paths and polite HTTP helpers."""
import math
import os
import time

import numpy as np
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_STATIC = os.path.join(ROOT, "data", "static")
SITE_DATA = os.path.join(ROOT, "site", "data")
CACHE = os.path.join(ROOT, "cache")

# ---- Grid definition (EPSG:4326 lattice, row 0 = north edge) -----------------
LAT_MIN, LON_MIN, LON_MAX = 44.5, -124.2, -121.4
DLAT, DLON = 0.009, 0.0125                      # ~1 km at 45.7 N
NCOLS = int(round((LON_MAX - LON_MIN) / DLON))  # 224
NROWS = int(math.ceil((47.0 - LAT_MIN) / DLAT - 1e-9))  # 278
LAT_MAX = round(LAT_MIN + NROWS * DLAT, 6)      # 47.002
BBOX = [LON_MIN, LAT_MIN, LON_MAX, LAT_MAX]     # [W, S, E, N]
KM_PER_DEG_LAT = 111.32
CELL_KM_Y = DLAT * KM_PER_DEG_LAT
CELL_KM_X = DLON * KM_PER_DEG_LAT * math.cos(math.radians((LAT_MIN + LAT_MAX) / 2))

USER_AGENT = "pnw-mushroom-map/1.0 (static GitHub Pages project; ben@sternthal.org)"


def cell_centers():
    """1-D arrays of cell-center latitudes (north->south) and longitudes (west->east)."""
    lats = LAT_MAX - (np.arange(NROWS) + 0.5) * DLAT
    lons = LON_MIN + (np.arange(NCOLS) + 0.5) * DLON
    return lats, lons


def cell_index(lat, lon):
    """Vectorised (row, col) for lat/lon arrays; out-of-bbox -> -1."""
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    row = np.floor((LAT_MAX - lat) / DLAT).astype(np.int64)
    col = np.floor((lon - LON_MIN) / DLON).astype(np.int64)
    bad = (row < 0) | (row >= NROWS) | (col < 0) | (col >= NCOLS)
    row[bad] = -1
    col[bad] = -1
    return row, col


def grid_meta():
    return {
        "crs": "EPSG:4326",
        "bbox": BBOX,
        "rows": NROWS,
        "cols": NCOLS,
        "dlat": DLAT,
        "dlon": DLON,
        "row0": "north",
        "cell_km": [round(CELL_KM_X, 3), round(CELL_KM_Y, 3)],
    }


class FetchError(RuntimeError):
    pass


def _request(url, params=None, tries=5, timeout=90, backoff=2.0, stream=False):
    last = None
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout,
                             headers={"User-Agent": USER_AGENT}, stream=stream)
            if r.status_code == 200:
                return r
            if r.status_code == 429:                     # Open-Meteo minutely quota: wait it out
                last = f"HTTP 429: {r.text[:120]}"
                print(f"  rate limited, sleeping 65s ({attempt + 1}/{tries})")
                time.sleep(65)
                continue
            if r.status_code in (500, 502, 503, 504):
                last = f"HTTP {r.status_code}: {r.text[:200]}"
            else:
                raise FetchError(f"{url} -> HTTP {r.status_code}: {r.text[:300]}")
        except (requests.ConnectionError, requests.Timeout) as e:
            last = repr(e)
        wait = backoff ** attempt
        print(f"  retry {attempt + 1}/{tries} after {wait:.0f}s ({last})")
        time.sleep(wait)
    raise FetchError(f"{url} failed after {tries} tries: {last}")


def http_json(url, params=None, **kw):
    return _request(url, params, **kw).json()


def http_bytes(url, params=None, **kw):
    return _request(url, params, **kw).content


def save_npz(path, **arrays):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, **arrays)
