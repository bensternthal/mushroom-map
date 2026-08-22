"""Scoring engine shared by daily_update.py (live) and build_season.py (playback).

Pipeline for one target day:
  weather grids [nd, 26, 29]  --bilinear-->  cell grids  --lapse-->  conditions(day)  --model-->  scores
"""
import datetime as dt
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C                      # noqa: E402
import model as M                       # noqa: E402

WX_STEP = 0.1
WX_LATS = np.round(np.arange(C.LAT_MIN, 47.0 + 1e-9, WX_STEP), 4)      # 26
WX_LONS = np.round(np.arange(C.LON_MIN, C.LON_MAX + 1e-9, WX_STEP), 4)  # 29
WX_POINTS = [(la, lo) for la in WX_LATS for lo in WX_LONS]              # 754
NS = len(M.SPECIES_ORDER)


def season_year(day):
    return day.year if day.month >= 8 else day.year - 1


class Static:
    """Committed static layers (data/static)."""

    def __init__(self):
        L = np.load(os.path.join(C.DATA_STATIC, "layers.npz"))
        self.elev = np.maximum(L["elev"].astype(np.float32), 0)
        self.tree_pct = L["tree_pct"]
        self.tree_frac = self.tree_pct.astype(np.float32) / 100
        self.eco_class = L["eco_class"]
        self.coast_km = L["coast_km"].astype(np.float32)
        self.ocean = L["ocean"]
        self.kde = np.load(os.path.join(C.DATA_STATIC, "inat_kde.npz"))["kde"].astype(np.float32) / 255
        self.habitat = np.stack([M.habitat(M.SPECIES[k], self.elev, self.tree_frac, self.eco_class, self.coast_km)
                                 for k in M.SPECIES_ORDER])
        self.prior = np.stack([M.prior_from_density(self.kde[i]) for i in range(NS)])


def parse_weather(results, expected_n, wx_lats=WX_LATS, wx_lons=WX_LONS):
    """Open-Meteo multi-point JSON -> (dates, P, TX, TN, E) with shapes [nd,nlat,nlon] / [nlat,nlon]."""
    if isinstance(results, dict):
        results = [results]
    if len(results) != expected_n:
        raise C.FetchError(f"weather: expected {expected_n} points, got {len(results)}")
    times = results[0]["daily"]["time"]
    nd = len(times)
    shape = (nd, len(wx_lats), len(wx_lons))
    P, TX, TN = (np.zeros(shape, np.float32) for _ in range(3))
    E = np.zeros(shape[1:], np.float32)
    for k, r in enumerate(results):
        i, j = divmod(k, len(wx_lons))
        dly = r["daily"]
        if dly["time"] != times:
            raise C.FetchError("weather: inconsistent date axis across points")
        for A, key in ((P, "precipitation_sum"), (TX, "temperature_2m_max"), (TN, "temperature_2m_min")):
            A[:, i, j] = np.array([v if v is not None else np.nan for v in dly[key]], np.float32)
        E[i, j] = r.get("elevation", 0.0) or 0.0
    for A in (P, TX, TN):                       # forward-fill the occasional null, then zero
        for t in range(1, nd):
            m = np.isnan(A[t])
            A[t][m] = A[t - 1][m]
        A[np.isnan(A)] = 0
    return [dt.date.fromisoformat(t) for t in times], P, TX, TN, E


def bilinear_to_cells(field, wx_lats=WX_LATS, wx_lons=WX_LONS):
    """field[..., len(wx_lats), len(wx_lons)] on a regular weather grid -> [..., NROWS, NCOLS]."""
    lats, lons = C.cell_centers()
    dlat, dlon = wx_lats[1] - wx_lats[0], wx_lons[1] - wx_lons[0]
    fi = np.clip((lats - wx_lats[0]) / dlat, 0, len(wx_lats) - 1 - 1e-6)
    fj = np.clip((lons - wx_lons[0]) / dlon, 0, len(wx_lons) - 1 - 1e-6)
    i0, j0 = np.floor(fi).astype(int), np.floor(fj).astype(int)
    wi, wj = (fi - i0).astype(np.float32)[:, None], (fj - j0).astype(np.float32)[None, :]
    f = np.asarray(field, np.float32)
    a = f[..., i0, :][..., :, j0]
    b = f[..., i0, :][..., :, j0 + 1]
    c = f[..., i0 + 1, :][..., :, j0]
    d = f[..., i0 + 1, :][..., :, j0 + 1]
    return a * (1 - wi) * (1 - wj) + b * (1 - wi) * wj + c * wi * (1 - wj) + d * wi * wj


class CellWeather:
    """Daily weather on the cell grid, temperatures lapse-adjusted to cell elevation."""

    def __init__(self, static, dates, P, TX, TN, E, wx_lats=WX_LATS, wx_lons=WX_LONS):
        self.dates = dates
        g = (wx_lats, wx_lons)
        self.P = bilinear_to_cells(P, *g)
        lapse = M.LAPSE_C_PER_KM * (static.elev - bilinear_to_cells(E, *g)) / 1000.0
        self.Tmean = bilinear_to_cells((TX + TN) / 2, *g) - lapse
        self.Tmin = bilinear_to_cells(TN, *g) - lapse
        p48 = self.P.copy()
        p48[1:] += self.P[:-1]
        self.event = p48 >= M.TRIGGER_MM_48H          # [nd, rows, cols] 48-h trigger events

    def conditions(self, day):
        """Conditions as of `day` using complete days strictly before it."""
        it = self.dates.index(day)
        if it < 30:
            raise ValueError(f"need 30 days of history before {day}")
        ss = dt.date(season_year(day), *M.SEASON_START)
        first = np.full(self.P.shape[1:], np.nan, np.float32)
        latest = np.full(self.P.shape[1:], np.nan, np.float32)
        for t in range(it):
            if self.dates[t] < ss:
                continue
            ev = self.event[t]
            ds = float((day - self.dates[t]).days)
            first[ev & np.isnan(first)] = ds
            latest[ev] = ds
        return dict(
            rain14=self.P[it - 14:it].sum(axis=0), rain30=self.P[it - 30:it].sum(axis=0),
            fc3=self.P[it:it + 3].sum(axis=0) if it + 1 <= len(self.dates) else np.zeros_like(first),
            tmean7=self.Tmean[it - 7:it].mean(axis=0),
            frost_days=(self.Tmin[it - M.FROST_WINDOW:it] < M.FROST_TEMP_C).sum(axis=0),
            ds_first=first, ds_latest=latest,
        )


def score_day(static, cond, day, recent_n):
    """Returns dict with scores [NS,rows,cols] uint8 plus factor grids."""
    doy = day.timetuple().tm_yday
    out = dict(scores=np.zeros((NS, C.NROWS, C.NCOLS), np.uint8), dyn=np.zeros((NS, C.NROWS, C.NCOLS), np.float32),
               factors={})
    for si, key in enumerate(M.SPECIES_ORDER):
        sp = M.SPECIES[key]
        rain = M.rain_score(sp, cond["rain14"], cond["rain30"])
        temp = M.temp_score(sp, cond["tmean7"], cond["frost_days"])
        pheno = M.phenology(sp, doy, static.elev, recent_n.get(key, 0))
        trig = M.trigger_score(sp, cond["ds_first"], cond["ds_latest"])
        dyn = M.dynamic_score(sp, rain, temp, pheno, trig)
        s = M.total_score(static.habitat[si], static.prior[si], dyn)
        s[static.ocean] = 0
        out["scores"][si] = s
        out["dyn"][si] = dyn
        out["factors"][key] = dict(rain=rain, temp=temp, pheno=pheno, trig=trig)
    return out


def species_summary(static, res):
    summ = {}
    for si, key in enumerate(M.SPECIES_ORDER):
        s = res["scores"][si]
        land = ~static.ocean & (static.habitat[si] > 0.05)
        f = res["factors"][key]
        summ[key] = {
            "score_max": int(s.max()), "score_mean_habitat": round(float(s[land].mean()), 1),
            "cells_ge15": int((s >= 15).sum()), "cells_ge50": int((s >= 50).sum()),
            "factor_means_habitat": {k: round(float(np.mean(np.asarray(v)[land]) if np.ndim(v) else float(v)), 2)
                                     for k, v in f.items()},
        }
    return summ


def count_recent(obs_dates_by_species, day, window=M.RECENT_DAYS):
    """obs_dates_by_species: {key: sorted np.array of date ordinals}."""
    lo, hi = day.toordinal() - window, day.toordinal()
    return {k: int(np.searchsorted(v, hi, side="right") - np.searchsorted(v, lo, side="left"))
            for k, v in obs_dates_by_species.items()}
