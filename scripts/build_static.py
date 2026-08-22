"""One-time static layer build.  Run locally; outputs in data/static/ are committed.

  python scripts/build_static.py [--refresh-inat] [--refresh-padus] [--no-preview]

Outputs
  data/static/layers.npz          elev (int16 m), tree_pct (uint8), eco_class (uint8),
                                  eco_l4 (uint8 index into grid.json l4_codes), coast_km (float16),
                                  ocean (bool)
  data/static/inat_kde.npz        kde (uint8 0-255, [species, rows, cols]) normalised density
  data/static/inat_historical.json  cached raw observations (id, date, lat, lon, obscured, grade)
  data/static/public_lands.geojson  simplified PAD-US federal + state forest polygons
  data/static/grid.json           grid definition + class tables + source notes
  data/static/preview/*.png       eyeball checks
"""
import argparse
import datetime as dt
import io
import json
import math
import os
import struct
import sys
import time
import zipfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C                      # noqa: E402
import model as M                       # noqa: E402

PREVIEW = os.path.join(C.DATA_STATIC, "preview")
N = C.NROWS * C.NCOLS


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def accumulate(row_of_y, col_of_x, values, weights_fn=None):
    """Average `values[H, W]` into grid cells given per-row cell rows and per-col cell cols.
    Returns (sum, count) flattened over the grid."""
    vr = row_of_y >= 0
    vc = col_of_x >= 0
    sub = values[vr][:, vc]
    idx = (row_of_y[vr][:, None] * C.NCOLS + col_of_x[vc][None, :]).ravel()
    s = np.bincount(idx, weights=sub.ravel().astype(np.float64), minlength=N)
    n = np.bincount(idx, minlength=N)
    return s, n


def stats(name, a):
    a = np.asarray(a, dtype=np.float64)
    print(f"  {name:12s} min {np.nanmin(a):8.2f}  p5 {np.nanpercentile(a, 5):8.2f}  "
          f"median {np.nanmedian(a):8.2f}  p95 {np.nanpercentile(a, 95):8.2f}  max {np.nanmax(a):8.2f}")


# ---------------------------------------------------------------------------
# 1. Elevation from AWS terrarium tiles (zoom 10)
# ---------------------------------------------------------------------------
def build_elevation(z=10):
    from PIL import Image
    n = 2 ** z

    def xt(lon):
        return int((lon + 180.0) / 360.0 * n)

    def yt(lat):
        r = math.radians(lat)
        return int((1.0 - math.log(math.tan(r) + 1 / math.cos(r)) / math.pi) / 2.0 * n)

    x0, x1 = xt(C.LON_MIN), xt(C.LON_MAX)
    y0, y1 = yt(C.LAT_MAX), yt(C.LAT_MIN)
    W, H = (x1 - x0 + 1) * 256, (y1 - y0 + 1) * 256
    mosaic = np.zeros((H, W), dtype=np.float32)
    tiles = [(x, y) for y in range(y0, y1 + 1) for x in range(x0, x1 + 1)]
    print(f"[elevation] {len(tiles)} terrarium tiles at z{z}")
    for i, (x, y) in enumerate(tiles):
        p = os.path.join(C.CACHE, "terrarium", str(z), str(x), f"{y}.png")
        if not os.path.exists(p):
            os.makedirs(os.path.dirname(p), exist_ok=True)
            data = C.http_bytes(f"https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png")
            with open(p, "wb") as f:
                f.write(data)
            time.sleep(0.05)
        im = np.asarray(Image.open(p).convert("RGB"), dtype=np.float32)
        e = im[..., 0] * 256.0 + im[..., 1] + im[..., 2] / 256.0 - 32768.0
        mosaic[(y - y0) * 256:(y - y0 + 1) * 256, (x - x0) * 256:(x - x0 + 1) * 256] = e
    px = x0 * 256 + np.arange(W) + 0.5
    py = y0 * 256 + np.arange(H) + 0.5
    lon = px / (256 * n) * 360.0 - 180.0
    lat = np.degrees(np.arctan(np.sinh(math.pi * (1 - 2 * py / (256 * n)))))
    row_of_y, _ = C.cell_index(lat, np.full_like(lat, C.LON_MIN + 1e-6))
    _, col_of_x = C.cell_index(np.full_like(lon, C.LAT_MIN + 1e-6), lon)
    s, cnt = accumulate(row_of_y, col_of_x, mosaic)
    elev = (s / np.maximum(cnt, 1)).reshape(C.NROWS, C.NCOLS)
    print(f"  pixels per cell ~{cnt[cnt > 0].mean():.0f}")
    stats("elev m", elev)
    return elev.astype(np.float32)


# ---------------------------------------------------------------------------
# 2. Forest fraction + ocean mask from ESA WorldCover 2021 (COG overview via Pillow)
# ---------------------------------------------------------------------------
WC_TILES = ["N42W126", "N42W123", "N45W126", "N45W123"]
WC_URL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{t}_Map.tif"


def build_worldcover(overview=3):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    tree = np.zeros(N)
    nodata = np.zeros(N)
    water = np.zeros(N)
    total = np.zeros(N)
    for t in WC_TILES:
        p = os.path.join(C.CACHE, f"ESA_WorldCover_10m_2021_v200_{t}_Map.tif")
        if not os.path.exists(p):
            print(f"[worldcover] downloading {t} (tens of MB, one time)")
            os.makedirs(C.CACHE, exist_ok=True)
            r = C._request(WC_URL.format(t=t), stream=True, timeout=600)
            with open(p, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        im = Image.open(p)
        im.seek(overview)
        a = np.asarray(im)
        size = a.shape[0]
        lat0 = int(t[1:3]) + 3
        lon0 = -int(t[4:7])
        step = 3.0 / size
        lat = lat0 - (np.arange(size) + 0.5) * step
        lon = lon0 + (np.arange(size) + 0.5) * step
        row_of_y, _ = C.cell_index(lat, np.full_like(lat, C.LON_MIN + 1e-6))
        _, col_of_x = C.cell_index(np.full_like(lon, C.LAT_MIN + 1e-6), lon)
        if not (row_of_y >= 0).any() or not (col_of_x >= 0).any():
            continue
        print(f"[worldcover] {t} overview {overview}: {size}x{size} px ({step * 111000:.0f} m)")
        for arr, cond in ((tree, a == 10), (nodata, a == 0), (water, a == 80)):
            s, _ = accumulate(row_of_y, col_of_x, cond.astype(np.float32))
            arr += s
        _, cnt = accumulate(row_of_y, col_of_x, a)
        total += cnt
    total = np.maximum(total, 1)
    tree_frac = (tree / total).reshape(C.NROWS, C.NCOLS)
    ocean_frac = (nodata / total).reshape(C.NROWS, C.NCOLS)
    water_frac = ((nodata + water) / total).reshape(C.NROWS, C.NCOLS)
    stats("tree_frac", tree_frac)
    print(f"  cells with >50% ocean(nodata): {(ocean_frac > 0.5).sum()}")
    return tree_frac.astype(np.float32), ocean_frac, water_frac


def shift(a, dr, dc):
    out = np.zeros_like(a)
    r0, r1 = max(dr, 0), a.shape[0] + min(dr, 0)
    c0, c1 = max(dc, 0), a.shape[1] + min(dc, 0)
    out[r0:r1, c0:c1] = a[r0 - dr:r1 - dr, c0 - dc:c1 - dc]
    return out


def build_coast(ocean_frac, water_frac):
    """Ocean = water connected to the bbox west edge (incl. bays, lower Columbia estuary) west of -123.4.
    Inside this bbox WorldCover codes the nearshore ocean as class 80 (water); nodata (0) is only far offshore."""
    lats, lons = C.cell_centers()
    water = water_frac > 0.5
    seed = (ocean_frac > 0.5) | (water & (np.arange(C.NCOLS)[None, :] == 0))   # open water on the west edge
    fillable = seed | (water & (lons[None, :] < -123.4))
    cur = seed.copy()
    for _ in range(2000):
        nxt = (cur | shift(cur, 1, 0) | shift(cur, -1, 0) | shift(cur, 0, 1) | shift(cur, 0, -1)) & fillable
        if np.array_equal(nxt, cur):
            break
        cur = nxt
    ocean = cur
    rows, cols = np.indices(ocean.shape)
    oy = (rows[ocean] + 0.5) * C.CELL_KM_Y
    ox = (cols[ocean] + 0.5) * C.CELL_KM_X
    ly = (rows.ravel() + 0.5) * C.CELL_KM_Y
    lx = (cols.ravel() + 0.5) * C.CELL_KM_X
    dist = np.full(N, 60.0, dtype=np.float32)
    step = 1500
    for i in range(0, N, step):
        d2 = (ly[i:i + step, None] - oy[None, :]) ** 2 + (lx[i:i + step, None] - ox[None, :]) ** 2
        dist[i:i + step] = np.sqrt(d2.min(axis=1))
    dist = np.minimum(dist, 60.0).reshape(C.NROWS, C.NCOLS)
    dist[ocean] = 0.0
    print(f"[coast] ocean cells {ocean.sum()}, land cells within 8 km of coast: {((dist <= 8) & ~ocean).sum()}")
    return ocean, dist


# ---------------------------------------------------------------------------
# 3. EPA Level IV ecoregions (shapefile in USGS Albers -> project cell centres forward)
# ---------------------------------------------------------------------------
ECO_URL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/{s}/{s}_eco_l4.zip"


def albers_forward(lat, lon, a=6378137.0, f=1 / 298.257222101,
                   lat1=29.5, lat2=45.5, lat0=23.0, lon0=-96.0):
    e2 = 2 * f - f * f
    e = math.sqrt(e2)

    def m(p):
        return np.cos(p) / np.sqrt(1 - e2 * np.sin(p) ** 2)

    def q(p):
        s = np.sin(p)
        return (1 - e2) * (s / (1 - e2 * s * s) - (1 / (2 * e)) * np.log((1 - e * s) / (1 + e * s)))

    p1, p2, p0 = map(math.radians, (lat1, lat2, lat0))
    m1, m2 = m(p1), m(p2)
    q0, q1, q2 = q(p0), q(p1), q(p2)
    n = (m1 * m1 - m2 * m2) / (q2 - q1)
    Cc = m1 * m1 + n * q1
    rho0 = a * math.sqrt(Cc - n * q0) / n
    phi = np.radians(np.asarray(lat, dtype=np.float64))
    lam = np.radians(np.asarray(lon, dtype=np.float64))
    rho = a * np.sqrt(Cc - n * q(phi)) / n
    th = n * (lam - math.radians(lon0))
    return rho * np.sin(th), rho0 - rho * np.cos(th)


def read_dbf(data):
    nrec, hlen, rlen = struct.unpack("<IHH", data[4:12])
    fields, p = [], 32
    while data[p] != 0x0D:
        b = data[p:p + 32]
        fields.append((b[:11].split(b"\0")[0].decode(), b[16]))
        p += 32
    out = []
    for i in range(nrec):
        r = data[hlen + i * rlen + 1: hlen + (i + 1) * rlen]
        rec, p = {}, 0
        for name, ln in fields:
            rec[name] = r[p:p + ln].decode("latin1").strip()
            p += ln
        out.append(rec)
    return out


def read_shp_rings(data):
    """Yield list-of-rings (Nx2 arrays) per record for polygon shapefiles."""
    pos, out = 100, []
    while pos < len(data):
        _, clen = struct.unpack(">ii", data[pos:pos + 8])
        pos += 8
        end = pos + clen * 2
        stype = struct.unpack("<i", data[pos:pos + 4])[0]
        rings = []
        if stype in (5, 15, 25):
            nparts, npts = struct.unpack("<ii", data[pos + 36:pos + 44])
            parts = struct.unpack(f"<{nparts}i", data[pos + 44:pos + 44 + 4 * nparts])
            pts = np.frombuffer(data, dtype="<f8", count=npts * 2, offset=pos + 44 + 4 * nparts).reshape(npts, 2)
            b = list(parts) + [npts]
            rings = [pts[b[i]:b[i + 1]] for i in range(nparts)]
        out.append(rings)
        pos = end
    return out


def build_ecoregions():
    import shapely
    from shapely import STRtree
    from shapely.geometry import Polygon
    polys, attrs = [], []
    for s in ("or", "wa"):
        p = os.path.join(C.CACHE, f"{s}_eco_l4.zip")
        if not os.path.exists(p):
            os.makedirs(C.CACHE, exist_ok=True)
            with open(p, "wb") as f:
                f.write(C.http_bytes(ECO_URL.format(s=s)))
        z = zipfile.ZipFile(p)
        shp = z.read(f"{s}_eco_l4.shp")
        dbf = read_dbf(z.read(f"{s}_eco_l4.dbf"))
        recs = read_shp_rings(shp)
        assert len(recs) == len(dbf), (len(recs), len(dbf))
        for rings, att in zip(recs, dbf):
            for ring in rings:
                if len(ring) < 4:
                    continue
                g = Polygon(ring)
                if not g.is_valid:
                    g = g.buffer(0)
                if g.is_empty:
                    continue
                polys.append(g)
                attrs.append(att)
    print(f"[ecoregions] {len(polys)} rings from OR+WA shapefiles")
    lats, lons = C.cell_centers()
    LAT, LON = np.meshgrid(lats, lons, indexing="ij")
    x, y = albers_forward(LAT.ravel(), LON.ravel())
    pts = shapely.points(x, y)
    tree = STRtree(polys)
    pi, gi = tree.query(pts, predicate="within")
    areas = np.array([g.area for g in polys])
    # nested rings (holes are separate outer rings of the enclosed polygon): keep the smallest match
    best = np.full(N, -1, dtype=np.int64)
    best_area = np.full(N, np.inf)
    for p_, g_ in zip(pi, gi):
        if areas[g_] < best_area[p_]:
            best_area[p_] = areas[g_]
            best[p_] = g_
    l4_codes = sorted({a["US_L4CODE"] for a in attrs})
    l4_names = {a["US_L4CODE"]: a["US_L4NAME"] for a in attrs}
    eco_l4 = np.zeros(N, dtype=np.uint8)          # 0 = none
    eco_class = np.zeros(N, dtype=np.uint8)
    hit = best >= 0
    codes = np.array([attrs[b]["US_L4CODE"] if b >= 0 else "" for b in best])
    for i, code in enumerate(l4_codes):
        sel = codes == code
        eco_l4[sel] = i + 1
        eco_class[sel] = M.eco_class_id(code)
    print(f"  cells matched {hit.sum()} / {N}")
    present = {}
    for code in codes[hit]:
        present[code] = present.get(code, 0) + 1
    for code, cnt in sorted(present.items(), key=lambda kv: -kv[1]):
        print(f"   {code:4s} {l4_names[code][:45]:45s} {cnt:6d} cells -> {M.ECO_L4_TO_CLASS.get(code, 'other')}")
    return eco_class.reshape(C.NROWS, C.NCOLS), eco_l4.reshape(C.NROWS, C.NCOLS), \
        [""] + l4_codes, [""] + [l4_names[c] for c in l4_codes]


# ---------------------------------------------------------------------------
# 4. iNaturalist historical observations + KDE
# ---------------------------------------------------------------------------
INAT = "https://api.inaturalist.org/v1/observations"
HIST_PATH = os.path.join(C.DATA_STATIC, "inat_historical.json")


def fetch_inat_all(taxa):
    out, id_above = [], 0
    while True:
        d = C.http_json(INAT, dict(
            taxon_id=",".join(map(str, taxa)), nelat=C.LAT_MAX, nelng=C.LON_MAX, swlat=C.LAT_MIN, swlng=C.LON_MIN,
            quality_grade="research,needs_id", per_page=200, order_by="id", order="asc", id_above=id_above))
        res = d["results"]
        for r in res:
            g = r.get("geojson")
            if not g or not r.get("observed_on"):
                continue
            lon, lat = g["coordinates"]
            out.append([r["id"], r["observed_on"], round(lat, 4), round(lon, 4),
                        1 if (r.get("obscured") or r.get("geoprivacy") == "obscured"
                              or r.get("taxon_geoprivacy") == "obscured") else 0,
                        r["quality_grade"][0]])
        if len(res) < 200:
            break
        id_above = res[-1]["id"]
        time.sleep(1.0)
    return out


def build_inat(refresh):
    if os.path.exists(HIST_PATH) and not refresh:
        hist = json.load(open(HIST_PATH))
        print(f"[inat] using cached {HIST_PATH} (fetched {hist['fetched']})")
    else:
        hist = {"fetched": dt.date.today().isoformat(), "species": {}}
        for key, sp in M.SPECIES.items():
            obs = fetch_inat_all(sp["inat_taxa"])
            hist["species"][key] = obs
            print(f"[inat] {key}: {len(obs)} observations ({sum(o[4] for o in obs)} obscured)")
            time.sleep(1.0)
        os.makedirs(C.DATA_STATIC, exist_ok=True)
        json.dump(hist, open(HIST_PATH, "w"), separators=(",", ":"))
    kde = np.zeros((len(M.SPECIES_ORDER), C.NROWS, C.NCOLS), dtype=np.float32)
    for si, key in enumerate(M.SPECIES_ORDER):
        obs = np.array(hist["species"][key], dtype=object)
        if len(obs) == 0:
            continue
        r, c = C.cell_index(obs[:, 2].astype(float), obs[:, 3].astype(float))
        ok = r >= 0
        cnt = np.bincount(r[ok] * C.NCOLS + c[ok], minlength=N).reshape(C.NROWS, C.NCOLS).astype(np.float32)
        d = gaussian_blur(cnt, M.KDE_SIGMA_KM / C.CELL_KM_Y, M.KDE_SIGMA_KM / C.CELL_KM_X)
        nz = d[d > 0]
        scale = np.percentile(nz, 99) if nz.size else 1.0
        kde[si] = np.clip(d / max(scale, 1e-9), 0, 1)
        print(f"  {key:18s} obs in bbox {ok.sum():5d}  kde>0.25 cells {(kde[si] > 0.25).sum()}")
    return hist, kde


def gaussian_blur(a, sr, sc):
    def kern(s):
        r = int(math.ceil(3 * s))
        x = np.arange(-r, r + 1)
        k = np.exp(-0.5 * (x / s) ** 2)
        return k / k.sum()
    kr, kc = kern(sr), kern(sc)
    out = np.apply_along_axis(lambda v: np.convolve(v, kc, mode="same"), 1, a)
    out = np.apply_along_axis(lambda v: np.convolve(v, kr, mode="same"), 0, out)
    return out


# ---------------------------------------------------------------------------
# 5. PAD-US public land (federal USFS/BLM + state forests)
# ---------------------------------------------------------------------------
PADUS = "https://services.arcgis.com/v01gqwM5QqNysAAi/ArcGIS/rest/services/Manager_Name_PADUS/FeatureServer/0/query"
PADUS_WHERE = ("FeatClass = 'Fee' AND ((Mang_Name IN ('USFS','BLM')) OR "
               "(Mang_Type = 'STAT' AND Des_Tp IN ('FORE','SRMA')))")
PL_PATH = os.path.join(C.DATA_STATIC, "public_lands.geojson")


def build_padus(refresh):
    if os.path.exists(PL_PATH) and not refresh:
        print(f"[padus] using cached {PL_PATH}")
        return
    import shapely
    from shapely.geometry import shape, box, mapping
    feats, offset = [], 0
    while True:
        d = C.http_json(PADUS, dict(
            where=PADUS_WHERE, geometry=f"{C.LON_MIN},{C.LAT_MIN},{C.LON_MAX},{C.LAT_MAX}",
            geometryType="esriGeometryEnvelope", inSR=4326, spatialRel="esriSpatialRelIntersects",
            outFields="Unit_Nm,Mang_Name,Mang_Type,Des_Tp,Pub_Access", outSR=4326,
            maxAllowableOffset=0.001, resultOffset=offset, resultRecordCount=1000, f="geojson"))
        feats += d.get("features", [])
        print(f"[padus] fetched {len(feats)} features")
        if not d.get("properties", {}).get("exceededTransferLimit"):
            break
        offset += 1000
        time.sleep(0.5)
    bbox = box(*C.BBOX)
    out, managers = [], {}
    for f in feats:
        try:
            g = shape(f["geometry"])
        except Exception:
            continue
        if not g.is_valid:
            g = g.buffer(0)
        g = g.intersection(bbox).simplify(0.0015, preserve_topology=True)
        parts = [p for p in getattr(g, "geoms", [g]) if p.geom_type == "Polygon" and p.area >= 3e-5]
        if not parts:
            continue
        g = shapely.MultiPolygon(parts) if len(parts) > 1 else parts[0]
        pr = f["properties"]
        if (pr.get("Unit_Nm") or "").strip().lower() in ("not in ownership",):
            continue
        managers[pr["Mang_Name"]] = managers.get(pr["Mang_Name"], 0) + 1
        out.append({"type": "Feature", "properties": {
            "name": pr.get("Unit_Nm") or "", "mgr": pr.get("Mang_Name"), "type": pr.get("Des_Tp"),
            "access": pr.get("Pub_Access")}, "geometry": round_geom(mapping(g))})
    fc = {"type": "FeatureCollection", "features": out,
          "source": "USGS PAD-US 4.1 (Manager Name fee layer); simplified ~150 m"}
    json.dump(fc, open(PL_PATH, "w"), separators=(",", ":"))
    print(f"[padus] wrote {len(out)} polygons, {os.path.getsize(PL_PATH) / 1e6:.2f} MB, managers {managers}")


def round_geom(g):
    def rnd(c):
        if isinstance(c[0], (float, int)):
            return [round(c[0], 4), round(c[1], 4)]
        return [rnd(x) for x in c]
    return {"type": g["type"], "coordinates": rnd(g["coordinates"])}


# ---------------------------------------------------------------------------
# previews
# ---------------------------------------------------------------------------
def preview(name, arr, cmap="viridis", vmin=None, vmax=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    os.makedirs(PREVIEW, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 7 * (C.LAT_MAX - C.LAT_MIN) / (C.LON_MAX - C.LON_MIN) * 1.4))
    im = ax.imshow(arr, extent=[C.LON_MIN, C.LON_MAX, C.LAT_MIN, C.LAT_MAX], cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect=1.4)
    ax.set_title(name)
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.savefig(os.path.join(PREVIEW, f"{name}.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh-inat", action="store_true")
    ap.add_argument("--refresh-padus", action="store_true")
    ap.add_argument("--no-preview", action="store_true")
    args = ap.parse_args()
    os.makedirs(C.DATA_STATIC, exist_ok=True)

    elev = build_elevation()
    tree_frac, ocean_frac, water_frac = build_worldcover()
    ocean, coast_km = build_coast(ocean_frac, water_frac)
    eco_class, eco_l4, l4_codes, l4_names = build_ecoregions()
    hist, kde = build_inat(args.refresh_inat)
    build_padus(args.refresh_padus)

    C.save_npz(os.path.join(C.DATA_STATIC, "layers.npz"),
               elev=np.round(elev).astype(np.int16),
               tree_pct=np.round(tree_frac * 100).astype(np.uint8),
               eco_class=eco_class, eco_l4=eco_l4,
               coast_km=coast_km.astype(np.float16), ocean=ocean)
    C.save_npz(os.path.join(C.DATA_STATIC, "inat_kde.npz"), kde=np.round(kde * 255).astype(np.uint8))
    json.dump({
        "grid": C.grid_meta(),
        "eco_classes": M.ECO_CLASSES, "eco_labels": M.ECO_LABELS,
        "l4_codes": l4_codes, "l4_names": l4_names,
        "species": M.SPECIES_ORDER,
        "inat_fetched": hist["fetched"],
        "inat_counts": {k: len(v) for k, v in hist["species"].items()},
        "sources": {
            "elevation": "AWS Terrain Tiles (terrarium) z10, mean per cell",
            "forest": "ESA WorldCover 2021 v200 class 10 fraction (~75 m overview)",
            "ecoregions": "EPA Level IV ecoregions (OR, WA)",
            "inat": "iNaturalist research+needs_id observations, all years",
            "public_lands": "USGS PAD-US 4.1 fee layer",
        },
    }, open(os.path.join(C.DATA_STATIC, "grid.json"), "w"), indent=1)

    if not args.no_preview:
        preview("elevation", elev, "terrain", 0, 2500)
        preview("tree_fraction", tree_frac, "Greens", 0, 1)
        preview("eco_class", eco_class, "tab20", 0, 19)
        preview("coast_km", np.minimum(coast_km, 30), "Blues_r", 0, 30)
        for si, key in enumerate(M.SPECIES_ORDER):
            preview(f"kde_{key}", kde[si], "magma", 0, 1)
            h = M.habitat(M.SPECIES[key], elev, tree_frac, eco_class, coast_km)
            preview(f"habitat_{key}", h, "YlGn", 0, 1)
            print(f"  habitat {key:18s} mean {h.mean():.3f}  cells>0.5 {(h > 0.5).sum():6d}")
        print(f"previews in {PREVIEW}")
    print("done.")


if __name__ == "__main__":
    main()
