"""Langkah 2: medan domain 2D Carlisle dari LiDAR EA (WCS), area 2D usulan, garis sungai, nDSM untuk tinggi bangunan."""
import numpy as np, rasterio, geopandas as gpd, requests, time, io
from rasterio.transform import from_origin
from concurrent.futures import ThreadPoolExecutor
from shapely.geometry import box, Point, LineString
from shapely.ops import nearest_points
from pathlib import Path

R = Path(__file__).resolve().parents[2]; RAW, K, OUT = R / "raw", R / "kerja", R / "paket"
(OUT / "02_medan").mkdir(parents=True, exist_ok=True); (RAW / "lidar").mkdir(exist_ok=True)
CRS = "EPSG:27700"
X0, X1, Y0, Y1 = 334500, 349500, 550500, 561500          # kotak domain 2D (BNG, m)
SVC = {
    "dtm": ("lidar-composite-digital-terrain-model-dtm-1m", "13787b9a-26a4-4775-8523-806d13af58fc__Lidar_Composite_Elevation_DTM_1m"),
    "dsm": ("lidar-composite-digital-surface-model-first-return-dsm-1m", "df4e3ec3-315e-48aa-aaaf-b5ae74d7b2bb__Lidar_Composite_Elevation_FZ_DSM_1m"),
}

def unduh(jenis, x, y):
    f = RAW / "lidar" / f"{jenis}_{x}_{y}.tif"
    if f.exists() and f.stat().st_size > 1000:
        return f
    svc, cid = SVC[jenis]
    url = (f"https://environment.data.gov.uk/spatialdata/{svc}/wcs?service=WCS&version=2.0.1&request=GetCoverage"
           f"&CoverageId={cid}&format=image/tiff&subset=E({x},{x+1000})&subset=N({y},{y+1000})")
    for i in range(5):
        try:
            r = requests.get(url, timeout=180)
            if r.ok and r.headers.get("content-type", "").startswith("image/tiff"):
                f.write_bytes(r.content); return f
        except Exception:
            pass
        time.sleep(5 * (i + 1))
    return None

tiles = [(x, y) for x in range(X0, X1, 1000) for y in range(Y0, Y1, 1000)]
with ThreadPoolExecutor(6) as ex:
    hasil = list(ex.map(lambda a: (a[0], a[1:], unduh(*a)), [(j, x, y) for j in SVC for x, y in tiles]))
gagal = [h for h in hasil if h[2] is None]
print("ubin:", len(hasil), "gagal:", len(gagal))

W, H = X1 - X0, Y1 - Y0
mos = {j: np.full((H, W), np.nan, "float32") for j in SVC}
for j, (x, y), f in hasil:
    if f is None: continue
    with rasterio.open(f) as d:
        a = d.read(1).astype("float32"); a[a < -1000] = np.nan
        c0, r0 = int(d.transform.c - X0), int(Y1 - d.transform.f)
        h, w = a.shape
        mos[j][r0:r0 + h, c0:c0 + w] = a
dtm1, dsm1 = mos["dtm"], mos["dsm"]
print("celah DTM %:", round(100 * np.isnan(dtm1).mean(), 3))

# DTM 2 m (rata-rata 2x2) untuk HEC-RAS
d2 = np.nanmean(dtm1.reshape(H // 2, 2, W // 2, 2), axis=(1, 3)).astype("float32")
prof = dict(driver="GTiff", count=1, dtype="float32", crs=CRS, nodata=-9999, compress="deflate", predictor=3,
            tiled=True, blockxsize=512, blockysize=512, BIGTIFF="IF_SAFER")
with rasterio.open(OUT / "02_medan" / "dtm_2m_carlisle.tif", "w", height=H // 2, width=W // 2,
                   transform=from_origin(X0, Y1, 2, 2), **prof) as d:
    d.write(np.where(np.isnan(d2), -9999, d2), 1)
# nDSM 1 m (tinggi objek di atas tanah) untuk tinggi bangunan
ndsm = np.clip(dsm1 - dtm1, 0, 120).astype("float32")
with rasterio.open(K / "ndsm_1m.tif", "w", height=H, width=W, transform=from_origin(X0, Y1, 1, 1), **prof) as d:
    d.write(np.where(np.isnan(ndsm), -9999, ndsm), 1)
np.save(K / "dtm1.npy", dtm1)

# --- garis sungai (OS OpenMap Local) dan area 2D usulan ---
dom = box(X0, Y0, X1, Y1)
oml = next((RAW / "oml").rglob("data"))
sw = gpd.read_file(oml / "NY_SurfaceWater_Line.shp", bbox=(X0, Y0, X1, Y1)).to_crs(CRS)
sw = sw.clip(dom); sw.to_file(OUT / "02_medan" / "garis_sungai_os.gpkg", driver="GPKG")
swa = gpd.read_file(oml / "NY_SurfaceWater_Area.shp", bbox=(X0, Y0, X1, Y1)).to_crs(CRS).clip(dom)
swa.to_file(OUT / "02_medan" / "badan_air_os.gpkg", driver="GPKG")

rfo = gpd.read_file(R / "probe" / "rfo_carlisle.gpkg")
ol = rfo[rfo.name.isin(["Carlisle (including Crosby-on-Eden) 05_12_2015", "2005 Carlisle Event"])].union_all()
pos = gpd.read_file(OUT / "01_das" / "pos_debit.gpkg")
koridor = [ol.buffer(400)]
for _, p in pos[pos.pos != "Eden_Sheepmount"].iterrows():
    if dom.contains(p.geometry):
        a, b = nearest_points(p.geometry, ol)
        koridor.append(LineString([a, b]).buffer(300)); koridor.append(p.geometry.buffer(300))
area2d = gpd.GeoSeries(koridor, crs=CRS).union_all().intersection(dom.buffer(-50))
from shapely.geometry import Polygon, MultiPolygon
def tanpa_lubang(g):
    if isinstance(g, Polygon): return Polygon(g.exterior)
    return MultiPolygon([Polygon(p.exterior) for p in g.geoms])
area2d = tanpa_lubang(area2d.buffer(50).buffer(-50))
if isinstance(area2d, MultiPolygon):
    area2d = max(area2d.geoms, key=lambda g: g.area)
gpd.GeoDataFrame({"nama": ["Area2D_Carlisle_usulan"], "luas_km2": [round(area2d.area / 1e6, 1)]},
                 geometry=[area2d], crs=CRS).to_file(OUT / "02_medan" / "area_2d_usulan.gpkg", driver="GPKG")
print("area 2D usulan km2:", round(area2d.area / 1e6, 1))
