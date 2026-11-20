"""Langkah 1: DEM DAS (Copernicus 30 m -> EPSG:27700) dan sub-DAS dari batas resmi NRFA.

Sub-DAS HEC-HMS = daerah tangkapan tiap pos debit dikurangi daerah tangkapan pos di hulunya.
Batas resmi NRFA dipakai karena delineasi dari DEM 30 m (DSM, belum dikoreksi jembatan/tanggul) meleset 15-80%.
"""
import json, numpy as np, rasterio, geopandas as gpd, pandas as pd, requests
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.features import geometry_mask
from shapely.geometry import Point
from pathlib import Path

R = Path(__file__).resolve().parents[2]; RAW, K, OUT = R / "raw", R / "kerja", R / "paket"
(OUT / "01_das").mkdir(parents=True, exist_ok=True); K.mkdir(exist_ok=True)
CRS = "EPSG:27700"
X0, X1, Y0, Y1, RES = 295000, 395000, 490000, 585000, 30

# --- DEM 30 m untuk parameter sub-DAS (kemiringan) ---
src = [rasterio.open(p) for p in sorted(RAW.glob("Copernicus_*.tif"))]
mos, tr = merge(src)
W, H = (X1 - X0) // RES, (Y1 - Y0) // RES
dem = np.full((H, W), np.nan, "float32"); tr27 = from_origin(X0, Y1, RES, RES)
reproject(mos[0], dem, src_transform=tr, src_crs=src[0].crs, dst_transform=tr27, dst_crs=CRS,
          resampling=Resampling.bilinear, dst_nodata=np.nan)
prof = dict(driver="GTiff", height=H, width=W, count=1, dtype="float32", crs=CRS, transform=tr27,
            nodata=-9999, compress="deflate", tiled=True)
with rasterio.open(K / "dem30_27700.tif", "w", **prof) as d:
    d.write(np.where(np.isnan(dem), -9999, dem).astype("float32"), 1)

# --- pos debit (urut dari hulu ke hilir) dan pos tepat di hulunya ---
POS = {
    "76005": ("Eden_TempleSowerby", []),
    "76003": ("Eamont_Udford", []),
    "76017": ("Eden_GreatCorby", ["76005", "76003"]),
    "76008": ("Irthing_Greenholme", []),
    "76022": ("Caldew_Cummersdale", []),
    "76020": ("Petteril_NewbigginBridge", []),
    "76007": ("Eden_Sheepmount", ["76017", "76008", "76022", "76020"]),
}
WFS = ("https://hydrologygeoserver.ceh.ac.uk/geoserver/ows?service=WFS&version=2.0.0&typename=NRFA:NRFA_catchments"
       "&REQUEST=getfeature&CQL_FILTER=ID={}&outputformat=application/json&srsname=epsg:27700")
cat, info = {}, {}
for sid in POS:
    cat[sid] = gpd.read_file(WFS.format(sid)).to_crs(CRS).geometry.union_all()
    info[sid] = requests.get("https://nrfaapps.ceh.ac.uk/nrfa/ws/station-info",
                             params={"format": "json-object", "station": sid, "fields": "all"}, timeout=60).json()["data"][0]

# --- sub-DAS ---
baris = []
for sid, (nama, hulu) in POS.items():
    g = cat[sid]
    for h in hulu:
        g = g.difference(cat[h])
    nama_sub = "Lokal_Carlisle" if sid == "76007" else ("Eden_Tengah" if sid == "76017" else nama)
    baris.append({"sub_das": nama_sub, "outlet_pos": nama, "nrfa": sid, "geometry": g.buffer(0)})
sub = gpd.GeoDataFrame(baris, crs=CRS).explode(index_parts=False)
sub = sub[sub.area > 1e6].dissolve("sub_das", as_index=False, aggfunc="first")   # buang serpihan tepi batas
sub["luas_km2"] = (sub.area / 1e6).round(1)

gy, gx = np.gradient(np.where(np.isnan(dem), np.nanmean(dem), dem), RES)
slope_pct = (np.hypot(gx, gy) * 100).astype("float32")
sl, el = [], []
for geom in sub.geometry:
    m = ~geometry_mask([geom], (H, W), tr27)
    sl.append(float(np.nanmean(slope_pct[m]))); el.append(float(np.nanmean(dem[m])))
sub["kemiringan_pct"] = np.round(sl, 1); sub["elevasi_rata_m"] = np.round(el, 0)
# panjang aliran terpanjang: hukum Hack (nilai awal; dikoreksi saat kalibrasi)
sub["panjang_aliran_km"] = (1.4 * (sub.luas_km2 / 2.59) ** 0.568 * 1.609).round(1)
sub.to_file(OUT / "01_das" / "sub_das.gpkg", driver="GPKG")
gpd.GeoDataFrame({"nama": ["DAS_Eden_ke_Sheepmount"], "luas_km2": [round(cat["76007"].area / 1e6, 1)]},
                 geometry=[cat["76007"]], crs=CRS).to_file(OUT / "01_das" / "das_eden.gpkg", driver="GPKG")

HYD = "https://environment.data.gov.uk/hydrology/id/stations.json"
pos = []
for sid, (nama, _) in POS.items():
    it = requests.get(HYD, params={"nrfaStationID": sid}, timeout=60).json()["items"]
    pos.append({"pos": nama, "nrfa": sid, "ea_id": it[0]["notation"] if it else None,
                "luas_nrfa_km2": info[sid]["catchment-area"], "x": info[sid]["easting"], "y": info[sid]["northing"]})
pdf = pd.DataFrame(pos)
gpd.GeoDataFrame(pdf, geometry=[Point(a, b) for a, b in zip(pdf.x, pdf.y)], crs=CRS).to_file(
    OUT / "01_das" / "pos_debit.gpkg", driver="GPKG")
pdf.to_csv(K / "pos_debit.csv", index=False)
print(sub.drop(columns="geometry").to_string()); print(pdf.to_string())
print("Jumlah luas sub-DAS:", sub.luas_km2.sum(), "| NRFA Sheepmount:", info["76007"]["catchment-area"])
