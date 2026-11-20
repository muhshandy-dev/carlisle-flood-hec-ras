"""Langkah 4: tutupan lahan -> Manning n; tekstur tanah -> HSG -> CN; bangunan + tinggi; jalan; outline genangan."""
import numpy as np, pandas as pd, geopandas as gpd, rasterio, requests, time
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.features import rasterize, geometry_mask
from pathlib import Path

R = Path(__file__).resolve().parents[2]; RAW, K, OUT = R / "raw", R / "kerja", R / "paket"
O4 = OUT / "04_lahan_tanah_bangunan"; O4.mkdir(parents=True, exist_ok=True); CRS = "EPSG:27700"

def unduh(url, f, timeout=900):
    if f.exists() and f.stat().st_size > 10000: return f
    for i in range(4):
        try:
            r = requests.get(url, timeout=timeout, stream=True)
            if r.ok:
                with open(f, "wb") as o:
                    for ch in r.iter_content(1 << 20): o.write(ch)
                return f
        except Exception: pass
        time.sleep(10)
    raise RuntimeError(url)

# ---------- grid 10 m DAS dan grid 2 m domain 2D ----------
das = gpd.read_file(OUT / "01_das" / "das_eden.gpkg").geometry.iloc[0]
bx = [int(v // 1000 * 1000) for v in das.bounds]; X0, Y0, X1, Y1 = bx[0], bx[1], bx[2] + 1000, bx[3] + 1000
RS = 10; W, H = (X1 - X0) // RS, (Y1 - Y0) // RS; TR = from_origin(X0, Y1, RS, RS)
D0, D1, E0, E1 = 334500, 349500, 550500, 561500     # domain 2D

def ke_grid(src_path, tr, shape, resamp=Resampling.nearest, dtype="uint8"):
    out = np.zeros(shape, dtype)
    with rasterio.open(src_path) as s:
        reproject(rasterio.band(s, 1), out, dst_transform=tr, dst_crs=CRS, resampling=resamp)
    return out

# ---------- ESA WorldCover 2021 ----------
wc_tiles = []
for t in ["N54W003", "N54W006"]:
    f = RAW / f"wc_{t}.tif"
    unduh(f"https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{t}_Map.tif", f)
    wc_tiles.append(f)
from rasterio.merge import merge
src = [rasterio.open(f) for f in wc_tiles]
mos, mtr = merge(src, bounds=(-3.4, 54.3, -2.0, 55.25))
wc = np.zeros((H, W), "uint8")
reproject(mos[0], wc, src_transform=mtr, src_crs=src[0].crs, dst_transform=TR, dst_crs=CRS, resampling=Resampling.mode)
KELAS = {10: "Pohon", 20: "Semak", 30: "Rumput", 40: "Pertanian", 50: "Terbangun", 60: "Lahan terbuka",
         80: "Air", 90: "Lahan basah", 95: "Mangrove", 100: "Lumut/lichen"}
MANNING = {10: 0.12, 20: 0.07, 30: 0.035, 40: 0.04, 50: 0.10, 60: 0.03, 80: 0.035, 90: 0.07, 95: 0.12, 100: 0.04}
# CN AMC II per HSG (A,B,C,D) — tabel TR-55 yang disederhanakan
CN = {10: (30, 55, 70, 77), 20: (35, 56, 70, 77), 30: (39, 61, 74, 80), 40: (64, 75, 82, 85),
      50: (77, 85, 90, 92), 60: (77, 86, 91, 94), 80: (100, 100, 100, 100), 90: (85, 85, 85, 85),
      95: (85, 85, 85, 85), 100: (49, 69, 79, 84)}
pd.DataFrame([{"kode_worldcover": k, "kelas": v, "manning_n": MANNING[k], "CN_A": CN[k][0], "CN_B": CN[k][1],
               "CN_C": CN[k][2], "CN_D": CN[k][3]} for k, v in KELAS.items()]).to_csv(O4 / "tabel_manning_cn.csv", index=False)

# raster tutupan lahan 2 m untuk domain 2D (HEC-RAS: Land Cover layer)
wd = np.zeros(((E1 - E0) // 2, (D1 - D0) // 2), "uint8")
reproject(mos[0], wd, src_transform=mtr, src_crs=src[0].crs, dst_transform=from_origin(D0, E1, 2, 2), dst_crs=CRS,
          resampling=Resampling.nearest)
# perbaiki dengan OS OpenMap: bangunan = 50 (Terbangun), badan air = 80
oml = next((RAW / "oml").rglob("data"))
bang_all = gpd.read_file(oml / "NY_Building.shp", bbox=(D0, E0, D1, E1)).to_crs(CRS)
air = gpd.read_file(OUT / "02_medan" / "badan_air_os.gpkg")
trd = from_origin(D0, E1, 2, 2)
wd[rasterize([(g, 1) for g in air.geometry], wd.shape, transform=trd).astype(bool)] = 80
wd[rasterize([(g, 1) for g in bang_all.geometry], wd.shape, transform=trd).astype(bool)] = 50
prof = dict(driver="GTiff", count=1, crs=CRS, compress="deflate", tiled=True, blockxsize=512, blockysize=512)
with rasterio.open(O4 / "tutupan_lahan_2m_domain2d.tif", "w", height=wd.shape[0], width=wd.shape[1],
                   transform=trd, dtype="uint8", nodata=0, **prof) as d:
    d.write(wd, 1)
    d.write_colormap(1, {10: (0, 100, 0, 255), 20: (255, 187, 34, 255), 30: (255, 255, 76, 255), 40: (240, 150, 255, 255),
                         50: (250, 0, 0, 255), 60: (180, 180, 180, 255), 80: (0, 100, 200, 255), 90: (0, 150, 160, 255),
                         95: (0, 207, 117, 255), 100: (250, 230, 160, 255)})

# ---------- tekstur tanah SoilGrids (5-15 cm) -> HSG ----------
def soilgrids(var):
    f = RAW / f"sg_{var}.tif"
    url = (f"https://maps.isric.org/mapserv?map=/map/{var}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
           f"&COVERAGEID={var}_5-15cm_mean&FORMAT=image/tiff&SUBSET=long(-3.4,-2.0)&SUBSET=lat(54.3,55.25)"
           "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326")
    unduh(url, f, timeout=600)
    return ke_grid(f, TR, (H, W), Resampling.bilinear, "float32") / 10.0     # g/kg -> %
clay, sand, silt = soilgrids("clay"), soilgrids("sand"), soilgrids("silt")
def hsg(c, s, si):
    # kelas tekstur USDA (ringkas) -> HSG: A pasir/pasir berlempung; B lempung berpasir/lempung/lempung berdebu;
    # C lempung liat berpasir; D liat & lempung berliat
    g = np.full(c.shape, 2, "uint8")                                   # default B
    g[(s >= 85) | ((s >= 70) & (c < 15))] = 1                           # A
    g[(c >= 20) & (c < 35) & (s >= 45)] = 3                             # C: sandy clay loam
    g[(c >= 27) & ~((c < 35) & (s >= 45))] = 4                          # D: clay loam, silty clay loam, clay
    g[(c + s + si) < 1] = 0
    return g
HSG = hsg(clay, sand, silt)
cn = np.zeros((H, W), "float32")
for k, v in CN.items():
    for i in range(4):
        cn[(wc == k) & (HSG == i + 1)] = v[i]
cn[(cn == 0) & (wc > 0)] = 75
with rasterio.open(O4 / "CN_10m_DAS.tif", "w", height=H, width=W, transform=TR, dtype="float32", nodata=0, **prof) as d:
    d.write(cn, 1)
with rasterio.open(O4 / "HSG_10m_DAS.tif", "w", height=H, width=W, transform=TR, dtype="uint8", nodata=0, **prof) as d:
    d.write(HSG, 1)

# parameter per sub-DAS: CN rata-rata, % kedap air, Manning lahan, lag SCS
sub = gpd.read_file(OUT / "01_das" / "sub_das.gpkg")
rows = []
for r in sub.itertuples():
    m = ~geometry_mask([r.geometry], (H, W), TR)
    cn_ = float(np.nanmean(cn[m])); kedap = float(100 * np.mean(wc[m] == 50))
    S = 1000 / cn_ - 10                                        # inci
    L_ft = r.panjang_aliran_km * 3280.84; Y = max(r.kemiringan_pct, 0.5)
    lag_jam = L_ft ** 0.8 * (S + 1) ** 0.7 / (1900 * Y ** 0.5)  # rumus lag SCS
    tc = lag_jam / 0.6
    porsi = {KELAS[k]: round(100 * float(np.mean(wc[m] == k)), 1) for k in KELAS if np.mean(wc[m] == k) > 0.005}
    rows.append({"sub_das": r.sub_das, "luas_km2": r.luas_km2, "CN_rata": round(cn_, 1), "kedap_air_pct": round(kedap, 1),
                 "Ia_awal_mm": round(0.2 * S * 25.4, 1), "kemiringan_pct": r.kemiringan_pct,
                 "panjang_aliran_km": r.panjang_aliran_km, "Tc_jam": round(tc, 2), "lag_SCS_jam": round(lag_jam, 2),
                 "lag_SCS_menit": round(lag_jam * 60), "tutupan_lahan_pct": porsi})
par = pd.DataFrame(rows); par.to_csv(O4 / "parameter_awal_hec_hms_per_subdas.csv", index=False)
print(par.drop(columns="tutupan_lahan_pct").to_string())

# ---------- bangunan + tinggi dari nDSM ----------
dom_bang = bang_all.copy()
with rasterio.open(K / "ndsm_1m.tif") as d:
    nd = d.read(1); ndtr = d.transform
dom_bang = dom_bang[dom_bang.area >= 10].reset_index(drop=True)
lab = rasterize([(g, i + 1) for i, g in enumerate(dom_bang.geometry)], nd.shape, transform=ndtr, dtype="int32")
v = nd.ravel(); l = lab.ravel(); ok = (l > 0) & (v > -1)
df = pd.DataFrame({"l": l[ok], "v": v[ok]})
med = df.groupby("l").v.median(); p90 = df.groupby("l").v.quantile(0.9)
dom_bang["tinggi_m"] = dom_bang.index.map(lambda i: med.get(i + 1, np.nan)).astype(float).round(1)
dom_bang["tinggi_p90_m"] = dom_bang.index.map(lambda i: p90.get(i + 1, np.nan)).astype(float).round(1)
dom_bang["tinggi_m"] = dom_bang.tinggi_m.where(dom_bang.tinggi_m >= 2.5, np.nan)   # nilai terlalu rendah = salah ukur
dom_bang["tinggi_m"] = dom_bang.tinggi_m.fillna(dom_bang.tinggi_m.median()).round(1)
dom_bang["jumlah_lantai_perkiraan"] = np.maximum(1, np.round(dom_bang.tinggi_m / 3.0)).astype(int)
dom_bang["luas_m2"] = dom_bang.area.round(1)
# elevasi tanah di bangunan (untuk kedalaman di bangunan)
dtm1 = np.load(K / "dtm1.npy")
zz = pd.DataFrame({"l": l[ok], "z": dtm1.ravel()[ok]}).groupby("l").z.min()
dom_bang["elev_tanah_min_mAOD"] = dom_bang.index.map(lambda i: zz.get(i + 1, np.nan)).astype(float).round(2)
penting = gpd.read_file(oml / "NY_ImportantBuilding.shp", bbox=(D0, E0, D1, E1)).to_crs(CRS)
dom_bang = gpd.sjoin(dom_bang, penting[["CLASSIFICA", "geometry"]] if "CLASSIFICA" in penting else penting[["geometry"]],
                     how="left", predicate="intersects").drop(columns="index_right", errors="ignore")
dom_bang = dom_bang[~dom_bang.index.duplicated()]
dom_bang.rename(columns={"CLASSIFICA": "fungsi_penting"}).to_file(O4 / "bangunan_dengan_tinggi.gpkg", driver="GPKG")
print("bangunan:", len(dom_bang), "| median tinggi:", dom_bang.tinggi_m.median())
jalan = gpd.read_file(oml / "NY_Road.shp", bbox=(D0, E0, D1, E1)).to_crs(CRS)
jalan.to_file(O4 / "jalan_os.gpkg", driver="GPKG")

# ---------- outline genangan resmi per kejadian (validasi) ----------
rfo = gpd.read_file(R / "probe" / "rfo_carlisle.gpkg")
O5 = OUT / "05_validasi"; O5.mkdir(exist_ok=True)
peta = {"desmond_2015": ["Carlisle (including Crosby-on-Eden) 05_12_2015"], "jan_2005": ["2005 Carlisle Event"],
        "ciara_2020": ["Carlisle Storm Ciara 09-02-2020"]}
for k, names in peta.items():
    g = rfo[rfo["name"].isin(names)]
    g.to_file(O5 / f"genangan_resmi_{k}.gpkg", driver="GPKG")
    print(k, len(g), round(g.area.sum() / 1e6, 2), "km2")
