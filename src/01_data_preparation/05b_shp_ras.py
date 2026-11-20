"""Langkah 5b: file siap-impor RAS Mapper (shapefile + .prj) dan DTM dengan tanggul 'dibakar'."""
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from shapely.geometry import LineString, Point, MultiLineString
from shapely.ops import nearest_points, linemerge
from pyproj import CRS as PCRS
from pathlib import Path

R = Path(__file__).resolve().parents[2]; P = R / "paket"; OR = P / "07_input_hec_ras"; S = OR / "shp_untuk_RAS_Mapper"
S.mkdir(exist_ok=True); CRS = "EPSG:27700"
prj = PCRS.from_epsg(27700).to_wkt(version="WKT1_ESRI")
(OR / "British_National_Grid_EPSG27700.prj").write_text(prj)

area = gpd.read_file(P / "02_medan" / "area_2d_usulan.gpkg").geometry.iloc[0]
pos = gpd.read_file(P / "01_das" / "pos_debit.gpkg").set_index("pos")
# potong koridor Caldew dan Irthing tepat di pos debit, supaya batas hulu = lokasi pos
from shapely.geometry import box
cal, irt = pos.loc["Caldew_Cummersdale"].geometry, pos.loc["Irthing_Greenholme"].geometry
CUT = {"BC_Caldew": (LineString([(335000, cal.y - 100), (341000, cal.y - 100)]), box(335000, 550000, 341000, cal.y - 100)),
       "BC_Irthing": (LineString([(irt.x + 100, 557000), (irt.x + 100, 561500)]), box(irt.x + 100, 557000, 350000, 561500))}
for nm, (ln, bx) in CUT.items():
    area = area.difference(bx)
from shapely.geometry import Polygon, MultiPolygon
if isinstance(area, MultiPolygon): area = max(area.geoms, key=lambda g: g.area)
area = Polygon(area.exterior)
gpd.GeoDataFrame({"nama": ["Area2D_Carlisle_usulan"], "luas_km2": [round(area.area / 1e6, 1)]}, geometry=[area], crs=CRS
                 ).to_file(P / "02_medan" / "area_2d_usulan.gpkg", driver="GPKG")
ext = area.exterior
air = gpd.read_file(P / "02_medan" / "badan_air_os.gpkg")
air_u = air[air.area > 2000].union_all()
tma = gpd.read_file(P / "03_hidromet" / "pos_tma_dan_debit.gpkg")

def potong_tepi(titik, panjang):
    """Segmen tepi area 2D sepanjang `panjang` m yang berpusat di titik terdekat pada tepi."""
    d0 = ext.project(titik); L = ext.length
    t = np.linspace(d0 - panjang / 2, d0 + panjang / 2, 60) % L
    return LineString([ext.interpolate(v) for v in t])

def silang_sungai(titik_acuan, rad=1500):
    """Titik tepi area 2D tempat sungai keluar/masuk, paling dekat dengan titik acuan."""
    seg = ext.intersection(air_u.buffer(5))
    kand = [g.centroid for g in getattr(seg, "geoms", [seg]) if not g.is_empty and g.centroid.distance(titik_acuan) < rad]
    return min(kand, key=lambda p: p.distance(titik_acuan)) if kand else nearest_points(ext, titik_acuan)[0]

bc = []
for nama, gauge, lebar in [("BC_Eden_hulu", "Eden_GreatCorby", 500), ("BC_Petteril", "Petteril_NewbigginBridge", 300)]:
    c = silang_sungai(pos.loc[gauge].geometry)
    bc.append({"nama": nama, "jenis": "Flow Hydrograph", "geometry": potong_tepi(c, lebar)})
for nm, (ln, bx) in CUT.items():
    g = ln.intersection(area.buffer(1))
    if g.geom_type == "MultiLineString": g = max(g.geoms, key=lambda q: q.length)
    bc.append({"nama": nm, "jenis": "Flow Hydrograph", "geometry": g})
# hilir: Eden keluar di tepi barat laut domain
c = silang_sungai(Point(336000, 561400), rad=3000)
bc.append({"nama": "BC_Eden_hilir", "jenis": "Normal Depth", "geometry": potong_tepi(c, 900)})
gbc = gpd.GeoDataFrame(bc, crs=CRS); gbc.to_file(S / "garis_batas_BC.shp")
print(gbc.assign(panjang=gbc.length.round(0), x=gbc.centroid.x.round(0), y=gbc.centroid.y.round(0)).drop(columns="geometry"))

gpd.GeoDataFrame({"nama": ["Carlisle"]}, geometry=[area], crs=CRS).to_file(S / "area_2d_Carlisle.shp")

# reference points: pos TMA di dalam area 2D yang punya datum
tma["dalam_area_2d"] = tma.within(area.buffer(100))
tma.to_file(P / "03_hidromet" / "pos_tma_dan_debit.gpkg", driver="GPKG")
rp = tma[tma.dalam_area_2d & tma.datum_tersedia & tma.within(area)][["nama", "sungai", "datum_mAOD", "geometry"]]
rp = rp.copy(); rp["nama"] = rp.nama.str.slice(0, 32)
rp.to_file(S / "reference_points_TMA.shp"); print("reference points:", list(rp.nama))

# reference line di Sheepmount: garis utara-selatan melintang dataran banjir Eden, dipotong area 2D
sm = pos.loc["Eden_Sheepmount"].geometry
garis = LineString([(sm.x, sm.y - 2000), (sm.x, sm.y + 2000)]).intersection(area)
if garis.geom_type == "MultiLineString":
    garis = min(garis.geoms, key=lambda g: g.distance(sm))
gpd.GeoDataFrame({"nama": ["Eden_Sheepmount"]}, geometry=[garis], crs=CRS).to_file(S / "reference_line_debit.shp")
print("reference line panjang m:", round(garis.length))

# breakline tanggul + DTM dengan crest tanggul dibakar
tg = gpd.read_file(OR / "tanggul_tembok_banjir_EA.gpkg")
tg = tg[tg.geometry.geom_type.isin(["LineString", "MultiLineString"])].copy()
tg = tg.explode(index_parts=False); tg = tg[tg.length > 5]
tg[["asset_id", "asset_sub_type", "crest_mAOD", "geometry"]].rename(columns={"asset_sub_type": "jenis", "crest_mAOD": "crest"}
    ).to_file(S / "breakline_tanggul.shp")
with rasterio.open(P / "02_medan" / "dtm_2m_carlisle.tif") as d:
    z = d.read(1); prof = d.profile; tr = d.transform
ok = tg[np.isfinite(tg.crest_mAOD)]
crest = rasterize([(g, float(c)) for g, c in zip(ok.geometry, ok.crest_mAOD)], z.shape, transform=tr, fill=-9999,
                  all_touched=True, dtype="float32", merge_alg=rasterio.enums.MergeAlg.replace)
naik = (crest > z) & (z > -1000)
z2 = np.where(naik, crest, z).astype("float32")
with rasterio.open(P / "02_medan" / "dtm_2m_carlisle_dengan_tanggul.tif", "w", **prof) as d:
    d.write(z2, 1)
print("sel dinaikkan:", int(naik.sum()), "| kenaikan median m:", round(float(np.median((crest - z)[naik])), 2))
