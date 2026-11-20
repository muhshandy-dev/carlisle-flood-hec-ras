"""Area 2D usulan: dataran rendah (tinggi di atas air terdekat <= 4 m) yang tersambung ke genangan 2005/2015."""
import numpy as np, rasterio, geopandas as gpd
from rasterio.features import rasterize, shapes
from rasterio.transform import from_origin
from scipy import ndimage
from shapely.geometry import shape, Polygon, MultiPolygon, box
from pathlib import Path
R = Path(__file__).resolve().parents[2]; OUT = R / "paket"; CRS = "EPSG:27700"
X0, X1, Y0, Y1, RS = 334500, 349500, 550500, 561500, 10
with rasterio.open(OUT / "02_medan" / "dtm_2m_carlisle.tif") as d:
    z = d.read(1, out_shape=((Y1 - Y0) // RS, (X1 - X0) // RS), resampling=rasterio.enums.Resampling.average).astype("float32")
tr = from_origin(X0, Y1, RS, RS); H, W = z.shape
air = gpd.read_file(OUT / "02_medan" / "badan_air_os.gpkg")
m_air = rasterize([(g, 1) for g in air.geometry], (H, W), transform=tr).astype(bool)
# hanya badan air memanjang (sungai), bukan kolam kecil
lab, n = ndimage.label(m_air); uk = ndimage.sum(m_air, lab, range(1, n + 1))
m_air = np.isin(lab, 1 + np.where(uk >= 200)[0])
_, (iy, ix) = ndimage.distance_transform_edt(~m_air, return_indices=True)
hand = z - z[iy, ix]
rfo = gpd.read_file(R / "probe" / "rfo_carlisle.gpkg")
ol = rfo[rfo["name"].isin(["Carlisle (including Crosby-on-Eden) 05_12_2015", "2005 Carlisle Event"])]
m_ol = rasterize([(g.buffer(300), 1) for g in ol.geometry], (H, W), transform=tr).astype(bool)
rendah = ((hand <= 4) & ~np.isnan(z)) | m_ol | m_air
lab, n = ndimage.label(rendah)
ids = np.unique(lab[m_ol & (lab > 0)])
m = np.isin(lab, ids)
m = ndimage.binary_closing(m, iterations=5); m = ndimage.binary_fill_holes(m)
m = ndimage.binary_opening(m, iterations=3)
poly = [shape(s) for s, v in shapes(m.astype("uint8"), mask=m, transform=tr)]
g = max(poly, key=lambda p: p.area).buffer(20).buffer(-20).simplify(10)
g = Polygon(g.exterior).intersection(box(X0, Y0, X1, Y1).buffer(-30))
if isinstance(g, MultiPolygon): g = max(g.geoms, key=lambda p: p.area)
gpd.GeoDataFrame({"nama": ["Area2D_Carlisle_usulan"], "luas_km2": [round(g.area / 1e6, 1)]}, geometry=[g], crs=CRS
                 ).to_file(OUT / "02_medan" / "area_2d_usulan.gpkg", driver="GPKG")
print("luas km2", round(g.area / 1e6, 1), "titik", len(g.exterior.coords))
