import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import sys, rasterio, numpy as np, geopandas as gpd
from shapely.geometry import box
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from rasterio.features import rasterize
U, BG, TAG = sys.argv[1], sys.argv[2], sys.argv[3]
B=(337500,552500,348000,560500); R=4.0
def baca(p):
    with rasterio.open(p) as s:
        w=from_bounds(*B,s.transform); h=int((B[3]-B[1])/R); wd=int((B[2]-B[0])/R)
        a=s.read(1,window=w,out_shape=(h,wd),resampling=Resampling.average,masked=True,boundless=True)
        return a, rasterio.transform.from_bounds(*B,wd,h), s.crs
dtm,tr,crs=baca(PAKET+'02_medan/dtm_2m_carlisle_dengan_tanggul.tif')
wse,_,_=baca(U+'WSE (Max).Terrain.dtm_2m_carlisle_dengan_tanggul.tif')
dep,_,_=baca(U+'Depth (Max).Terrain.dtm_2m_carlisle_dengan_tanggul.tif')
bg=gpd.read_file(BG).to_crs(crs); bg=bg[bg.intersects(box(*B))].copy()
fp=rasterize([(g,1) for g in bg.geometry],dtm.shape,transform=tr).astype(bool)
d=np.ma.filled(dep,0).astype('float32'); w=np.ma.filled(wse,np.nan).astype('float32')
basah=(d>=0.05)&~np.isnan(w)&~fp
ter=np.ma.filled(dtm,np.nan).astype('float32'); ter[basah]=np.maximum(ter[basah],w[basah])
prof=dict(driver='GTiff',height=ter.shape[0],width=ter.shape[1],count=1,dtype='float32',crs=crs,transform=tr,nodata=-9999,compress='deflate')
with rasterio.open(f'terrain_banjir_{TAG}_4m.tif','w',**prof) as o: o.write(np.nan_to_num(ter,nan=-9999),1)
with rasterio.open(f'kedalaman_{TAG}_4m.tif','w',**prof) as o: o.write(np.where(d>=0.05,d,-9999).astype('float32'),1)
bg['tinggi_m']=bg['tinggi_m'].fillna(6).clip(lower=3)
keep=['tinggi_m','kedalaman_maks_m','kelas_bahaya','fungsi_penting','geometry']
bg[bg.kedalaman_maks_m>=0.15][keep].to_file(f'bangunan_tergenang_{TAG}.gpkg',driver='GPKG')
bg[bg.kedalaman_maks_m<0.15][keep].to_file(f'bangunan_aman_{TAG}.gpkg',driver='GPKG')
print(TAG, ter.shape, int((bg.kedalaman_maks_m>=0.15).sum()), len(bg), float(d.max()))
