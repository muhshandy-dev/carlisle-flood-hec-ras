import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import sys, numpy as np, pandas as pd, geopandas as gpd, rasterio, h5py, pyvista as pv
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from rasterio.features import rasterize
from scipy.spatial import cKDTree
from shapely.geometry import box
HDF, TAG, JUDUL, NF = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
P=PAKET+''
B=(338300,554000,345300,559200); R=4.0; VE=3.0
with rasterio.open(P+'02_medan/dtm_2m_carlisle_dengan_tanggul.tif') as s:
    h=int((B[3]-B[1])/R); w=int((B[2]-B[0])/R)
    dtm=s.read(1,window=from_bounds(*B,s.transform),out_shape=(h,w),resampling=Resampling.average,masked=True).filled(np.nan)
dtm=np.nan_to_num(dtm,nan=np.nanmin(dtm)); tr=rasterio.transform.from_bounds(*B,w,h)
f=h5py.File(HDF,'r'); bb='Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'
t=pd.to_datetime([x.decode() for x in f[bb+'Time Date Stamp'][:]],format='%d%b%Y %H:%M:%S')
ws=f[bb+'2D Flow Areas/Carlisle/Water Surface'][:]; n=ws.shape[1]
xy=f['Geometry/2D Flow Areas/Carlisle/Cells Center Coordinate'][:n]
area=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/area_2d_Carlisle.shp').geometry[0]
inside=rasterize([(area,1)],(h,w),transform=tr).astype(bool)
cols,rows=np.meshgrid(np.arange(w),np.arange(h)); X=B[0]+(cols+0.5)*R; Y=B[3]-(rows+0.5)*R
idx=np.zeros((h,w),int); _,k=cKDTree(xy).query(np.c_[X[inside],Y[inside]]); idx[inside]=k
bg=gpd.read_file(P+'04_lahan_tanah_bangunan/bangunan_dengan_tinggi.gpkg').to_crs(27700)
bg=bg[bg.intersects(box(*B))].reset_index(drop=True); bg['tinggi_m']=bg['tinggi_m'].fillna(6).clip(3,40)
lab=rasterize([(g,i+1) for i,g in enumerate(bg.geometry)],(h,w),transform=tr,dtype='int32'); bm=lab>0
ht=np.zeros((h,w)); ht[bm]=bg['tinggi_m'].values[lab[bm]-1]
zl=(dtm+ht)*VE
# warna medan: abu-abu-hijau lembut dengan bayangan
gy,gx=np.gradient(dtm,R); shade=np.clip(0.75+0.25*(gx-gy),0.4,1)
base=np.dstack([0.78*shade,0.80*shade,0.74*shade])
land=pv.StructuredGrid(X,Y,zl); water=pv.StructuredGrid(X,Y,(dtm-1)*VE)
ti=np.linspace(12*3600,(t[-1]-t[0]).total_seconds(),NF); t0=(t-t[0]).total_seconds().values
pl=pv.Plotter(off_screen=True,window_size=(1280,720)); pl.set_background('#dfe8f0',top='#9fb8d0')
def rgbflat(a): return (np.clip(a,0,1)*255).astype(np.uint8).reshape(-1,3,order='F')
land.point_data['rgb']=rgbflat(base); water.point_data['rgb']=rgbflat(np.zeros((h,w,3)))
pl.add_mesh(land,scalars='rgb',rgb=True,smooth_shading=False,lighting=True,ambient=0.35,diffuse=0.7)
pl.add_mesh(water,scalars='rgb',rgb=True,opacity=0.88,lighting=True,ambient=0.4,specular=0.3)
cx,cy=(B[0]+B[2])/2,(B[1]+B[3])/2
txt=pl.add_text(' ',position=(20,680),font_size=14,color='black',font='times')
pl.add_text(JUDUL,position='lower_left',font_size=11,color='black',font='times')
pl.add_text('Biru: genangan | Merah: bangunan tergenang | Vertikal x3',position=(760,12),font_size=9,color='black',font='times')
BL=['','Jan','Feb','Mar','Apr','Mei','Jun','Jul','Agu','Sep','Okt','Nov','Des']
import os; os.makedirs('fr_'+TAG,exist_ok=True)
for fi,tv in enumerate(ti):
    j=min(np.searchsorted(t0,tv,side='right')-1,len(t0)-2); a=(tv-t0[j])/(t0[j+1]-t0[j]); wv=(1-a)*ws[j]+a*ws[j+1]
    wse=np.where(inside,wv[idx],-999); d=wse-dtm; wet=inside&(d>=0.05)&~bm
    zw=np.where(wet,wse,dtm-1.0)*VE; water.points[:,2]=zw.ravel(order='F')
    c=np.zeros((h,w,3)); dd=np.clip(d/3,0,1)
    c[...,0]=0.55-0.5*dd; c[...,1]=0.75-0.5*dd; c[...,2]=0.95-0.4*dd
    water.point_data['rgb']=rgbflat(c)
    db=np.zeros(len(bg)+1); np.maximum.at(db,lab[bm],np.where(inside,d,0)[bm])
    col=base.copy(); st=db[lab]
    col[bm]=(0.93,0.93,0.91); col[bm&(st>=0.15)]=(0.88,0.25,0.2); col[bm&(st>=1.5)]=(0.6,0.07,0.07)
    land.point_data['rgb']=rgbflat(col)
    ang=np.deg2rad(215+25*fi/max(NF-1,1)); dist=4300
    pl.camera.position=(cx+dist*np.sin(ang),cy+dist*np.cos(ang),dist*0.42+40*VE); pl.camera.focal_point=(cx,cy+300,20*VE); pl.camera.up=(0,0,1)
    tm=t[0]+pd.Timedelta(seconds=float(tv))
    lbl=f'{tm.day:02d} {BL[tm.month]} {tm.year}  pukul {tm:%H:%M}' if TAG!='t100' else f'Jam ke-{tv/3600:04.1f} sejak awal hujan rencana'
    txt.SetInput(lbl+f'    |    bangunan tergenang: {int((db[1:]>=0.15).sum()):,}'.replace(',','.'))
    txt.Modified(); pl.render(); pl.screenshot(f'fr_{TAG}/{fi:04d}.png')
print('selesai',NF)
