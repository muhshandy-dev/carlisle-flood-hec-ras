import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import sys, numpy as np, pandas as pd, geopandas as gpd, rasterio, h5py
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from rasterio.features import rasterize
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif'],'font.size':10})
HDF, TAG, JUDUL, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
P=PAKET+''
B=(337500,552500,348000,560500); R=4.0
with rasterio.open(P+'02_medan/dtm_2m_carlisle_dengan_tanggul.tif') as s:
    h=int((B[3]-B[1])/R); w=int((B[2]-B[0])/R)
    dtm=s.read(1,window=from_bounds(*B,s.transform),out_shape=(h,w),resampling=Resampling.average,masked=True).filled(np.nan)
tr=rasterio.transform.from_bounds(*B,w,h)
f=h5py.File(HDF,'r'); bb='Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'
t=pd.to_datetime([x.decode() for x in f[bb+'Time Date Stamp'][:]],format='%d%b%Y %H:%M:%S')
ws=f[bb+'2D Flow Areas/Carlisle/Water Surface'][:]; n=ws.shape[1]
g='Geometry/2D Flow Areas/Carlisle/'; xy=f[g+'Cells Center Coordinate'][:n]
q=np.abs(f[bb+'2D Flow Areas/Carlisle/Boundary Conditions/BC_Eden_hilir'][:][:,1])
area=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/area_2d_Carlisle.shp').geometry[0]
inside=rasterize([(area,1)],(h,w),transform=tr).astype(bool)
cols,rows=np.meshgrid(np.arange(w),np.arange(h)); X=B[0]+(cols+0.5)*R; Y=B[3]-(rows+0.5)*R
idx=np.full((h,w),-1); _,k=cKDTree(xy).query(np.c_[X[inside],Y[inside]]); idx[inside]=k
gy,gx=np.gradient(np.nan_to_num(dtm,nan=np.nanmean(dtm)),R); hs=np.clip(0.62+0.35*(gx-gy),0,1)
bg=gpd.read_file(P+'04_lahan_tanah_bangunan/bangunan_dengan_tinggi.gpkg').to_crs(27700)
from shapely.geometry import box; bg=bg[bg.intersects(box(*B))].reset_index(drop=True)
lab=rasterize([(gm,i+1) for i,gm in enumerate(bg.geometry)],(h,w),transform=tr,dtype='int32')
bmask=lab>0; nb=len(bg)
cm=LinearSegmentedColormap.from_list('air',['#c6e2f5','#4a90c8','#08306b'])
# interpolasi 15 menit
tt=pd.date_range(t[0],t[-1],freq='15min'); ti=(tt-t[0]).total_seconds().values; t0=(t-t[0]).total_seconds().values
dmax_b=np.zeros(nb+1)
fig=plt.figure(figsize=(12.8,7.2),dpi=100); ax=fig.add_axes([0.01,0.02,0.70,0.90]); axq=fig.add_axes([0.75,0.58,0.23,0.28])
for fi,tv in enumerate(ti):
    j=np.searchsorted(t0,tv,side='right')-1; j=min(j,len(t0)-2); a=(tv-t0[j])/(t0[j+1]-t0[j])
    wv=(1-a)*ws[j]+a*ws[j+1]
    d=np.zeros((h,w),'float32'); d[inside]=wv[idx[inside]]-dtm[inside]; d[~np.isfinite(d)]=0
    db=np.where(bmask,0,d)
    # status bangunan: kedalaman di sekeliling tapak (piksel tapak pakai sel terdekat)
    dbld=np.zeros(nb+1); np.maximum.at(dbld,lab[bmask],d[bmask]); dmax_b=np.maximum(dmax_b,dbld)
    rgb=np.dstack([hs]*3)*0.92+0.05
    wet=(db>=0.05)&inside
    c=cm(np.clip(db/3,0,1))[...,:3]; al=np.clip(0.55+db/4,0,0.92)[...,None]
    rgb[wet]=(rgb[wet]*(1-al[wet])+c[wet]*al[wet])
    st=dbld[lab]; rgb[bmask]=(0.86,0.86,0.86); rgb[bmask&(st>=0.15)]=(0.85,0.2,0.18); rgb[bmask&(st>=1.5)]=(0.55,0.05,0.05)
    ax.clear(); ax.imshow(rgb,extent=[B[0],B[2],B[1],B[3]]); ax.set_xticks([]); ax.set_yticks([])
    BL=['','Jan','Feb','Mar','Apr','Mei','Jun','Jul','Agu','Sep','Okt','Nov','Des']
    lbl=(f'{tt[fi].day:02d} {BL[tt[fi].month]} {tt[fi].year}  pukul {tt[fi]:%H:%M}') if TAG!='t100' else f'Jam ke-{(tt[fi]-tt[0]).total_seconds()/3600:04.1f} sejak awal hujan rencana'
    ax.set_title(f'{JUDUL}   |   {lbl}',fontsize=14,loc='left')
    nb_now=int(((dbld[1:]>=0.15)).sum())
    ax.text(0.01,0.02,f'Bangunan tergenang saat ini: {nb_now:,}'.replace(',','.')+f'\nBangunan pernah tergenang: {int((dmax_b[1:]>=0.15).sum()):,}'.replace(',','.'),
        transform=ax.transAxes,fontsize=11,va='bottom',bbox=dict(fc='white',ec='0.6',alpha=0.85))
    axq.clear(); axq.plot(t,q,color='#1f4e79',lw=1.4); axq.axvline(tt[fi],color='#c0392b',lw=1.2)
    axq.set_title('Debit keluar hilir (m³/s)',fontsize=10); axq.tick_params(labelsize=8); axq.grid(alpha=.3)
    axq.xaxis.set_major_locator(matplotlib.dates.DayLocator() if TAG!='t100' else matplotlib.dates.HourLocator(byhour=[0,12])); axq.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m') if TAG!='t100' else matplotlib.dates.DateFormatter('%H:%M'))
    if fi==0:
        lg=fig.add_axes([0.75,0.18,0.23,0.30]); lg.axis('off')
        import matplotlib.patches as mp
        lg.legend(handles=[mp.Patch(color=cm(0.05),label='Genangan dangkal'),mp.Patch(color=cm(0.5),label='Genangan sedang'),mp.Patch(color=cm(1.0),label='Genangan > 3 m'),
            mp.Patch(color=(0.86,0.86,0.86),label='Bangunan kering'),mp.Patch(color=(0.85,0.2,0.18),label='Bangunan tergenang ≥ 0,15 m'),mp.Patch(color=(0.55,0.05,0.05),label='Bangunan tergenang > 1,5 m')],
            loc='upper left',fontsize=10,frameon=False)
        fig.text(0.75,0.08,'Model: HEC-RAS 6.6 (2D)\nSel 30 m, ditampilkan pada DTM 4 m',fontsize=9,color='0.3')
    fig.savefig(f'fr_{TAG}/{fi:04d}.png')
print(len(ti),'frame')
