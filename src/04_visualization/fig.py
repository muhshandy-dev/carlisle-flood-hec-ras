import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import numpy as np, pandas as pd, geopandas as gpd, rasterio, h5py
from rasterio.enums import Resampling
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','Times New Roman','DejaVu Serif'],'font.size':9})
P=PAKET+''; F='fig/'
# ---------- Gambar 1: lokasi ----------
das=gpd.read_file(P+'01_das/das_eden.gpkg'); sub=gpd.read_file(P+'01_das/sub_das.gpkg').to_crs(27700); das=das.to_crs(27700)
area=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/area_2d_Carlisle.shp')
bc=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/garis_batas_BC.shp')
pos=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/reference_points_TMA.shp')
pq=gpd.read_file(P+'01_das/pos_debit.gpkg').to_crs(27700)
fig,ax=plt.subplots(1,2,figsize=(10,4.6),gridspec_kw={'width_ratios':[1,1.25]})
sub.plot(ax=ax[0],column=sub.columns[0],cmap='Pastel1',edgecolor='0.35',lw=0.6)
das.boundary.plot(ax=ax[0],color='k',lw=1.2); area.boundary.plot(ax=ax[0],color='#c0392b',lw=1.2)
pq.plot(ax=ax[0],color='#1f4e79',markersize=14,zorder=5)
for _,r in sub.iterrows():
    c=r.geometry.representative_point(); ax[0].annotate(str(r[sub.columns[0]]).replace('_','\n'),(c.x,c.y),ha='center',fontsize=6)
ax[0].set_title('(a) DAS Eden hingga Sheepmount dan 7 sub-DAS',fontsize=9); ax[0].set_xticks([]); ax[0].set_yticks([]); ax[0].set_xlabel(''); ax[0].set_ylabel('')
with rasterio.open(P+'02_medan/dtm_2m_carlisle_dengan_tanggul.tif') as s:
    f=8; z=s.read(1,out_shape=(s.height//f,s.width//f),resampling=Resampling.average,masked=True); b=s.bounds
gy,gx=np.gradient(z.filled(np.nan),16); hs=np.clip(0.55+0.5*(gx-gy),0,1)
ax[1].imshow(z,cmap='gist_earth',extent=[b.left,b.right,b.bottom,b.top],vmin=-30,vmax=160)
ax[1].imshow(hs,cmap='gray',alpha=0.35,extent=[b.left,b.right,b.bottom,b.top])
area.boundary.plot(ax=ax[1],color='k',lw=1); bc.plot(ax=ax[1],color='#c0392b',lw=3)
pos.plot(ax=ax[1],color='#1f4e79',markersize=12,zorder=5)
for _,r in bc.iterrows():
    c=r.geometry.centroid; ax[1].annotate(r['nama'].replace('BC_',''),(c.x,c.y),xytext=(4,-10) if 'hilir' in r['nama'] else (4,4),textcoords='offset points',fontsize=7,color='#c0392b')
ax[1].set_xlim(334500,349500); ax[1].set_ylim(550500,561500); ax[1].set_title('(b) Domain 2D, garis batas aliran, dan pos TMA',fontsize=9)
ax[1].set_xlabel('Easting (m, EPSG:27700)'); ax[1].set_ylabel('Northing (m)'); ax[1].tick_params(labelsize=7)
fig.tight_layout(); fig.savefig(F+'g1_lokasi.png',dpi=200); plt.close(fig)
# ---------- Gambar 2: kalibrasi HMS ----------
from hecdss import HecDss
d=HecDss(HMS+'Run_Desmond.dss')
def ts(p):
    t=d.get(p); return pd.Series(np.array(t.values),index=pd.to_datetime(t.times))
fig,ax=plt.subplots(1,3,figsize=(10,2.9)); hms={}
for a,(el,nm) in zip(ax,[('Eden_TempleSowerby','Eden di Temple Sowerby'),('Irthing_Greenholme','Irthing di Greenholme'),('J_Sheepmount','Eden di Sheepmount')]):
    s=ts(f'//{el}/FLOW/01Dec2015/15Minute/RUN:Run_Desmond/'); o=ts(f'//{el}/FLOW-OBSERVED/01Dec2015/15Minute/RUN:Run_Desmond/')
    o=o.where(o>-1e30); m=o.notna()&s.notna(); oo,ss=o[m],s[m]
    nse=1-((oo-ss)**2).sum()/((oo-oo.mean())**2).sum()
    hms[nm]=dict(Qobs=oo.max(),Qsim=ss.max(),dQ=100*(ss.max()-oo.max())/oo.max(),dt=(ss.idxmax()-oo.idxmax()).total_seconds()/3600,NSE=nse,bias=100*(ss.sum()-oo.sum())/oo.sum())
    a.plot(oo.index,oo,'k',lw=1.2,label='Observasi'); a.plot(ss.index,ss,color='#c0392b',lw=1.1,label='HEC-HMS')
    a.set_title(nm,fontsize=9); a.set_ylabel('Debit (m$^3$/s)'); a.grid(alpha=.3); a.tick_params(labelsize=7)
    a.xaxis.set_major_locator(matplotlib.dates.DayLocator()); a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m'))
ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(F+'g2_hms.png',dpi=200); plt.close(fig)
print(pd.DataFrame(hms).T.round(3))
# ---------- Gambar 3: RAS TMA (6 pos dalam kota) + debit ----------
f=h5py.File(RAS+'Carlisle.p01.hdf','r'); b='Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'
t=pd.to_datetime([x.decode() for x in f[b+'Time Date Stamp'][:]],format='%d%b%Y %H:%M:%S')
ws=f[b+'2D Flow Areas/Carlisle/Water Surface'][:]; g='Geometry/2D Flow Areas/Carlisle/'
xy=f[g+'Cells Center Coordinate'][:ws.shape[1]]; zmin=f[g+'Cells Minimum Elevation'][:ws.shape[1]]
obs=pd.read_csv(P+'03_hidromet/tma_desmond_2015_15menit_mAOD.csv',parse_dates=[0],index_col=0)
q=np.abs(f[b+'2D Flow Areas/Carlisle/Boundary Conditions/BC_Eden_hilir'][:][:,1])
qo=pd.read_csv(P+'03_hidromet/debit_desmond_2015_15menit_m3s.csv',parse_dates=[0],index_col=0)['Eden_Sheepmount']
w0=pd.Timestamp('2015-12-04 12:00')
names=['Sheepmount','Sands Centre','Linstock','Denton Holme','Botcherby Bridge']
fig,ax=plt.subplots(2,3,figsize=(10,5)); ax=ax.ravel()
m=t>=w0; qq=qo[(qo.index>=w0)&(qo.index<=t[-1])]; ax[0].plot(qq.index,qq,'k',lw=1.2,label='Observasi'); ax[0].plot(t[m],q[m],color='#1f4e79',lw=1.1,label='HEC-RAS')
ax[0].set_title('Debit Eden di Sheepmount',fontsize=9); ax[0].set_ylabel('m$^3$/s'); ax[0].legend(fontsize=7)
for a,nm in zip(ax[1:],names):
    r=pos[pos.nama==nm].iloc[0]; dd=np.hypot(xy[:,0]-r.geometry.x,xy[:,1]-r.geometry.y); c=np.where(dd<45)[0]; k=c[np.argmin(zmin[c])]
    o=obs[nm]; o=o[(o.index>=w0)&(o.index<=t[-1])]; a.plot(o.index,o,'k',lw=1.2); a.plot(t[m],ws[m,k],color='#1f4e79',lw=1.1)
    a.set_title('TMA '+nm,fontsize=9); a.set_ylabel('mAOD')
for a in ax: a.grid(alpha=.3); a.tick_params(labelsize=7); a.xaxis.set_major_locator(matplotlib.dates.DayLocator()); a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%d/%m'))
fig.tight_layout(); fig.savefig(F+'g3_ras_validasi.png',dpi=200); plt.close(fig)
# ---------- Gambar 4: kecocokan genangan ----------
from rasterio.features import rasterize
import matplotlib.patches as mp
D=RAS+'Demond/Depth (Max).Terrain.dtm_2m_carlisle_dengan_tanggul.tif'
with rasterio.open(D) as s:
    f8=2; dep=s.read(1,out_shape=(s.height//f8,s.width//f8),resampling=Resampling.nearest,masked=True).filled(0)
    tr=s.transform*s.transform.scale(f8,f8); crs=s.crs; bb=s.bounds
ol=gpd.read_file(P+'05_validasi/genangan_resmi_desmond_2015.gpkg').to_crs(crs)
riv=gpd.read_file(P+'02_medan/badan_air_os.gpkg').to_crs(crs)
ar=rasterize([(area.geometry[0],1)],dep.shape,transform=tr).astype(bool)
ob=rasterize([(g,1) for g in ol.geometry],dep.shape,transform=tr).astype(bool)&ar
rv=rasterize([(g,1) for g in riv.geometry],dep.shape,transform=tr).astype(bool)
sm=(dep>=0.15)&ar
rgb=np.ones(dep.shape+(3,)); rgb[~ar]=(0.9,0.9,0.9); rgb[ob&sm]=(0.12,0.47,0.71); rgb[sm&~ob]=(0.84,0.37,0.31); rgb[ob&~sm]=(0.95,0.74,0.22); rgb[rv&ar]=(0.6,0.6,0.6)
hull=ol.union_all().convex_hull.intersection(area.geometry[0])
fig,a=plt.subplots(figsize=(8,5.4)); a.imshow(rgb,extent=[bb.left,bb.right,bb.bottom,bb.top])
gpd.GeoSeries([hull],crs=crs).boundary.plot(ax=a,color='k',lw=0.8,ls='--')
a.legend(handles=[mp.Patch(color=(0.12,0.47,0.71),label='Tergenang pada model dan peta resmi'),mp.Patch(color=(0.84,0.37,0.31),label='Tergenang hanya pada model'),
 mp.Patch(color=(0.95,0.74,0.22),label='Tergenang hanya pada peta resmi'),mp.Patch(color=(0.6,0.6,0.6),label='Badan air permanen (tidak dinilai)'),mp.Patch(facecolor='white',edgecolor='0.6',label='Tidak tergenang'),
 matplotlib.lines.Line2D([],[],color='k',ls='--',lw=0.8,label='Zona cakupan peta resmi')],fontsize=7,loc='lower right')
a.set_xlim(334500,349500); a.set_ylim(550500,561500); a.tick_params(labelsize=7); a.set_xlabel('Easting (m)'); a.set_ylabel('Northing (m)')
fig.tight_layout(); fig.savefig(F+'g4_genangan.png',dpi=200); plt.close(fig)
