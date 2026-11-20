import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import h5py, numpy as np, pandas as pd, geopandas as gpd, json, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
hdf, out = sys.argv[1], sys.argv[2]
P=PAKET+''
f=h5py.File(hdf,'r'); b='Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/'
t=pd.to_datetime([s.decode() for s in f[b+'Time Date Stamp'][:]],format='%d%b%Y %H:%M:%S')
ws=f[b+'2D Flow Areas/Carlisle/Water Surface'][:]
g='Geometry/2D Flow Areas/Carlisle/'
xy=f[g+'Cells Center Coordinate'][:]; zmin=f[g+'Cells Minimum Elevation'][:]
n=ws.shape[1]; xy=xy[:n]; zmin=zmin[:n]
obs=pd.read_csv(P+'03_hidromet/tma_desmond_2015_15menit_mAOD.csv',parse_dates=[0],index_col=0)
pos=gpd.read_file(P+'07_input_hec_ras/shp_untuk_RAS_Mapper/reference_points_TMA.shp')
batas={'Greenholme','Cummersdale','Newbiggin Bridge'}
rows=[]; ser={}
w0=pd.Timestamp('2015-12-04 12:00')
for _,r in pos.iterrows():
    nm=r['nama']
    if nm not in obs.columns: continue
    d=np.hypot(xy[:,0]-r.geometry.x, xy[:,1]-r.geometry.y); c=np.where(d<45)[0]
    if len(c)==0: c=[d.argmin()]
    k=c[np.argmin(zmin[c])]
    sim=pd.Series(ws[:,k],index=t); o=obs[nm].reindex(t)
    m=(t>=w0)&o.notna().values
    if m.sum()<10: continue
    oo,ss=o[m],sim[m]
    nse=1-((oo-ss)**2).sum()/((oo-oo.mean())**2).sum()
    rows.append(dict(pos=nm,sungai=r['sungai'],dekat_batas=nm in batas,puncak_obs=round(oo.max(),2),puncak_sim=round(ss.max(),2),
        selisih_puncak_m=round(ss.max()-oo.max(),2),selisih_waktu_jam=round((ss.idxmax()-oo.idxmax()).total_seconds()/3600,1),NSE=round(nse,2),
        lolos_puncak_0_3m=abs(ss.max()-oo.max())<=0.3))
    ser[nm]=(oo,ss)
df=pd.DataFrame(rows); df.to_csv(out+'/tma_desmond.csv',index=False); print(df.to_string())
# debit keluar ~ Sheepmount
q=f[b+'2D Flow Areas/Carlisle/Boundary Conditions/BC_Eden_hilir'][:]
qs=pd.Series(np.abs(q[:,1]),index=t)
qo=pd.read_csv(P+'03_hidromet/debit_desmond_2015_15menit_m3s.csv',parse_dates=[0],index_col=0)['Eden_Sheepmount'].reindex(t)
m=(t>=w0)&qo.notna().values; oo,ss=qo[m],qs[m]
nse=1-((oo-ss)**2).sum()/((oo-oo.mean())**2).sum()
dq=dict(puncak_obs=round(oo.max(),1),puncak_sim=round(ss.max(),1),selisih_pct=round(100*(ss.max()-oo.max())/oo.max(),1),
    selisih_waktu_jam=round((ss.idxmax()-oo.idxmax()).total_seconds()/3600,1),NSE=round(nse,3))
print('Debit Sheepmount (aliran keluar BC_Eden_hilir):',dq); json.dump({k:float(v) for k,v in dq.items()},open(out+'/debit_sheepmount.json','w'),indent=1)
names=list(ser)[:12]; fig,ax=plt.subplots(4,3,figsize=(14,12)); ax=ax.ravel()
for i,nm in enumerate(names):
    oo,ss=ser[nm]; ax[i].plot(oo.index,oo,'k',lw=1.5,label='Observasi'); ax[i].plot(ss.index,ss,'#d1495b',lw=1.3,label='HEC-RAS')
    ax[i].set_title(nm,fontsize=10); ax[i].tick_params(labelsize=7); ax[i].grid(alpha=.3)
for j in range(len(names),12): ax[j].axis('off')
ax[0].legend(); fig.suptitle('TMA (mAOD) Storm Desmond 2015: observasi vs HEC-RAS 6.6'); fig.tight_layout(); fig.savefig(out+'/grafik_tma_desmond.png',dpi=110)
fig,a=plt.subplots(figsize=(9,4)); a.plot(oo.index,oo,'k',label='Observasi Sheepmount'); a.plot(ss.index,ss,'#d1495b',label='HEC-RAS (aliran keluar hilir)')
a.set_ylabel('Debit (m³/s)'); a.grid(alpha=.3); a.legend(); fig.tight_layout(); fig.savefig(out+'/grafik_debit_sheepmount.png',dpi=110)
