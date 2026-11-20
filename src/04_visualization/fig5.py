import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import rasterio, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.patches as mp
from matplotlib.colors import ListedColormap, BoundaryNorm
plt.rcParams.update({'font.family':'serif','font.serif':['Liberation Serif','Times New Roman'],'font.size':9})
Q=OUTDIR+'qgis_3d/'
fig,ax=plt.subplots(1,2,figsize=(10,4.1))
cm=ListedColormap(['#9ecae1','#4292c6','#08306b']); nm=BoundaryNorm([0.05,0.76,1.5,99],3)
for a,t,ttl in zip(ax,['desmond','t100'],['(a) Storm Desmond 2015','(b) Skenario kala ulang 100 tahun']):
    with rasterio.open(Q+f'kedalaman_{t}_4m.tif') as s: d=s.read(1,masked=True); b=s.bounds
    with rasterio.open(Q+f'terrain_banjir_{t}_4m.tif') as s: z=s.read(1,masked=True)
    gy,gx=np.gradient(z.filled(np.nan)); hs=np.clip(0.6+0.15*(gx-gy),0,1)
    e=[b.left,b.right,b.bottom,b.top]; a.imshow(hs,cmap='gray',extent=e,vmin=0,vmax=1); a.imshow(d,cmap=cm,norm=nm,extent=e)
    a.set_title(ttl,fontsize=9); a.set_xticks([]); a.set_yticks([])
fig.legend(handles=[mp.Patch(color=c,label=l) for c,l in zip(['#9ecae1','#4292c6','#08306b'],['< 0,76 m','0,76–1,5 m','> 1,5 m'])],loc='lower center',ncol=3,fontsize=8)
fig.tight_layout(rect=[0,0.07,1,1]); fig.savefig('fig/g5_banding.png',dpi=200)
