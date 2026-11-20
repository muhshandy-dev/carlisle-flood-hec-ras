import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import sys, pandas as pd, numpy as np
from hecdss import HecDss
dss=sys.argv[1]; t0,t1=sys.argv[2],sys.argv[3]
d=HecDss(dss); paths=d.get_catalog().uncondensed_paths
peta={'BC_Eden_hulu':'J_GREATCORBY','BC_Irthing':'IRTHING_GREENHOLME','BC_Caldew':'CALDEW_CUMMERSDALE','BC_Petteril':'PETTERIL_NEWBIGGINBRIDGE'}
ser={}
for bc,el in peta.items():
    p=[x for x in paths if f'/{el}/FLOW/' in x.upper() and 'RUN_T100' in x.upper()]
    assert p, (el, [x for x in paths if el in x.upper()][:5])
    ts=d.get(p[0]); s=pd.Series(np.array(ts.values),index=pd.to_datetime(ts.times)); ser[bc]=s
    print(bc,p[0],'puncak %.1f m3/s pada %s'%(s.max(),s.idxmax()))
js=[x for x in paths if '/J_SHEEPMOUNT/FLOW/' in x.upper() and 'RUN_T100' in x.upper()]
if js:
    ts=d.get(js[0]); s=pd.Series(np.array(ts.values),index=pd.to_datetime(ts.times)); print('J_Sheepmount puncak %.1f pada %s'%(s.max(),s.idxmax()))
df=pd.read_csv('hidrograf_T100_skala_NRFA.csv',index_col=0,parse_dates=True).loc[t0:t1]; print(len(df), df.index[0], df.index[-1])
df.to_csv('hidrograf_T100_hms.csv')
mun=open(str(ROOT/'hec_models'/'hec-ras'/'Carlisle.u01'),newline='').read().split('\r\n')
i=[k for k,l in enumerate(mun) if l.startswith('Met Point Raster Parameters')][0]; tail=[l for l in mun[i:] if l!='']
SA='Carlisle'
loc=lambda bc: f"Boundary Location={'':16},{'':16},{'':8},{'':8},{'':16},{SA:<16},{'':16},{bc:<32},{'':32}"
L=["Flow Title=T100","Program Version=6.60","Use Restart= 0 ","Restart Filename="]
slopes={'BC_Eden_hulu':0.001,'BC_Irthing':0.001,'BC_Caldew':0.002,'BC_Petteril':0.002}
start=pd.Timestamp(t0).strftime('%d%b%Y').upper()+','+pd.Timestamp(t0).strftime('%H%M')
for bc,sl in slopes.items():
    v=df[bc].values; fv=[f"{x:8.2f}" for x in v]
    L+=[loc(bc),"Interval=15MIN",f"Flow Hydrograph= {len(v)} "]+[''.join(fv[k:k+10]) for k in range(0,len(fv),10)]
    L+=["Stage Hydrograph TW Check=0",f"Flow Hydrograph Slope= {sl} ","DSS Path=","Use DSS=False","Use Fixed Start Time=True",
        f"Fixed Start Date/Time={start}","Is Critical Boundary=False","Critical Boundary Flow="]
L+=[loc('BC_Eden_hilir'),"Friction Slope=0.0007,0"]+tail
open(OUTDIR+'Carlisle.u02','w',newline='').write('\r\n'.join(L)+'\r\n'); print('u02 ok', start)
