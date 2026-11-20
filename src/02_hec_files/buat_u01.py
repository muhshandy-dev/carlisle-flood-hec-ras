import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import pandas as pd
# Template: contoh proyek resmi HEC-RAS 6.6 (2D Unsteady Flow Hydraulics/Muncie/Muncie.u01)
mun=open(os.environ.get('MUNCIE_U01','Muncie.u01'),newline='').read().split('\r\n')
i=[k for k,l in enumerate(mun) if l.startswith('Met Point Raster Parameters')][0]
tail=[l for l in mun[i:] if l!='']
d=pd.read_csv(PAKET+'07_input_hec_ras/hidrograf_masuk_desmond_2015_15menit.csv',parse_dates=[0],index_col=0)
s=d.loc['2015-12-04 00:00':'2015-12-07 12:00']
SA='Carlisle'
def loc(bc): return f"Boundary Location={'':16},{'':16},{'':8},{'':8},{'':16},{SA:<16},{'':16},{bc:<32},{'':32}"
L=["Flow Title=Desmond_2015","Program Version=6.60","Use Restart= 0 ","Restart Filename="]
slopes={'BC_Eden_hulu':0.001,'BC_Irthing':0.001,'BC_Caldew':0.002,'BC_Petteril':0.002}
for bc,sl in slopes.items():
    v=s[bc].values; assert len(v)==337
    L+=[loc(bc),"Interval=15MIN",f"Flow Hydrograph= {len(v)} "]
    fv=[f"{x:8.2f}" for x in v]; assert all(len(t)==8 for t in fv)
    L+=[''.join(fv[k:k+10]) for k in range(0,len(fv),10)]
    L+=["Stage Hydrograph TW Check=0",f"Flow Hydrograph Slope= {sl} ","DSS Path=","Use DSS=False",
        "Use Fixed Start Time=True","Fixed Start Date/Time=04DEC2015,0000","Is Critical Boundary=False","Critical Boundary Flow="]
L+=[loc('BC_Eden_hilir'),"Friction Slope=0.0007,0"]
L+=tail
open('Carlisle.u01','w',newline='').write('\r\n'.join(L)+'\r\n')
print(len(L))
