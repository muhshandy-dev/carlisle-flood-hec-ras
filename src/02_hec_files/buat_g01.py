import os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PAKET = os.environ.get('PAKET_DIR', str(ROOT / 'paket')) + '/'          # keluaran skrip 01_data_preparation
RAS = os.environ.get('RAS_DIR', 'C:/RAS_Carlisle66') + '/'              # folder proyek HEC-RAS
HMS = os.environ.get('HMS_DIR', str(Path.home() / 'Documents' / 'Carlisle_Eden')) + '/'  # folder proyek HEC-HMS
OUTDIR = os.environ.get('OUT_DIR', str(ROOT / 'out')) + '/'
import geopandas as gpd
from shapely.geometry import Point
src=RAS+'Carlisle.g01'
shp=PAKET+'07_input_hec_ras/shp_untuk_RAS_Mapper/'
lines=open(src,newline='').read().split('\r\n')
if lines[-1]=='': lines=lines[:-1]
p=gpd.read_file(shp+'area_2d_Carlisle.shp').geometry[0]
ring=list(p.exterior.coords)
if ring[0]!=ring[-1]: ring.append(ring[0])
f=lambda v: f"{v:16.4f}"
assert all(len(f(v))==16 for xy in ring for v in xy)
NAME='Carlisle'
blk=[f"Storage Area={NAME:<16},,", f"Storage Area Surface Line= {len(ring)} "]
blk+= [f(x)+f(y)+' '*16 for x,y in ring]
blk+=["Storage Area Type= 1 ","Storage Area Area=","Storage Area Min Elev=","Storage Area Is2D=-1",
"Storage Area Point Generation Data=,,30,30","Storage Area Mannings=0.06",
"2D Cell Volume Filter Tolerance=0.003","2D Cell Minimum Area Fraction=0.01","2D Face Profile Filter Tolerance=0.003",
"2D Face Area Elevation Profile Filter Tolerance=0.003","2D Face Area Elevation Conveyance Ratio=0.02",
"2D Face Area Laminar Depth=0.2","2D Multiple Face Mann n=-1","2D Composite LC=-1",""]
bc=gpd.read_file(shp+'garis_batas_BC.shp')
for _,r in bc.iterrows():
    c=list(r.geometry.coords); m=r.geometry.interpolate(0.5,normalized=True)
    d=max(p.exterior.distance(Point(q)) for q in [c[0],c[-1]])
    print(r['nama'],len(c),'pts, len %.0f m'%r.geometry.length,'end-dist-to-perim %.1f'%d, 'crosses' ,r.geometry.intersects(p.exterior))
    blk+=[f"BC Line Name={r['nama']:<32}",f"BC Line Storage Area={NAME:<16}",
      f"BC Line Start Position= {c[0][0]} , {c[0][1]} ",f"BC Line Middle Position= {m.x} , {m.y} ",
      f"BC Line End Position= {c[-1][0]} , {c[-1][1]} ",f"BC Line Arc= {len(c)} "]
    vals=[f(v) for xy in c for v in xy]
    blk+=[''.join(vals[i:i+4]) for i in range(0,len(vals),4)]
    blk+=["BC Line Text Position= 1.79769313486232E+308 , 1.79769313486232E+308 "]
xmin,ymin,xmax,ymax=p.buffer(500).bounds
out=[]
for i,l in enumerate(lines):
    if l.startswith('Viewing Rectangle='): l=f"Viewing Rectangle= {xmin} , {xmax} , {ymax} , {ymin} "
    out.append(l)
    if i==3: out+=blk
open('Carlisle.g01','w',newline='').write('\r\n'.join(out)+'\r\n')
print(len(ring),'perimeter pts')
