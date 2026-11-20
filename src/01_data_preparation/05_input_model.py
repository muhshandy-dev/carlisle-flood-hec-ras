"""Langkah 5: tanggul (crest), input siap-pakai HEC-HMS dan HEC-RAS, kemiringan batas hilir."""
import json, numpy as np, pandas as pd, geopandas as gpd, rasterio
from shapely.geometry import LineString
from pathlib import Path

R = Path(__file__).resolve().parents[2]; P = R / "paket"; O3 = P / "03_hidromet"
OH, OR = P / "06_input_hec_hms", P / "07_input_hec_ras"
for d in (OH, OR): d.mkdir(exist_ok=True)
KEJ = ["desmond_2015", "ciara_2020", "jan_2005"]

# ---------- tanggul & tembok banjir EA (crest level) ----------
dfc = gpd.read_file(R / "probe" / "defences_raw.gpkg")
area = gpd.read_file(P / "02_medan" / "area_2d_usulan.gpkg").geometry.iloc[0]
dfc = dfc[dfc.intersects(area.buffer(200))].copy()
dfc["crest_mAOD"] = pd.to_numeric(dfc.effective_cl, errors="coerce").fillna(pd.to_numeric(dfc.actual_ucl, errors="coerce"))
with rasterio.open(P / "02_medan" / "dtm_2m_carlisle.tif") as d:
    def maks_dtm(g):
        if g.geom_type not in ("LineString", "MultiLineString"): return np.nan
        pts = [g.interpolate(t) for t in np.arange(0, g.length, 2.0)]
        v = np.array([x[0] for x in d.sample([(p.x, p.y) for p in pts])], float); v[v < -1000] = np.nan
        return float(np.nanmax(v)) if np.isfinite(v).any() else np.nan
    dfc["dtm2m_maks_mAOD"] = dfc.geometry.apply(maks_dtm)
dfc["selisih_crest_dtm_m"] = (dfc.crest_mAOD - dfc.dtm2m_maks_mAOD).round(2)
dfc["perlu_modifikasi_medan"] = dfc.selisih_crest_dtm_m > 0.3
dfc["mulai"] = pd.to_datetime(dfc.asset_start_date, errors="coerce").dt.date.astype(str)
keep = ["asset_id", "asset_name", "asset_sub_type", "water_course_name", "bank", "crest_mAOD", "dtm2m_maks_mAOD",
        "selisih_crest_dtm_m", "perlu_modifikasi_medan", "design_sop", "current_condition", "mulai", "geometry"]
dfc[[c for c in keep if c in dfc]].to_file(OR / "tanggul_tembok_banjir_EA.gpkg", driver="GPKG")
print("aset pertahanan:", len(dfc), "| perlu modifikasi medan (crest > DTM 2 m + 0,3 m):", int(dfc.perlu_modifikasi_medan.sum()))

# ---------- HEC-HMS: debit awal (baseflow) per sub-DAS per kejadian ----------
sub = pd.read_csv(P / "04_lahan_tanah_bangunan" / "parameter_awal_hec_hms_per_subdas.csv")
luas = dict(zip(sub.sub_das, sub.luas_km2))
peta_pos = {"Eden_TempleSowerby": "Eden_TempleSowerby", "Eamont_Udford": "Eamont_Udford", "Irthing_Greenholme": "Irthing_Greenholme",
            "Caldew_Cummersdale": "Caldew_Cummersdale", "Petteril_NewbigginBridge": "Petteril_NewbigginBridge"}
bf = []
for k in KEJ:
    q = pd.read_csv(O3 / f"debit_{k}_15menit_m3s.csv", index_col=0, parse_dates=True)
    q0 = q.iloc[:4].mean()
    est = {}
    for s, p in peta_pos.items():
        est[s] = q0.get(p, np.nan)
    # pos tanpa data -> dari debit spesifik (m3/s/km2) pos terdekat yang ada
    spes = np.nanmedian([est[s] / luas[s] for s in est if np.isfinite(est[s])])
    for s in est:
        if not np.isfinite(est[s]): est[s] = spes * luas[s]
    gc = q0.get("Eden_GreatCorby", np.nan)
    est["Eden_Tengah"] = max(gc - est["Eden_TempleSowerby"] - est["Eamont_Udford"], spes * luas["Eden_Tengah"]) if np.isfinite(gc) else spes * luas["Eden_Tengah"]
    sm = q0.get("Eden_Sheepmount", np.nan)
    sisa = sm - gc - est["Irthing_Greenholme"] - est["Caldew_Cummersdale"] - est["Petteril_NewbigginBridge"] if np.isfinite(sm) and np.isfinite(gc) else np.nan
    est["Lokal_Carlisle"] = sisa if np.isfinite(sisa) and sisa > 0 else spes * luas["Lokal_Carlisle"]
    for s, v in est.items():
        bf.append({"kejadian": k, "sub_das": s, "debit_awal_m3s": round(float(v), 2), "konstanta_resesi": 0.85,
                   "rasio_ambang_ke_puncak": 0.2})
pd.DataFrame(bf).to_csv(OH / "baseflow_awal_per_kejadian.csv", index=False)

# ---------- HEC-HMS: lag routing awal dari selisih waktu puncak observasi ----------
lag = []
for k in KEJ:
    q = pd.read_csv(O3 / f"debit_{k}_15menit_m3s.csv", index_col=0, parse_dates=True)
    t = {c: q[c].idxmax() for c in q if q[c].notna().mean() > 0.9}
    for a, b, nama in [("Eden_TempleSowerby", "Eden_GreatCorby", "R1_TempleSowerby_ke_GreatCorby"),
                       ("Eamont_Udford", "Eden_GreatCorby", "R2_Udford_ke_GreatCorby"),
                       ("Eden_GreatCorby", "Eden_Sheepmount", "R3_GreatCorby_ke_Sheepmount")]:
        if a in t and b in t:
            lag.append({"kejadian": k, "reach": nama, "selisih_puncak_jam": round((t[b] - t[a]).total_seconds() / 3600, 2)})
lag = pd.DataFrame(lag); lag.to_csv(OH / "lag_reach_dari_observasi.csv", index=False); print(lag.to_string())

# ---------- HEC-HMS: hujan per sub-DAS (format tempel: tanggal, jam, nilai) ----------
for k in KEJ:
    hy = pd.read_csv(O3 / f"hujan_{k}_15menit_per_subdas_mm.csv", index_col=0, parse_dates=True).drop(columns="DAS_Eden")
    out = hy.copy(); out.insert(0, "Tanggal", out.index.strftime("%d%b%Y").str.upper()); out.insert(1, "Jam", out.index.strftime("%H:%M"))
    out.to_csv(OH / f"hujan_{k}_15menit_untuk_HMS.csv", index=False)
    q = pd.read_csv(O3 / f"debit_{k}_15menit_m3s.csv", index_col=0, parse_dates=True)
    q.insert(0, "Tanggal", q.index.strftime("%d%b%Y").str.upper()); q.insert(1, "Jam", q.index.strftime("%H:%M"))
    q.to_csv(OH / f"debit_observasi_{k}_untuk_HMS.csv", index=False)
import shutil
shutil.copy(O3 / "hujan_rencana_ABM_24jam_per_jam.csv", OH / "hujan_rencana_ABM_24jam_per_jam.csv")

# ---------- HEC-RAS: hidrograf masuk observasi untuk kalibrasi hidrolika ----------
INFLOW = {"BC_Eden_hulu": "Eden_GreatCorby", "BC_Irthing": "Irthing_Greenholme", "BC_Caldew": "Caldew_Cummersdale",
          "BC_Petteril": "Petteril_NewbigginBridge"}
for k in KEJ:
    q = pd.read_csv(O3 / f"debit_{k}_15menit_m3s.csv", index_col=0, parse_dates=True)
    bc = pd.DataFrame({b: q[p] if p in q else np.nan for b, p in INFLOW.items()}, index=q.index).interpolate(limit=8)
    bc.insert(0, "jam_ke", ((bc.index - bc.index[0]).total_seconds() / 3600).round(2))
    bc.round(3).to_csv(OR / f"hidrograf_masuk_{k}_15menit.csv")
    print(k, "kosong per BC (%):", (bc.isna().mean() * 100).round(1).to_dict())

# ---------- HEC-RAS: kemiringan muka air untuk batas hilir (Normal Depth) ----------
lv = pd.read_csv(O3 / "tma_desmond_2015_15menit_mAOD.csv", index_col=0, parse_dates=True)
pos = gpd.read_file(O3 / "pos_tma_dan_debit.gpkg").set_index("nama")
a, b = "Sands Centre", "Sheepmount"
jarak = pos.loc[a].geometry.distance(pos.loc[b].geometry) * 1.25      # faktor kelokan sungai ~1,25
s_puncak = (lv[a].max() - lv[b].max()) / jarak
s_rendah = (lv[a].iloc[:8].mean() - lv[b].iloc[:8].mean()) / jarak
info = {"kemiringan_muka_air_puncak_Desmond": round(float(s_puncak), 5), "kemiringan_muka_air_awal": round(float(s_rendah), 5),
        "sumber": f"{a} -> {b}, jarak sungai ~{jarak/1000:.1f} km",
        "saran_batas_hilir": "Normal Depth dengan friction slope ~ nilai di atas (0,0005-0,001); uji sensitivitas +/-50%"}
json.dump(info, open(OR / "batas_hilir_normal_depth.json", "w"), indent=1); print(info)
