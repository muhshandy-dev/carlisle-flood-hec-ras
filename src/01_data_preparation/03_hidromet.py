"""Langkah 3: hujan kejadian per sub-DAS, debit & TMA observasi, analisis frekuensi hujan (gaya SNI 2415:2016),
hujan rencana ABM-Mononobe, dan frekuensi debit puncak NRFA sebagai pembanding."""
import io, time, json, numpy as np, pandas as pd, geopandas as gpd, requests
from concurrent.futures import ThreadPoolExecutor
from scipy import stats
from shapely.geometry import Point, MultiPoint
from shapely.ops import voronoi_diagram
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

R = Path(__file__).resolve().parents[2]; K, OUT = R / "kerja", R / "paket" / "03_hidromet"
OUT.mkdir(parents=True, exist_ok=True); CRS = "EPSG:27700"
HYD = "https://environment.data.gov.uk/hydrology"
KEJADIAN = {  # jendela unduh (UTC); hujan mulai sebelum banjir
    "desmond_2015": ("2015-12-03", "2015-12-10"),
    "ciara_2020":   ("2020-02-07", "2020-02-14"),
    "jan_2005":     ("2005-01-05", "2005-01-12"),
}

def get(url, **p):
    for i in range(5):
        try:
            r = requests.get(url, params=p, timeout=300)
            if r.ok: return r
        except Exception:
            pass
        time.sleep(4 * (i + 1))
    return None

def bacaan(measure, a, b):
    r = get(f"{HYD}/id/measures/{measure}/readings.csv", **{"mineq-date": a, "max-date": b, "_limit": 2000000})
    if r is None or not r.text.strip(): return pd.Series(dtype=float)
    d = pd.read_csv(io.StringIO(r.text))
    if d.empty: return pd.Series(dtype=float)
    d = d[d.get("quality", "Good").fillna("Good") != "Missing"] if "quality" in d else d
    return pd.Series(d.value.values, index=pd.to_datetime(d.dateTime, utc=True)).sort_index()

sub = gpd.read_file(R / "paket" / "01_das" / "sub_das.gpkg")
das = gpd.read_file(R / "paket" / "01_das" / "das_eden.gpkg").geometry.iloc[0]

# ---------- 1. stasiun hujan di dalam DAS (+5 km) ----------
c = gpd.GeoSeries([das.centroid], crs=CRS).to_crs(4326).iloc[0]
items = get(f"{HYD}/id/stations.json", observedProperty="rainfall", lat=c.y, long=c.x, dist=70, _limit=500).json()["items"]
st = []
for it in items:
    try: p = gpd.GeoSeries([Point(float(it["long"]), float(it["lat"]))], crs=4326).to_crs(CRS).iloc[0]
    except Exception: continue
    if das.buffer(5000).contains(p):
        st.append({"id": it["notation"], "nama": it.get("label"), "geometry": p})
st = gpd.GeoDataFrame(st, crs=CRS)
print("stasiun hujan:", len(st))

# ---------- 2. hujan harian penuh + 15 menit kejadian ----------
def unduh_st(sid):
    harian = bacaan(f"{sid}-rainfall-t-86400-mm-qualified", "1990-01-01", "2026-11-01")
    ev = {k: bacaan(f"{sid}-rainfall-t-900-mm-qualified", a, b) for k, (a, b) in KEJADIAN.items()}
    return sid, harian, ev
with ThreadPoolExecutor(6) as ex:
    hasil = list(ex.map(unduh_st, st.id))
HARIAN = {s: h for s, h, _ in hasil if len(h)}
EV = {k: {s: e[k] for s, _, e in hasil if len(e[k])} for k in KEJADIAN}
st["tahun_harian"] = st.id.map(lambda s: len(HARIAN[s]) / 365.25 if s in HARIAN else 0).round(1)
st.to_file(OUT / "stasiun_hujan.gpkg", driver="GPKG")

# bobot Thiessen per sub-DAS (dihitung ulang tiap kombinasi stasiun yang tersedia)
_cache = {}
def bobot(ids, poly):
    key = (tuple(sorted(ids)), poly.wkt[:60])
    if key in _cache: return _cache[key]
    pts = st.set_index("id").loc[list(ids)].geometry
    if len(pts) == 1:
        w = pd.Series([1.0], index=pts.index)
    else:
        vor = voronoi_diagram(MultiPoint(list(pts)), envelope=poly.buffer(50000))
        w = {}
        for cell in vor.geoms:
            sid = pts.index[[cell.contains(p) for p in pts]][0]
            w[sid] = cell.intersection(poly).area
        w = pd.Series(w); w = w / w.sum()
    _cache[key] = w
    return w

def rata_area(df, poly):
    """df: kolom = stasiun, baris = waktu. Rata-rata Thiessen dengan stasiun yang tersedia per baris."""
    out = np.full(len(df), np.nan)
    pola = df.notna()
    for kombinasi, idx in pola.groupby(list(pola.columns)).groups.items():
        ids = [s for s, ada in zip(df.columns, kombinasi) if ada]
        if not ids: continue
        w = bobot(ids, poly)
        out[df.index.get_indexer(idx)] = (df.loc[idx, ids] * w[ids].values).sum(axis=1).values
    return pd.Series(out, index=df.index)

# ---------- 3. hujan kejadian per sub-DAS (15 menit) ----------
for k, (a, b) in KEJADIAN.items():
    df = pd.DataFrame(EV[k]).resample("15min").sum(min_count=1)
    df = df.loc[:, df.notna().mean() > 0.8]                       # stasiun dengan data >= 80%
    hy = pd.DataFrame({r.sub_das: rata_area(df, r.geometry) for r in sub.itertuples()})
    hy["DAS_Eden"] = rata_area(df, das)
    hy.index = hy.index.tz_convert(None); hy.index.name = "waktu_utc"
    hy.round(2).to_csv(OUT / f"hujan_{k}_15menit_per_subdas_mm.csv")
    print(k, "stasiun 15-mnt:", df.shape[1], "| total hujan DAS (mm):", round(hy.DAS_Eden.sum(), 1),
          "| maks sub-DAS:", hy.drop(columns="DAS_Eden").sum().round(0).to_dict())

# ---------- 4. debit & TMA observasi ----------
pos = pd.read_csv(K / "pos_debit.csv")
lvl = get(f"{HYD}/id/stations.json", observedProperty="waterLevel", lat=54.895, long=-2.94, dist=12, _limit=200).json()["items"]
tma_kota = {it["notation"]: (it["label"], it.get("riverName"), it.get("stationReference") or it.get("wiskiID"),
                            float(it["easting"]), float(it["northing"])) for it in lvl}
def datum(ref):
    r = get(f"https://environment.data.gov.uk/flood-monitoring/id/stations/{ref}/stageScale")
    try: return float(r.json()["items"]["datum"])
    except Exception: return np.nan
info_tma = []
for sid, (lab, riv, ref, x, y) in tma_kota.items():
    info_tma.append({"ea_id": sid, "nama": lab, "sungai": riv if isinstance(riv, str) else str(riv), "ref": ref,
                     "datum_mAOD": datum(ref) if ref else np.nan, "x": x, "y": y})
info_tma = pd.DataFrame(info_tma)
for _, p in pos.iterrows():
    if p.ea_id not in set(info_tma.ea_id):
        it = get(f"{HYD}/id/stations/{p.ea_id}.json").json()["items"][0]
        ref = it.get("stationReference") or it.get("wiskiID")
        info_tma.loc[len(info_tma)] = [p.ea_id, p.pos, it.get("riverName"), ref, datum(ref), p.x, p.y]
gpd.GeoDataFrame(info_tma, geometry=[Point(a, b) for a, b in zip(info_tma.x, info_tma.y)], crs=CRS).to_file(
    OUT / "pos_tma_dan_debit.gpkg", driver="GPKG")
ringkas = []
for k, (a, b) in KEJADIAN.items():
    Q = {p.pos: bacaan(f"{p.ea_id}-flow-i-900-m3s-qualified", a, b) for _, p in pos.iterrows()}
    L = {}
    for _, t in info_tma.iterrows():
        s = bacaan(f"{t.ea_id}-level-i-900-m-qualified", a, b)
        if len(s): L[f"{t.nama}"] = s + (t.datum_mAOD if np.isfinite(t.datum_mAOD) else 0)
    q = pd.DataFrame({n: s for n, s in Q.items() if len(s)}).resample("15min").mean()
    l = pd.DataFrame(L).resample("15min").mean()
    for d in (q, l): d.index = d.index.tz_convert(None); d.index.name = "waktu_utc"
    q.round(3).to_csv(OUT / f"debit_{k}_15menit_m3s.csv"); l.round(3).to_csv(OUT / f"tma_{k}_15menit_mAOD.csv")
    for n in q:
        if q[n].notna().any():
            ringkas.append({"kejadian": k, "pos": n, "debit_puncak_m3s": round(q[n].max(), 1), "waktu_puncak": q[n].idxmax()})
    for n in l:
        if l[n].notna().any():
            ringkas.append({"kejadian": k, "pos": n, "tma_puncak_mAOD": round(l[n].max(), 2), "waktu_puncak": l[n].idxmax()})
    # grafik
    hy = pd.read_csv(OUT / f"hujan_{k}_15menit_per_subdas_mm.csv", index_col=0, parse_dates=True)
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    ax[0].bar(hy.index, hy.DAS_Eden.resample("1h").sum().reindex(hy.index).fillna(0), width=0.04, color="#3b6fb6")
    ax[0].invert_yaxis(); ax[0].set_ylabel("Hujan DAS (mm/jam)")
    q.plot(ax=ax[1]); ax[1].set_ylabel("Debit (m³/s)"); ax[1].set_title(k)
    fig.tight_layout(); fig.savefig(OUT / f"grafik_{k}.png", dpi=110); plt.close(fig)
pd.DataFrame(ringkas).to_csv(OUT / "ringkasan_puncak_kejadian.csv", index=False)

# ---------- 5. hujan harian maksimum tahunan rata-rata DAS ----------
hd = pd.DataFrame(HARIAN).resample("1D").sum(min_count=1)
hd = hd.loc["1992":"2025"]
das_harian = rata_area(hd, das)
th = das_harian.groupby(das_harian.index.year)
lengkap = th.count() / th.size()
amax = th.max()[lengkap >= 0.9]
pd.DataFrame({"tahun": amax.index, "hujan_harian_maks_mm": amax.round(1).values,
              "tanggal": [das_harian[das_harian.index.year == y].idxmax().date() for y in amax.index]}
             ).to_csv(OUT / "hujan_harian_maks_tahunan_DAS.csv", index=False)
x = amax.values; n = len(x)
mean, sd = x.mean(), x.std(ddof=1); cs = stats.skew(x, bias=False); ck = stats.kurtosis(x, fisher=False, bias=False)
lx = np.log10(x); mean_l, sd_l, cs_l = lx.mean(), lx.std(ddof=1), stats.skew(lx, bias=False)
T = np.array([2, 5, 10, 25, 50, 100])
# Gumbel (metode momen); Log Normal; Log Pearson III
K_g = -(np.sqrt(6) / np.pi) * (0.5772 + np.log(np.log(T / (T - 1))))
R_g = mean + K_g * sd
R_ln = 10 ** (mean_l + stats.norm.ppf(1 - 1 / T) * sd_l)
R_lp3 = 10 ** (mean_l + stats.pearson3.ppf(1 - 1 / T, cs_l) * sd_l)
R_n = mean + stats.norm.ppf(1 - 1 / T) * sd
# uji Smirnov-Kolmogorov dan Chi-kuadrat
a_g, b_g = mean - 0.5772 * sd * np.sqrt(6) / np.pi, sd * np.sqrt(6) / np.pi
cdf = {
    "Normal": lambda v: stats.norm.cdf(v, mean, sd),
    "Log Normal": lambda v: stats.norm.cdf(np.log10(v), mean_l, sd_l),
    "Gumbel": lambda v: stats.gumbel_r.cdf(v, a_g, b_g),
    "Log Pearson III": lambda v: stats.pearson3.cdf(np.log10(v), cs_l, mean_l, sd_l),
}
G = int(1 + 3.322 * np.log10(n)); uji = []
for nm, F in cdf.items():
    xs = np.sort(x); Fe = np.arange(1, n + 1) / (n + 1)
    dmax = np.max(np.abs(F(xs) - Fe)); dkr = 1.36 / np.sqrt(n)
    tepi = np.quantile(x, np.linspace(0, 1, G + 1)); tepi[0], tepi[-1] = -np.inf, np.inf
    obs = np.histogram(x, tepi)[0]; exp = n * np.diff(np.clip(F(np.where(np.isinf(tepi), np.sign(tepi) * 1e9, tepi)), 0, 1))
    chi = float(np.sum((obs - exp) ** 2 / np.maximum(exp, 1e-9))); dk = G - 2 - 1
    uji.append({"distribusi": nm, "D_maks": round(dmax, 3), "D_kritis_5pct": round(dkr, 3), "lolos_SK": dmax < dkr,
                "chi2": round(chi, 2), "chi2_kritis_5pct": round(stats.chi2.ppf(0.95, max(dk, 1)), 2),
                "lolos_chi2": chi < stats.chi2.ppf(0.95, max(dk, 1))})
uji = pd.DataFrame(uji)
syarat = pd.DataFrame([
    {"distribusi": "Normal", "syarat": "Cs ≈ 0; Ck ≈ 3", "terpenuhi": abs(cs) < 0.3 and abs(ck - 3) < 0.5},
    {"distribusi": "Log Normal", "syarat": "Cs ≈ 3Cv + Cv³", "terpenuhi": abs(cs - (3 * sd / mean + (sd / mean) ** 3)) < 0.3},
    {"distribusi": "Gumbel", "syarat": "Cs ≈ 1,14; Ck ≈ 5,4", "terpenuhi": abs(cs - 1.14) < 0.3 and abs(ck - 5.4) < 1.0},
    {"distribusi": "Log Pearson III", "syarat": "selain di atas", "terpenuhi": True},
])
uji = uji.merge(syarat, on="distribusi")
pilih = uji[uji.lolos_SK & uji.lolos_chi2 & uji.terpenuhi]
pilih = pilih.distribusi.iloc[0] if len(pilih) else uji.sort_values("D_maks").distribusi.iloc[0]
fr = pd.DataFrame({"kala_ulang_tahun": T, "Normal_mm": R_n.round(1), "LogNormal_mm": R_ln.round(1),
                   "Gumbel_mm": R_g.round(1), "LogPearsonIII_mm": R_lp3.round(1)})
fr["dipilih"] = pilih
kol = {"Normal": "Normal_mm", "Log Normal": "LogNormal_mm", "Gumbel": "Gumbel_mm", "Log Pearson III": "LogPearsonIII_mm"}[pilih]
fr["hujan_rencana_24jam_mm"] = fr[kol]
fr.to_csv(OUT / "hujan_rencana_kala_ulang.csv", index=False)
uji.to_csv(OUT / "uji_kecocokan_distribusi.csv", index=False)
json.dump({"n_tahun": int(n), "tahun_awal": int(amax.index.min()), "tahun_akhir": int(amax.index.max()),
           "rata2": round(mean, 2), "sd": round(sd, 2), "Cs": round(cs, 3), "Ck": round(ck, 3), "Cv": round(sd / mean, 3),
           "dipilih": pilih}, open(OUT / "statistik_hujan.json", "w"), indent=1)
print(uji.to_string()); print(fr.to_string())

# ---------- 6. hujan rencana per jam: Alternating Block Method + intensitas Mononobe (24 jam) ----------
abm = {}
for _, rw in fr.iterrows():
    R24 = rw.hujan_rencana_24jam_mm; t = np.arange(1, 25)
    I = R24 / 24 * (24 / t) ** (2 / 3); kum = I * t; blok = np.diff(np.r_[0, kum])
    urut = np.argsort(-blok); pola = np.zeros(24); tengah = 11
    for i, j in enumerate(urut):
        pos_ = tengah + ((i + 1) // 2) * (1 if i % 2 else -1)
        pola[pos_] = blok[j]
    abm[f"T{int(rw.kala_ulang_tahun)}_mm"] = pola.round(2)
abm = pd.DataFrame(abm, index=pd.Index(range(1, 25), name="jam_ke"))
abm.to_csv(OUT / "hujan_rencana_ABM_24jam_per_jam.csv")

# ---------- 7. frekuensi debit puncak NRFA (pembanding hasil HEC-HMS) ----------
qf = []
for _, p in pos.iterrows():
    js = get("https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series", format="json-object", **{"data-type": "amax-flow"},
             station=str(p.nrfa)).json()
    v = np.array(js["data-stream"][1::2], float); v = v[np.isfinite(v)]
    if len(v) < 10: continue
    c_, loc, sc = stats.genextreme.fit(v)
    lv = np.log10(v); lp3 = 10 ** (lv.mean() + stats.pearson3.ppf(1 - 1 / T, stats.skew(lv, bias=False)) * lv.std(ddof=1))
    for t_, g_, l_ in zip(T, stats.genextreme.ppf(1 - 1 / T, c_, loc, sc), lp3):
        qf.append({"pos": p.pos, "n_tahun": len(v), "kala_ulang_tahun": t_, "Q_GEV_m3s": round(g_, 1), "Q_LP3_m3s": round(l_, 1)})
pd.DataFrame(qf).to_csv(OUT / "frekuensi_debit_puncak_NRFA.csv", index=False)
print(pd.DataFrame(qf).query("pos=='Eden_Sheepmount'").to_string())
