"""Emulator sederhana HEC-HMS (SCS-CN + SCS UH + Recession + Lag) untuk mencari titik awal kalibrasi,
plus koreksi hujan per sub-DAS ke total hujan DAS resmi NRFA (CEH-GEAR, data-type 'cdr')."""
import numpy as np, pandas as pd, requests, json
from scipy.optimize import minimize
from pathlib import Path
R = Path(__file__).resolve().parents[2]; H3 = R / "paket" / "03_hidromet"
DT = 0.25  # jam
SUB = pd.read_csv(R / "paket/04_lahan_tanah_bangunan/parameter_awal_hec_hms_per_subdas.csv").set_index("sub_das")
A = SUB.luas_km2.to_dict()
EV = {"desmond_2015": ("2015-12-03", "2015-12-09"), "ciara_2020": ("2020-02-07", "2020-02-13"), "jan_2005": ("2005-01-05", "2005-01-11")}
NRFA = {"Eden_TempleSowerby": "76005", "Eamont_Udford": "76003", "Irthing_Greenholme": "76008", "Caldew_Cummersdale": "76022",
        "Petteril_NewbigginBridge": "76020", "GC": "76017", "SM": "76007"}
AREA_N = {"76005": 616.4, "76003": 396.2, "76008": 334.6, "76022": 244.0, "76020": 145.82, "76017": 1373.0, "76007": 2286.5}

def cdr(sid):
    js = requests.get("https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series", params={"format": "json-object", "data-type": "cdr",
                      "station": sid}, timeout=120).json(); s = js["data-stream"]
    return pd.Series(s[1::2], index=pd.to_datetime(s[0::2]), dtype=float)
CDR = {k: cdr(v) for k, v in NRFA.items()}

faktor = {}
for ev, (a, b) in EV.items():
    h = pd.read_csv(H3 / f"hujan_{ev}_15menit_per_subdas_mm.csv", index_col=0, parse_dates=True)
    tot = {k: CDR[k][a:b].sum() for k in CDR}
    t = {k: tot[k] for k in ["Eden_TempleSowerby", "Eamont_Udford", "Irthing_Greenholme", "Caldew_Cummersdale", "Petteril_NewbigginBridge"]}
    aN = AREA_N
    t["Eden_Tengah"] = (tot["GC"] * aN["76017"] - tot["Eden_TempleSowerby"] * aN["76005"] - tot["Eamont_Udford"] * aN["76003"]) / (aN["76017"] - aN["76005"] - aN["76003"])
    t["Lokal_Carlisle"] = (tot["SM"] * aN["76007"] - tot["GC"] * aN["76017"] - tot["Irthing_Greenholme"] * aN["76008"]
                           - tot["Caldew_Cummersdale"] * aN["76022"] - tot["Petteril_NewbigginBridge"] * aN["76020"]) / (
                           aN["76007"] - aN["76017"] - aN["76008"] - aN["76022"] - aN["76020"])
    faktor[ev] = {s: round(float(t[s] / h[s].sum()), 3) for s in t}
json.dump(faktor, open(R / "kerja" / "faktor_hujan_nrfa.json", "w"), indent=1)
print(json.dumps(faktor, indent=1))

# ---------------- emulator ----------------
_xr = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5])
_yr = np.array([0, .03, .1, .19, .31, .47, .66, .82, .93, .99, 1, .99, .93, .86, .78, .68, .56, .46, .39, .33, .28, .207, .147, .107, .077, .055, .04, .029, .021, .015, .011, .005, 0])
def uh(area, lag_h):
    tp = lag_h + DT / 2; t = np.arange(0, 5 * tp + DT, DT)
    q = np.interp(t / tp, _xr, _yr); q = q / (q.sum() * DT * 3600) * area * 1e6 * 0.001   # m3/s per 1 mm
    return q
def subbasin(p, area, cn, ia, lag_h, q0, k, ratio):
    S = 25400 / cn - 254; P = np.cumsum(p); Pe = np.where(P > ia, (P - ia) ** 2 / (P - ia + S), 0)
    ex = np.diff(np.r_[0, Pe]); d = np.convolve(ex, uh(area, lag_h))[:len(p)]
    kq = k ** (DT / 24); qb = q0 * kq ** np.arange(len(p)); tot = d + qb
    ip = int(np.argmax(tot)); thr = ratio * tot[ip]
    j = next((i for i in range(ip, len(tot)) if tot[i] < thr), None)
    if j is not None:
        n_ = kq ** np.arange(len(tot) - j); tot[j:] = tot[j] * n_ + np.maximum(d[j:] - d[j] * n_, 0)
    return tot
def lagr(q, lag_h):
    n = int(round(lag_h / DT)); return np.r_[np.full(n, q[0]), q[:len(q) - n]] if n > 0 else q
def nse(o, s): m = ~np.isnan(o); return 1 - np.sum((o[m] - s[m]) ** 2) / np.sum((o[m] - o[m].mean()) ** 2)

ev = "desmond_2015"
h = pd.read_csv(H3 / f"hujan_{ev}_15menit_per_subdas_mm.csv", index_col=0, parse_dates=True)
q = pd.read_csv(H3 / f"debit_{ev}_15menit_m3s.csv", index_col=0, parse_dates=True).reindex(h.index)
bf = pd.read_csv(R / "paket/06_input_hec_hms/baseflow_awal_per_kejadian.csv").query("kejadian==@ev").set_index("sub_das").debit_awal_m3s
P = {s: h[s].values * faktor[ev][s] for s in faktor[ev]}
IA, RATIO = 10.0, 0.1
hasil = {}
def fit_head(s, obs):
    o = q[obs].values
    def f(x):
        cn, lag, k = x
        if not (40 <= cn <= 98 and 1 <= lag <= 20 and 0.7 <= k <= 0.98): return 9
        sim = subbasin(P[s], A[s], cn, IA, lag, bf[s], k, RATIO)
        vb = (np.nansum(sim) - np.nansum(o)) / np.nansum(o)
        return -(nse(o, sim)) + 2 * abs(vb)
    best = min((minimize(f, x0, method="Nelder-Mead", options={"maxiter": 600}) for x0 in [(70, 6, .85), (85, 5, .9), (60, 8, .8)]), key=lambda r: r.fun)
    cn, lag, k = best.x; sim = subbasin(P[s], A[s], cn, IA, lag, bf[s], k, RATIO); o = q[obs].values
    hasil[s] = {"CN": round(cn, 1), "lag_menit": int(round(lag * 60)), "k": round(k, 3), "NSE": round(nse(o, sim), 3),
                "puncak_sim": round(sim.max(), 1), "puncak_obs": round(np.nanmax(o), 1),
                "bias_vol_pct": round(100 * (np.nansum(sim) - np.nansum(o)) / np.nansum(o), 1),
                "selisih_jam": (h.index[np.argmax(sim)] - h.index[np.nanargmax(o)]).total_seconds() / 3600}
    return sim
SIM = {}
for s, o in [("Eden_TempleSowerby", "Eden_TempleSowerby"), ("Irthing_Greenholme", "Irthing_Greenholme"),
             ("Caldew_Cummersdale", "Caldew_Cummersdale"), ("Petteril_NewbigginBridge", "Petteril_NewbigginBridge")]:
    SIM[s] = fit_head(s, o)
# Eamont: tidak ada observasi 2015 -> pakai CN, k Temple Sowerby; lag awal dari GIS (dikalibrasi lagi dengan 2005)
eam = dict(cn=hasil["Eden_TempleSowerby"]["CN"], lag=SUB.loc["Eamont_Udford", "lag_SCS_menit"] / 60, k=hasil["Eden_TempleSowerby"]["k"])
# Great Corby: cari CN & lag Eden_Tengah + lag R1, R2
og = q["Eden_GreatCorby"].values
def f_gc(x):
    cn, lag, r1, r2 = x
    if not (40 <= cn <= 98 and 1 <= lag <= 20 and 0 <= r1 <= 12 and 0 <= r2 <= 12): return 9
    sE = subbasin(P["Eamont_Udford"], A["Eamont_Udford"], eam["cn"], IA, eam["lag"], bf["Eamont_Udford"], eam["k"], RATIO)
    sT = subbasin(P["Eden_Tengah"], A["Eden_Tengah"], cn, IA, lag, bf["Eden_Tengah"], eam["k"], RATIO)
    sim = lagr(SIM["Eden_TempleSowerby"], r1) + lagr(sE, r2) + sT
    return -nse(og, sim) + 2 * abs((np.nansum(sim) - np.nansum(og)) / np.nansum(og))
b = min((minimize(f_gc, x0, method="Nelder-Mead", options={"maxiter": 800}) for x0 in [(75, 7, 6, 5), (85, 5, 4, 4), (65, 9, 7, 6)]), key=lambda r: r.fun)
cnT, lagT, r1, r2 = b.x
sE = subbasin(P["Eamont_Udford"], A["Eamont_Udford"], eam["cn"], IA, eam["lag"], bf["Eamont_Udford"], eam["k"], RATIO)
sT = subbasin(P["Eden_Tengah"], A["Eden_Tengah"], cnT, IA, lagT, bf["Eden_Tengah"], eam["k"], RATIO)
GC = lagr(SIM["Eden_TempleSowerby"], r1) + lagr(sE, r2) + sT
hasil["Eamont_Udford"] = {"CN": eam["cn"], "lag_menit": int(round(eam["lag"] * 60)), "k": eam["k"], "catatan": "tanpa observasi 2015"}
hasil["Eden_Tengah"] = {"CN": round(cnT, 1), "lag_menit": int(round(lagT * 60)), "k": eam["k"]}
hasil["J_GreatCorby"] = {"NSE": round(nse(og, GC), 3), "puncak_sim": round(GC.max(), 1), "puncak_obs": round(np.nanmax(og), 1),
                         "bias_vol_pct": round(100 * (np.nansum(GC) - np.nansum(og)) / np.nansum(og), 1),
                         "selisih_jam": (h.index[np.argmax(GC)] - h.index[np.nanargmax(og)]).total_seconds() / 3600}
hasil["R1"] = {"lag_menit": int(round(r1 * 60))}; hasil["R2"] = {"lag_menit": int(round(r2 * 60))}
# Sheepmount: CN & lag Lokal + lag R3
osm = q["Eden_Sheepmount"].values
def f_sm(x):
    cn, lag, r3 = x
    if not (40 <= cn <= 98 and 1 <= lag <= 20 and 0 <= r3 <= 10): return 9
    sL = subbasin(P["Lokal_Carlisle"], A["Lokal_Carlisle"], cn, IA, lag, bf["Lokal_Carlisle"], 0.85, RATIO)
    sim = lagr(GC, r3) + SIM["Irthing_Greenholme"] + SIM["Caldew_Cummersdale"] + SIM["Petteril_NewbigginBridge"] + sL
    return -nse(osm, sim) + 2 * abs((np.nansum(sim) - np.nansum(osm)) / np.nansum(osm))
b = min((minimize(f_sm, x0, method="Nelder-Mead", options={"maxiter": 600}) for x0 in [(70, 7, 3.5), (85, 5, 2), (60, 9, 5)]), key=lambda r: r.fun)
cnL, lagL, r3 = b.x
sL = subbasin(P["Lokal_Carlisle"], A["Lokal_Carlisle"], cnL, IA, lagL, bf["Lokal_Carlisle"], 0.85, RATIO)
SM = lagr(GC, r3) + SIM["Irthing_Greenholme"] + SIM["Caldew_Cummersdale"] + SIM["Petteril_NewbigginBridge"] + sL
hasil["Lokal_Carlisle"] = {"CN": round(cnL, 1), "lag_menit": int(round(lagL * 60)), "k": 0.85}
hasil["R3"] = {"lag_menit": int(round(r3 * 60))}
hasil["J_Sheepmount"] = {"NSE": round(nse(osm, SM), 3), "puncak_sim": round(SM.max(), 1), "puncak_obs": round(np.nanmax(osm), 1),
                         "bias_vol_pct": round(100 * (np.nansum(SM) - np.nansum(osm)) / np.nansum(osm), 1),
                         "selisih_jam": (h.index[np.argmax(SM)] - h.index[np.nanargmax(osm)]).total_seconds() / 3600}
json.dump(hasil, open(R / "kerja" / "kalibrasi_emulator_desmond.json", "w"), indent=1)
for k_, v in hasil.items(): print(k_, v)
