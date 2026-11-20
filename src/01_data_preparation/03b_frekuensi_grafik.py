"""Perbaikan langkah 3: uji Chi-kuadrat yang benar, GEV debit dengan L-moment, grafik hujan-debit, cek kelengkapan pos."""
import json, numpy as np, pandas as pd, requests
from scipy import stats
from scipy.special import gamma as G_
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

R = Path(__file__).resolve().parents[2]; O = R / "paket" / "03_hidromet"
T = np.array([2, 5, 10, 25, 50, 100])

# ---------- hujan: analisis frekuensi ----------
am = pd.read_csv(O / "hujan_harian_maks_tahunan_DAS.csv")
x = am.hujan_harian_maks_mm.values.astype(float); n = len(x)
mean, sd = x.mean(), x.std(ddof=1); cs = stats.skew(x, bias=False); ck = stats.kurtosis(x, fisher=False, bias=False)
lx = np.log10(x); ml, sl, csl = lx.mean(), lx.std(ddof=1), stats.skew(lx, bias=False)
a_g, b_g = mean - 0.5772 * sd * np.sqrt(6) / np.pi, sd * np.sqrt(6) / np.pi
DIST = {
    "Normal": (lambda v: stats.norm.cdf(v, mean, sd), lambda p: stats.norm.ppf(p, mean, sd)),
    "Log Normal": (lambda v: stats.norm.cdf(np.log10(v), ml, sl), lambda p: 10 ** stats.norm.ppf(p, ml, sl)),
    "Gumbel": (lambda v: stats.gumbel_r.cdf(v, a_g, b_g), lambda p: stats.gumbel_r.ppf(p, a_g, b_g)),
    "Log Pearson III": (lambda v: stats.pearson3.cdf(np.log10(v), csl, ml, sl),
                        lambda p: 10 ** stats.pearson3.ppf(p, csl, ml, sl)),
}
Gk = int(round(1 + 3.322 * np.log10(n))); dk = Gk - 2 - 1
chi_kr = stats.chi2.ppf(0.95, dk); d_kr = 1.36 / np.sqrt(n)
uji = []
for nm, (F, Q) in DIST.items():
    xs = np.sort(x); Fe = np.arange(1, n + 1) / (n + 1)
    dmax = float(np.max(np.abs(F(xs) - Fe)))
    # Chi-kuadrat: kelas dengan peluang teoretis sama (batas dari invers distribusi)
    batas = Q(np.linspace(0, 1, Gk + 1)[1:-1])
    obs = np.bincount(np.searchsorted(batas, x), minlength=Gk); exp = n / Gk
    chi = float(np.sum((obs - exp) ** 2 / exp))
    uji.append({"distribusi": nm, "D_maks": round(dmax, 3), "D_kritis_5pct": round(d_kr, 3), "lolos_SK": dmax < d_kr,
                "chi2": round(chi, 2), "chi2_kritis_5pct": round(chi_kr, 2), "lolos_chi2": chi < chi_kr})
uji = pd.DataFrame(uji)
uji["syarat_parameter"] = ["Cs≈0; Ck≈3", "Cs≈3Cv+Cv³", "Cs≈1,14; Ck≈5,4", "bebas (selain di atas)"]
cv = sd / mean
uji["parameter_terpenuhi"] = [abs(cs) < 0.3 and abs(ck - 3) < 0.5, abs(cs - (3 * cv + cv ** 3)) < 0.3,
                               abs(cs - 1.14) < 0.3 and abs(ck - 5.4) < 1.0, True]
ok = uji[uji.lolos_SK & uji.lolos_chi2 & uji.parameter_terpenuhi]
pilih = ok.sort_values("D_maks").distribusi.iloc[0] if len(ok) else uji.sort_values("D_maks").distribusi.iloc[0]
fr = pd.DataFrame({"kala_ulang_tahun": T})
for nm, (F, Q) in DIST.items():
    fr[nm.replace(" ", "") + "_mm"] = np.round(Q(1 - 1 / T), 1)
fr["dipilih"] = pilih; fr["hujan_rencana_24jam_mm"] = fr[pilih.replace(" ", "") + "_mm"]
fr.to_csv(O / "hujan_rencana_kala_ulang.csv", index=False); uji.to_csv(O / "uji_kecocokan_distribusi.csv", index=False)
json.dump({"n_tahun": n, "tahun_awal": int(am.tahun.min()), "tahun_akhir": int(am.tahun.max()), "rata2_mm": round(mean, 2),
           "sd_mm": round(sd, 2), "Cs": round(cs, 3), "Ck": round(ck, 3), "Cv": round(cv, 3), "kelas_chi2": Gk,
           "dk_chi2": dk, "distribusi_dipilih": pilih}, open(O / "statistik_hujan.json", "w"), indent=1)
print(uji.to_string()); print(fr.to_string())

# ABM + Mononobe 24 jam (dibuat ulang dari hujan terpilih)
abm = {}
for _, rw in fr.iterrows():
    R24 = rw.hujan_rencana_24jam_mm; t = np.arange(1, 25)
    blok = np.diff(np.r_[0, R24 / 24 * (24 / t) ** (2 / 3) * t])
    pola = np.zeros(24); tengah = 11
    for i, j in enumerate(np.argsort(-blok)):
        pola[tengah + ((i + 1) // 2) * (1 if i % 2 else -1)] = blok[j]
    abm[f"T{int(rw.kala_ulang_tahun)}_mm"] = pola.round(2)
pd.DataFrame(abm, index=pd.Index(range(1, 25), name="jam_ke")).to_csv(O / "hujan_rencana_ABM_24jam_per_jam.csv")

# ---------- debit: GEV (L-moment, Hosking 1990) + LP3 ----------
def lmom(v):
    v = np.sort(v); n = len(v); j = np.arange(1, n + 1)
    b0 = v.mean(); b1 = np.sum((j - 1) / (n - 1) * v) / n; b2 = np.sum((j - 1) * (j - 2) / ((n - 1) * (n - 2)) * v) / n
    l1, l2, l3 = b0, 2 * b1 - b0, 6 * b2 - 6 * b1 + b0
    return l1, l2, l3 / l2
def gev_lm(v):
    l1, l2, t3 = lmom(v); c = 2 / (3 + t3) - np.log(2) / np.log(3); k = 7.859 * c + 2.9554 * c ** 2
    a = l2 * k / ((1 - 2 ** -k) * G_(1 + k)); u = l1 - a * (1 - G_(1 + k)) / k
    return lambda F: u + a / k * (1 - (-np.log(F)) ** k)
pos = pd.read_csv(R / "kerja" / "pos_debit.csv"); qf = []
for _, p in pos.iterrows():
    js = requests.get("https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series", params={"format": "json-object",
                      "data-type": "amax-flow", "station": str(p.nrfa)}, timeout=120).json()
    v = np.array(js["data-stream"][1::2], float); v = v[np.isfinite(v) & (v > 0)]
    if len(v) < 10: continue
    q = gev_lm(v); lv = np.log10(v)
    lp3 = 10 ** stats.pearson3.ppf(1 - 1 / T, stats.skew(lv, bias=False), lv.mean(), lv.std(ddof=1))
    for t_, l_ in zip(T, lp3):
        qf.append({"pos": p.pos, "n_tahun": len(v), "kala_ulang_tahun": int(t_), "Q_GEV_m3s": round(float(q(1 - 1 / t_)), 1),
                   "Q_LP3_m3s": round(float(l_), 1), "Q_maks_tercatat_m3s": round(float(v.max()), 1)})
qf = pd.DataFrame(qf); qf.to_csv(O / "frekuensi_debit_puncak_NRFA.csv", index=False)
print(qf[qf.pos == "Eden_Sheepmount"].to_string())

# ---------- grafik hujan-debit & kelengkapan ----------
kel = []
for k in ["desmond_2015", "ciara_2020", "jan_2005"]:
    hy = pd.read_csv(O / f"hujan_{k}_15menit_per_subdas_mm.csv", index_col=0, parse_dates=True)
    q = pd.read_csv(O / f"debit_{k}_15menit_m3s.csv", index_col=0, parse_dates=True)
    lv = pd.read_csv(O / f"tma_{k}_15menit_mAOD.csv", index_col=0, parse_dates=True)
    for c in q: kel.append({"kejadian": k, "seri": f"debit {c}", "lengkap_pct": round(100 * q[c].notna().mean(), 1)})
    for c in lv: kel.append({"kejadian": k, "seri": f"TMA {c}", "lengkap_pct": round(100 * lv[c].notna().mean(), 1)})
    hj = hy.DAS_Eden.resample("1h").sum()
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    ax[0].bar(hj.index, hj.values, width=1 / 24, align="edge", color="#3b6fb6"); ax[0].invert_yaxis()
    ax[0].set_ylabel("Hujan rata-rata DAS\n(mm/jam)")
    q.plot(ax=ax[1], lw=1.3); ax[1].set_ylabel("Debit (m³/s)"); ax[1].grid(alpha=.3)
    ax[0].set_title(f"{k}: hujan rata-rata DAS Eden dan debit observasi (15 menit, UTC)")
    fig.tight_layout(); fig.savefig(O / f"grafik_{k}.png", dpi=110); plt.close(fig)
pd.DataFrame(kel).to_csv(O / "kelengkapan_data_kejadian.csv", index=False)
print(pd.DataFrame(kel).query("lengkap_pct < 90").to_string())
