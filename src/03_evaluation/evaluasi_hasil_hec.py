"""
Evaluasi hasil HEC-RAS / HEC-HMS terhadap target kalibrasi (konsep simulasi banjir), peta bahaya, bangunan terdampak,
dan file siap 3D. Studi kasus: Carlisle (Storm Desmond 2015, Storm Ciara 2020, Januari 2005).

Cara pakai (Windows, Python 3.10+):
    pip install geopandas rasterio rashdf matplotlib
    python evaluasi_hasil_hec.py --paket "D:/carlisle_paket" --kejadian desmond_2015 ^
        --hdf "D:/ras/Carlisle.p01.hdf" --kedalaman "D:/ras/Desmond/Depth (Max).tif" ^
        --kecepatan "D:/ras/Desmond/Velocity (Max).tif" --wse "D:/ras/Desmond/WSE (Max).tif" --keluar hasil_desmond

Semua argumen selain --paket dan --kejadian opsional; bagian yang datanya tidak diberikan akan dilewati.
Di RAS Mapper, beri nama Reference Point PERSIS sama dengan nama pos TMA (mis. "Sheepmount", "Sands Centre")
dan Reference Line PERSIS sama dengan nama pos debit (mis. "Eden_Sheepmount").
--hms: CSV debit simulasi HEC-HMS, kolom pertama waktu, kolom lain bernama pos debit (mis. "Eden_GreatCorby").
"""
import argparse, json, warnings
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# ---- target dari dokumen konsep (usulan; ubah bila tim menyepakati angka lain) ----
TARGET = {"selisih_puncak_tma_m": 0.3, "selisih_puncak_debit_pct": 20, "selisih_waktu_puncak_jam": 1.0,
          "nse_cukup": 0.5, "nse_baik": 0.7, "csi": 0.5, "properti_tergenang_pct": 20}
PROPERTI_TERCATAT = {"desmond_2015": 2128}      # Laporan Section 19 Carlisle 2015 (EA & Cumbria County Council)
# kelas kedalaman (m). Nilai default sering dikutip dari Perka BNPB 2/2012; VERIFIKASI dengan lampiran peraturan.
KELAS_KEDALAMAN = [(0.05, 0.76, 1, "Rendah"), (0.76, 1.5, 2, "Sedang"), (1.5, 99, 3, "Tinggi")]
AMBANG_BASAH = 0.15                              # m; sel dianggap tergenang untuk CSI

def nse(o, s):
    m = o.notna() & s.notna()
    o, s = o[m], s[m]
    return float(1 - ((o - s) ** 2).sum() / ((o - o.mean()) ** 2).sum()) if len(o) > 10 else np.nan

def bandingkan(obs, sim, satuan, jenis):
    sim = sim.reindex(obs.index, method="nearest", tolerance=pd.Timedelta("10min"))
    m = obs.notna() & sim.notna()
    if m.sum() < 10:
        return None
    o, s = obs[m], sim[m]
    r = {"puncak_obs": round(o.max(), 3), "puncak_sim": round(s.max(), 3),
         "selisih_waktu_puncak_jam": round((s.idxmax() - o.idxmax()).total_seconds() / 3600, 2), "NSE": round(nse(o, s), 3)}
    if jenis == "tma":
        r["selisih_puncak"] = round(s.max() - o.max(), 3)
        r["lolos_puncak"] = abs(r["selisih_puncak"]) <= TARGET["selisih_puncak_tma_m"]
    else:
        r["selisih_puncak"] = round(100 * (s.max() - o.max()) / o.max(), 1)
        r["lolos_puncak"] = abs(r["selisih_puncak"]) <= TARGET["selisih_puncak_debit_pct"]
    r["lolos_waktu"] = abs(r["selisih_waktu_puncak_jam"]) <= TARGET["selisih_waktu_puncak_jam"]
    r["kategori_NSE"] = "baik" if r["NSE"] >= TARGET["nse_baik"] else ("cukup" if r["NSE"] >= TARGET["nse_cukup"] else "kurang")
    r["satuan"] = satuan
    return r

def grafik(obs, sim, judul, ylabel, f):
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(obs.index, obs.values, "k", lw=1.6, label="Observasi")
    ax.plot(sim.index, sim.values, "#d1495b", lw=1.4, label="Simulasi")
    ax.set_title(judul); ax.set_ylabel(ylabel); ax.grid(alpha=.3); ax.legend()
    fig.autofmt_xdate(); fig.tight_layout(); fig.savefig(f, dpi=110); plt.close(fig)

def baca_raster(path, ref_profile=None):
    with rasterio.open(path) as d:
        a = d.read(1, masked=True).filled(np.nan).astype("float32"); prof = d.profile; tr = d.transform; crs = d.crs
    a[a < -1000] = np.nan
    return a, tr, crs, prof

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paket", required=True); ap.add_argument("--kejadian", required=True)
    ap.add_argument("--hdf"); ap.add_argument("--kedalaman"); ap.add_argument("--kecepatan"); ap.add_argument("--wse")
    ap.add_argument("--hms"); ap.add_argument("--keluar", default="hasil_evaluasi")
    a = ap.parse_args()
    P, K, O = Path(a.paket), a.kejadian, Path(a.keluar); O.mkdir(parents=True, exist_ok=True)
    H3 = P / "03_hidromet"
    # kejadian nyata punya observasi; skenario kala ulang (mis. --kejadian T100) tidak
    f_t, f_q = H3 / f"tma_{K}_15menit_mAOD.csv", H3 / f"debit_{K}_15menit_m3s.csv"
    tma_obs = pd.read_csv(f_t, index_col=0, parse_dates=True) if f_t.exists() else pd.DataFrame()
    q_obs = pd.read_csv(f_q, index_col=0, parse_dates=True) if f_q.exists() else pd.DataFrame()
    if tma_obs.empty and q_obs.empty:
        print(f"Tidak ada data observasi untuk '{K}' -> hanya peta bahaya, bangunan, dan 3D yang dihitung.")
    pos = gpd.read_file(H3 / "pos_tma_dan_debit.gpkg")
    pos_ok = set(pos[pos.dalam_area_2d & pos.datum_tersedia].nama)
    hasil, ringkas = [], {}

    # ---------- 1. HEC-RAS: TMA di reference point dan debit di reference line ----------
    if a.hdf:
        from rashdf import RasPlanHdf
        with RasPlanHdf(a.hdf) as h:
            try:
                st = h.reference_points_stage(use_names=True); st.index = pd.to_datetime(st.index)
            except Exception as e:
                st = pd.DataFrame(); print("Reference point tidak ditemukan:", e)
            try:
                fl = h.reference_lines_flow(use_names=True); fl.index = pd.to_datetime(fl.index)
            except Exception as e:
                fl = pd.DataFrame(); print("Reference line tidak ditemukan:", e)
        for n in st.columns:
            if n in tma_obs and n in pos_ok:
                r = bandingkan(tma_obs[n], st[n], "mAOD", "tma")
                if r: hasil.append({"model": "HEC-RAS", "jenis": "TMA", "pos": n, **r})
                grafik(tma_obs[n], st[n], f"TMA {n} - {K}", "mAOD", O / f"tma_{n.replace(' ', '_')}.png")
        for n in fl.columns:
            if n in q_obs:
                r = bandingkan(q_obs[n], fl[n].abs(), "m3/s", "debit")
                if r: hasil.append({"model": "HEC-RAS", "jenis": "Debit", "pos": n, **r})
                grafik(q_obs[n], fl[n].abs(), f"Debit {n} - {K}", "m³/s", O / f"debit_ras_{n}.png")
        tak_cocok = [n for n in list(st.columns) + list(fl.columns) if n not in tma_obs and n not in q_obs]
        if tak_cocok: print("Nama reference tanpa pasangan observasi (cek ejaan):", tak_cocok)

    # ---------- 2. HEC-HMS: debit di pos ----------
    if a.hms:
        qs = pd.read_csv(a.hms, index_col=0, parse_dates=True)
        for n in qs.columns:
            if n in q_obs:
                r = bandingkan(q_obs[n], qs[n], "m3/s", "debit")
                if r: hasil.append({"model": "HEC-HMS", "jenis": "Debit", "pos": n, **r})
                grafik(q_obs[n], qs[n], f"Debit HMS {n} - {K}", "m³/s", O / f"debit_hms_{n}.png")
    if hasil:
        df = pd.DataFrame(hasil); df.to_csv(O / "kalibrasi_per_pos.csv", index=False)
        print(df[["model", "jenis", "pos", "puncak_obs", "puncak_sim", "selisih_puncak", "selisih_waktu_puncak_jam", "NSE",
                  "lolos_puncak", "lolos_waktu", "kategori_NSE"]].to_string())
        for (m, j), g in df.groupby(["model", "jenis"]):
            ringkas[f"{m} {j}: pos lolos target puncak"] = f"{int(g.lolos_puncak.sum())}/{len(g)}"
            ringkas[f"{m} {j}: pos lolos target waktu"] = f"{int(g.lolos_waktu.sum())}/{len(g)}"
            ringkas[f"{m} {j}: NSE median"] = round(float(g.NSE.median()), 3)

    # ---------- 3. Sebaran genangan: CSI terhadap outline resmi ----------
    if a.kedalaman:
        dep, tr, crs, prof = baca_raster(a.kedalaman)
        area = gpd.read_file(P / "02_medan" / "area_2d_usulan.gpkg").to_crs(crs)
        di_area = rasterize([(g, 1) for g in area.geometry], dep.shape, transform=tr).astype(bool)
        f_ol = P / "05_validasi" / f"genangan_resmi_{K}.gpkg"
        if f_ol.exists():
            ol = gpd.read_file(f_ol).to_crs(crs)
            obs = rasterize([(g, 1) for g in ol.geometry], dep.shape, transform=tr).astype(bool) & di_area
            sim = (np.nan_to_num(dep) >= AMBANG_BASAH) & di_area
            sungai = gpd.read_file(P / "02_medan" / "badan_air_os.gpkg").to_crs(crs)
            alur = rasterize([(g, 1) for g in sungai.geometry], dep.shape, transform=tr).astype(bool)
            obs, sim = obs & ~alur, sim & ~alur                 # alur sungai tidak dihitung
            A, B, C = int((obs & sim).sum()), int((~obs & sim).sum()), int((obs & ~sim).sum())
            px = abs(tr.a * tr.e) / 1e6
            ringkas.update({"CSI": round(A / max(A + B + C, 1), 3), "hit_rate": round(A / max(A + C, 1), 3),
                            "false_alarm_ratio": round(B / max(A + B, 1), 3), "luas_genangan_sim_km2": round((A + B) * px, 2),
                            "luas_genangan_resmi_km2": round((A + C) * px, 2), "lolos_CSI": A / max(A + B + C, 1) >= TARGET["csi"]})
            rgb = np.ones(dep.shape + (3,)); rgb[obs & sim] = (0.1, 0.5, 0.9); rgb[~obs & sim] = (0.9, 0.3, 0.3)
            rgb[obs & ~sim] = (0.95, 0.75, 0.2); rgb[~di_area] = (0.85, 0.85, 0.85)
            fig, ax = plt.subplots(figsize=(10, 7)); b = rasterio.transform.array_bounds(*dep.shape, tr)
            ax.imshow(rgb, extent=[b[0], b[2], b[1], b[3]]); ax.set_title(
                f"{K}: biru = cocok, merah = simulasi saja, kuning = resmi saja | CSI {ringkas['CSI']}")
            fig.tight_layout(); fig.savefig(O / "peta_kecocokan_genangan.png", dpi=120); plt.close(fig)

        # ---------- 4. Peta bahaya ----------
        kls = np.zeros(dep.shape, "uint8")
        for lo, hi, v, _ in KELAS_KEDALAMAN:
            kls[(dep >= lo) & (dep < hi)] = v
        p2 = {**prof, "dtype": "uint8", "nodata": 0, "count": 1, "compress": "deflate"}
        with rasterio.open(O / "bahaya_kelas_kedalaman.tif", "w", **p2) as d: d.write(kls, 1)
        px = abs(tr.a * tr.e) / 1e6
        luas = [{"kelas": n, "batas_m": f"{lo}-{hi}", "luas_km2": round(float((kls == v).sum() * px), 3)} for lo, hi, v, n in KELAS_KEDALAMAN]
        if a.kecepatan:
            vel, trv, _, _ = baca_raster(a.kecepatan)
            if vel.shape != dep.shape:
                v2 = np.full(dep.shape, np.nan, "float32")
                reproject(vel, v2, src_transform=trv, src_crs=crs, dst_transform=tr, dst_crs=crs, resampling=Resampling.bilinear)
                vel = v2
            dv = np.nan_to_num(dep) * np.nan_to_num(vel); D, V = np.nan_to_num(dep), np.nan_to_num(vel)
            h = np.zeros(dep.shape, "uint8")               # kelas H1-H6 (Australian Rainfall & Runoff, Smith dkk. 2014)
            basah = D >= 0.05
            h[basah] = 6
            h[basah & (dv <= 4) & (D <= 4) & (V <= 4)] = 5
            h[basah & (dv <= 1.0) & (D <= 2) & (V <= 2)] = 4
            h[basah & (dv <= 0.6) & (D <= 1.2) & (V <= 2)] = 3
            h[basah & (dv <= 0.6) & (D <= 0.5) & (V <= 2)] = 2
            h[basah & (dv <= 0.3) & (D <= 0.3) & (V <= 2)] = 1
            with rasterio.open(O / "bahaya_kedalaman_x_kecepatan_H1_H6.tif", "w", **p2) as d: d.write(h, 1)
            for v in range(1, 7):
                luas.append({"kelas": f"H{v}", "batas_m": "DxV", "luas_km2": round(float((h == v).sum() * px), 3)})
        pd.DataFrame(luas).to_csv(O / "luas_per_kelas_bahaya.csv", index=False)

        # ---------- 5. Bangunan terdampak + file 3D ----------
        bg = gpd.read_file(P / "04_lahan_tanah_bangunan" / "bangunan_dengan_tinggi.gpkg").to_crs(crs)
        lab = rasterize([(g, i + 1) for i, g in enumerate(bg.geometry)], dep.shape, transform=tr, dtype="int32",
                        all_touched=True)
        d_ = np.nan_to_num(dep); m = lab > 0
        mx = pd.Series(d_[m]).groupby(lab[m]).max()
        bg["kedalaman_maks_m"] = [round(float(mx.get(i + 1, 0.0)), 2) for i in range(len(bg))]
        bg["kelas_bahaya"] = "Kering"
        for lo, hi, v, n in KELAS_KEDALAMAN:
            bg.loc[(bg.kedalaman_maks_m >= lo) & (bg.kedalaman_maks_m < hi), "kelas_bahaya"] = n
        warna = {"Kering": "#9e9e9e", "Rendah": "#2e7d32", "Sedang": "#f9a825", "Tinggi": "#c62828"}
        bg["warna_3d"] = bg.kelas_bahaya.map(warna)
        bg.to_file(O / "bangunan_terdampak.gpkg", driver="GPKG")
        out3d = bg[["tinggi_m", "kedalaman_maks_m", "kelas_bahaya", "warna_3d", "fungsi_penting", "geometry"]].copy()
        out3d["elev_dasar_mAOD"] = bg.get("elev_tanah_min_mAOD")
        out3d.to_crs(4326).to_file(O / "bangunan_3d_wgs84.geojson", driver="GeoJSON")
        tergenang = bg[bg.kedalaman_maks_m >= AMBANG_BASAH]
        rb = tergenang.groupby("kelas_bahaya").size().to_dict()
        ringkas["bangunan_tergenang_(>=0,15 m)"] = int(len(tergenang)); ringkas["bangunan_per_kelas"] = rb
        ringkas["bangunan_penting_tergenang"] = tergenang.fungsi_penting.dropna().value_counts().to_dict()
        if K in PROPERTI_TERCATAT:
            ref = PROPERTI_TERCATAT[K]; s = 100 * (len(tergenang) - ref) / ref
            ringkas["properti_tercatat_laporan"] = ref
            ringkas["selisih_jumlah_bangunan_pct"] = round(s, 1)
            ringkas["catatan_properti"] = ("Satu gedung OS dapat berisi banyak properti (flat/ruko); bandingkan sebagai "
                                           "indikasi, bukan angka pasti.")
        if a.wse:
            w, trw, _, pw = baca_raster(a.wse)
            if w.shape == dep.shape:
                w[np.nan_to_num(dep) < 0.05] = np.nan
                with rasterio.open(O / "muka_air_3d_mAOD.tif", "w", **{**pw, "dtype": "float32", "nodata": -9999,
                                                                      "compress": "deflate"}) as d:
                    d.write(np.where(np.isnan(w), -9999, w).astype("float32"), 1)

    json.dump(ringkas, open(O / "ringkasan.json", "w"), indent=1, ensure_ascii=False, default=str)
    print(json.dumps(ringkas, indent=1, ensure_ascii=False, default=str))

if __name__ == "__main__":
    main()
