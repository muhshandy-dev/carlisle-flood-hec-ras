# Data

Seluruh data berasal dari sumber terbuka dan diunduh pada November 2026.

## `inputs/`

| File | Isi |
|---|---|
| `hidrograf_masuk_desmond_2015_15menit.csv` | Debit terukur di 4 garis batas hulu HEC-RAS (m³/s, interval 15 menit) |
| `hidrograf_T100_hms_asli.csv` | Hidrograf T100 dari HEC-HMS sebelum penyesuaian |
| `hidrograf_T100_skala_NRFA.csv` | Hidrograf T100 setelah limpasan langsung dikali 1,708 (puncak gabungan 1.563,4 m³/s) |
| `debit_desmond_2015_15menit_m3s.csv` | Debit terukur di pos-pos DAS Eden |
| `tma_desmond_2015_15menit_mAOD.csv` | Tinggi muka air terukur (mAOD) |
| `hujan_harian_maks_tahunan_DAS.csv` | Seri hujan harian maksimum tahunan rata-rata DAS, 1992–2025 |
| `uji_kecocokan_distribusi.csv` | Uji Smirnov-Kolmogorov dan Chi-kuadrat (gaya SNI 2415:2016) |
| `hujan_rencana_kala_ulang.csv` | Hujan rencana 24 jam per kala ulang |
| `frekuensi_debit_puncak_NRFA.csv` | Debit kala ulang dari GEV L-moment, data NRFA |
| `batas_hilir_normal_depth.json` | Kemiringan untuk batas hilir normal depth |

## `gis/` (EPSG:27700)

Area 2D, garis batas aliran (BC), titik TMA, garis acuan debit, breakline tanggul, dan wilayah perapatan alur sungai.

## Stasiun NRFA

| ID | Nama |
|---|---|
| 76007 | Eden at Sheepmount |
| 76017 | Eden at Great Corby |
| 76005 | Eden at Temple Sowerby |
| 76003 | Eamont at Udford |
| 76008 | Irthing at Greenholme |
| 76022 | Caldew at Cummersdale |
| 76020 | Petteril at Newbiggin Bridge |

## Link API yang dipakai skrip

Link berikut menampilkan data mentah (JSON) bila dibuka di browser.

- Pos TMA di Carlisle: https://environment.data.gov.uk/hydrology/id/stations.json?observedProperty=waterLevel&lat=54.895&long=-2.94&dist=12
- Pos hujan di sekitar DAS Eden: https://environment.data.gov.uk/hydrology/id/stations.json?observedProperty=rainfall&lat=54.75&long=-2.75&dist=70
- Datum pos TMA (contoh Sheepmount): https://environment.data.gov.uk/flood-monitoring/id/stations/765512/stageScale
- Debit puncak tahunan NRFA 76007: https://nrfaapps.ceh.ac.uk/nrfa/ws/time-series?format=json-object&data-type=amax-flow&station=76007
- Recorded Flood Outlines (GeoPackage nasional, ±84 MB): https://environment.data.gov.uk/api/file/download?fileDataSetId=ed73f2e8-a3c2-44db-952d-6e359c7c3987&fileName=Recorded_Flood_Outlines.gpkg.zip
- ESA WorldCover petak N54W003: https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N54W003_Map.tif
- SoilGrids WCS: https://maps.isric.org/

## Lisensi data

Environment Agency, Ordnance Survey, dan UKCEH: Open Government Licence v3.0. ESA WorldCover dan SoilGrids: CC BY 4.0.
