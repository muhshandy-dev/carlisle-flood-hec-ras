# File model HEC

Folder ini berisi file teks proyek. File biner (terrain, land cover, HDF hasil, DSS) tidak disertakan karena besar dan dapat dibuat ulang.

## HEC-HMS 4.14 (`hec-hms/`)

| File | Isi |
|---|---|
| `Carlisle_Eden.hms` | Daftar komponen proyek |
| `Eden.basin` | 7 sub-DAS, 3 ruas penelusuran; parameter terkalibrasi (CN 60,6–64,6; lag 330–555 menit; Ia 10 mm; recession 0,85) |
| `Desmond.met`, `T100.met` | Hujan kejadian Desmond dan hujan rencana T100 (Specified Hyetograph) |
| `Desmond.control`, `T100.control` | Jendela waktu simulasi, interval 15 menit |
| `Carlisle_Eden.gage` | Pos hujan dan debit. Ubah path file DSS sesuai lokasi di komputer Anda |
| `Carlisle_Eden.run` | Run_Desmond dan Run_T100 |

## HEC-RAS 6.6 (`hec-ras/`)

| File | Isi |
|---|---|
| `Carlisle.prj` | Proyek |
| `Carlisle.g01` | Geometri: area 2D "Carlisle" (38,9 km²), mesh 30 m (42.634 sel), 5 garis BC |
| `Carlisle.u01` | Aliran Desmond 2015: debit terukur 15 menit + normal depth 0,0007 |
| `Carlisle.u02` | Aliran T100: hidrograf HEC-HMS yang disesuaikan ke frekuensi debit NRFA |
| `Carlisle.p01`, `Carlisle.p02` | Plan Desmond dan T100: diffusion wave, Δt 20 detik |
| `Carlisle.rasmap` | Konfigurasi RAS Mapper (proyeksi, terrain, land cover, peta hasil) |

Nilai Manning per kelas tutupan lahan: hutan 0,12; rumput 0,035; pertanian 0,04; terbangun 0,10; lahan terbuka 0,03; air 0,035; lahan basah 0,07; NoData 0,06.
