# data_map.md — data sources for nurture-oe

**Distilled from** `appendix_A_nurture_datasets.md` (2026-10-09) and `index.md` v1.10 (Project Knowledge), plus a directory inspection of the local data stores (2026-10-05). §1 re-synced to the PK snapshot on 2026-10-14. Project Knowledge is authoritative; this file is a working summary. Register IDs (D/G/O) refer to `index.md` §11.

**PK snapshot.** `docs/pk_snapshot/` holds read-only copies of `index-1-10.md`, `appendix_A_nurture_datasets.md` and `handoff.md`, made 2026-10-14. The originals are on OneDrive; never write there. Facts found in this repo after 2026-10-09 are in `docs/findings_for_pk.md`, not here, until they are folded into PK (next: `index.md` v1.11).

Paths are keys in `config/paths.yaml`. "Unverified" means nobody has checked it against the file yet, so check it before relying on it.

---

## 1. NURTURE example flight — 2026-02-17 (`nurture_example`)

Platform: NASA LaRC G-III N520NA. HALO and CloudCube share the airframe. NURTURE science flight #17, out of Goose Bay (MetNav header).

| File | Product | Coverage (UTC) |
|---|---|---|
| `nurture-HALO_G3_20260217_R0_L1.h5` | HALO Aerosol + Water Vapor L2 (HSRL/DIAL profiles) | 14:41:48–19:25:00 |
| `NURTURE-HALO-CloudAndSurface_G3_20260217_R0_L1.h5` | HALO Cloud & Surface L2 | same span |
| `NURTURE-CLOUDCUBE_G3_20260217_R0_S20260217a151300_E20260217a152010_KaWG.nc` | CloudCube L2, Ka/W/G | 15:13:00–15:20:10 |
| `NURTURE-CLOUDCUBE_G3_20260217_R0_S20260217a152158_E20260217a185937_KaWG.nc` | CloudCube L2, Ka/W/G (4.7 GB) | 15:21:58–18:59:37 — **STM segment granule** |
| `NURTURE-CLOUDCUBE_G3_20260217_R0_S20260217a190000_E20260217a192625_KaWG.nc` | CloudCube L2, Ka/W/G | 19:00:00–19:26:25 |
| `NURTURE-MetNav_G3_20260217_R0.ict` | G-III navigation/met (ICARTT, 1 Hz) | flight time |

In the HALO filenames, `L1` is the **flight number** (`flight_number_day`), not a processing level. Both HALO files are L2.

**The two HALO products are complements, not substitutes** (Appendix A §A.10.4, §A.10.6). Aerosol/WV is the source of lidar and DIAL profiles. Only Cloud & Surface carries the phase mask (incl. HOI), multilayer boundaries, the liquid extinction profile and the surface fields. They share one flight, time base and `z` grid.

### 1.1 HALO Aerosol + Water Vapor L2 (`halo_aerowv`) — the primary lidar/DIAL source

**Structure** (01-24 header dump, appendix §A.10; 02-17 values checked 2026-10-06, §A.10.5)
- HDF5, CF-1.10, R0. Root dims `time` and `z`. 02-17: **1,689 profiles × 845 bins**; 01-24 had 1,416 × 859. Root variables `lat`, `lon`, `time` (s since flight-date 00 UTC), `z` (m).
- `z` runs from −300 m to 12,351 m at 14.99 m, **ascending**, and is identical to C&S `CloudMask/Altitude`. It is a 15 m **sampling** grid; effective resolution is 30 m (3rd-order Bessel range filter, then interpolation `[O37]`). The header's `z:positive = "up "` has a trailing space.
- Groups `DataProducts`, `DataUncertainty`, `Nav_Data`, `State`, `UserInput` use `phony_dim_*` dimensions. Map them back to `time`/`z` by length and verify. `DataProducts/Altitude` exists on 02-17.
- Each 10 s profile is the **mean of twenty 0.5 s C&S profiles**, about 2.4 km along track `[O36]`. `Nav_Data/gps_time` (hours) equals `time` (seconds).
- Undefined within 0.5–2.5 km below the aircraft (product-dependent; `532_ext` 2.5 km).

**Name and attribute traps**
- Variable names start with digits (`532_ext`). Access them with explicit strings.
- Attribute names contain trailing spaces (`"Description "`, `"Units "`). Several variables carry both `Units ` and `units`, sometimes disagreeing (`Angstrom_Spherical`: `"sr"` vs `"1"`). Use CF `units` and warn on disagreement.
- `_FillValue = NaN` everywhere; no `mask_low`, no fill-reason flag `[G23]`.
- `State/Pressure` is in **atm**. `State_Type` scalar: 1 balloon sonde, 2 MERRA-2, 3 LaRC O₃ sonde, 4 dropsonde, 5 HRRR. 1064 products list dropsondes as state source, 532 products MERRA-2 `[G23]`.
- `UserInput/DEM_altitude`: long_name GMTED2010, global attribute GLOBE (cosmetic).
- Header ↔ README conflicts: appendix §A.10.3 (16 items).

**Roles and uncertainty** (appendix §A.10.2, §A.10.4, §A.10.6; `[D21]`, `[D22]`, `[G22]`)

| Variable(s) | Res. | `_unc`? | Role |
|---|---|---|---|
| `532_total_attn_bsc` (`1064_…`) | 10 s, 15 m | none | $\ln\beta'$ rows in liquid / optically thick cloud `[D21]` |
| `532_dep` (`1064_dep`) | 10 s | screened only | $\delta_v$ rows with $\beta'$ |
| `532_bsc`, `532_bsr` (+ screened) | 10 s | screened only | $\ln\beta_p$ rows where HSRL-valid (above cloud, thin ice) `[D21]`, switch `[O33]`; `bsr` = $\beta_p/\beta_m$ |
| `1064_bsc`, `1064_bsr` | 10 s | screened only | Fernald, normalized to 532 ⇒ not independent of 532 |
| `532_ext`, `532_Sa` | 60 s, 300 m running mean | yes | not rows alongside $\beta_p$ + $\beta'$; excluded `[O33]` |
| `1064_ext`, `532_extinction_from_backscatter` | | none | derived; exclude |
| `h2o_mmr_v` | 60 s, 315–705 m | yes | $\ln q_v$ rows out of cloud; in-liquid values not credible `[O32]` |
| `Aerosol_ID` | 10 s | — | class 1 = Ice (9 classes in header) |
| `/State` | (time, z) | none | MERRA-2 $T$, $p$ for $F(\mathbf{x})$ and $k_{\text{gas}}$ |

02-17 fill fractions in cloud (§A.10.5(b)): unscreened 532 $\beta'$/$\delta$/$\beta_p$ 0.91 liquid, 0.80 ice; 1064 0.17 / 0.48; screened fields and their `_unc` 0.055 / 0.50 (filled on identical pixels); `532_ext` 0.04 / 0.46. Cloud pixels there were defined from ≥50 % of the twenty C&S sub-profiles (liquid = class 0, ice = classes 1–2), >2.5 km below the aircraft: 10,395 liquid, 62,567 ice.

Running means are kept on the full grid ⇒ adjacent rows correlated `[O29]`. Readers return them as-is.

**Checks for task (a)** — done 2026-10-06 (in chat; to be reproduced as regression tests): variable diff vs 01-24 (partly), in-cloud population, `_cloud_screened_unc` masking, `z` ascending, C&S time-base match. **Still open:**
1. Full header dump and the complete diff against the 01-24 variable list.
2. Value of `State_Type` (expected 2 = MERRA-2). Units of every `State` variable.
3. `imu_pitch`/`imu_roll` during cloud legs: 10 s averages or samples? Compare with MetNav. `[G21]`, `[O28]`
4. NaN-pattern census per variable: near-field band, below-surface, above-cloud truncation of `h2o_mmr_v`. `[G23]`

*Task (a) run 2026-10-14: results in `docs/reports/task_a_readers.md`; new facts in `docs/findings_for_pk.md` until folded into v1.11.*

### 1.2 HALO Cloud & Surface L2 (`halo_cloudsurface`)

Appendix §A.2. 02-17 granule: **33,765 profiles at 0.5 s**, 845 bins identical to Aerosol/WV `z` (§A.2.4).

- **Only two profile datasets:** `DataProducts/CloudMask/Type` (profile × 845) and `DataProducts/WaterCloudProfile/ExtinctionProfile532` (profile × 501). Everything else is one value per profile or layer boundaries, stored as (N, 1) or (N, 10).
- `CloudMask/Type`: 0 Water dominant, 1 Ice/Water mix **or** Ice/HOI mix, 2 Ice dominant, 3 HOI, 4 Aerosol (probably incl. clear air: 64 % of non-NaN pixels), plus **undocumented class 5** (12,906 pixels on 02-17) `[G17]`, `[G23]`. PI says the mask is averaged to 10 s × 15 m, but it varies within 10 s windows `[G23]`.
- `WaterCloudProfile/Altitude`: 1.249 m spacing, $\Delta z \in [-374.7, +249.8]$ m relative to `Cloud Top Height` `[G13 closed]`. Profile units are m⁻¹ (attribute says m). Profiles stop at $\tau \approx 2$; they are Crosbie $\alpha^*$: first guess only, never $\mathbf{y}$.
- 1-D: `Cloud/{Cloud Top Height, Multilayer Cloud Top/Base Heights (N,10), IAB532, Sc532, LID532, WaterExtinction532}`, `Surface/{correctedReflectance, paraReflectance, surfaceDepol, surfaceMol, WindSpeed, altSurface}`, `Nav/{gps_lat, gps_lon, gps_alt, gps_time (hours)}`. Names contain spaces. Root datasets `000_Readme`, `000_revision_date`.
- Metadata conflicts awaiting the PI: appendix §A.3 (`LID532`, layer-array shape, `Sc532`/`IAB532` dependence).
- No uncertainties `[G12]`.

**Checks for task (a)** — done 2026-10-06: names, 2-D variables, time base (20 × 0.5 s per 10 s profile), `WaterCloudProfile/Altitude`. **Still open:** full header dump; the window convention of the 10 s profile relative to its 20 sub-profiles (start/mid/end) — needed for co-registration.

*Task (a) run 2026-10-14: results in `docs/reports/task_a_readers.md`; new facts in `docs/findings_for_pk.md` until folded into v1.11.*

### 1.3 CloudCube L2 (`cloudcube`)

Appendix §A.4 (readme + JPL email 2026-10-08). File headers first inspected 2026-10-14 (see `findings_for_pk.md`).

- NetCDF4, CF-1.9, L2, JPL. Variables `zm_{ka,w,g}` (dBZ), `katten_{ka,w,g}` (dB km⁻¹), `s0_{ka,w,g}` (dB), `range` (m from aircraft — slant range), `time` (s since 00 UTC, interval mid-time), `lat`, `lon`, `alt`, `roll`, `pitch`, `heading`, `incidence_angle`, `look_vector_ecef`, `landmask` (0 ocean, 1 river, 2 ice, 3 land), `nxyz`.
- **`zm` is NOT gas-corrected; `katten` is ONE-WAY** specific attenuation (JPL). $F(\mathbf{x})$ carries the full two-way gas path; display-only correction $Z_{gc}(r) = z_m + 2\int_0^r k_{\text{atten}}\,dr'$.
- **$N = 162$ independent samples per gate** (JPL); high-SNR random error ≈0.33 dB. SNR definition behind the sensitivities unknown `[G1]`, `[G24]`.
- Sensitivity at 1 km is a **global attribute per granule** (`radar_sensitivity_{ka,w,g}_dbz`). Appendix: Ka −19.26, W −32.64, G −37.58 dBZ. Read it from each file — the Ka value differs between 02-17 granules (finding 2026-10-14).
- **No uncertainty, MDS/SNR or QC fields** `[G24]`. **No Doppler for the whole campaign** `[G16 closed]`.
- Scope `[D12 amended 2026-10-09]`: W+G baseline; Ka an optional ice/snow row block above the melting layer per `[O15]`, conditional on the 02-17 Ka coverage check (task d). Readers load all three bands.

**Checks for task (a)**
1. Full header: dims, gate spacing, fill values, every attribute; is `zm` thresholded or noise-subtracted `[G24]`?
2. Time overlap with HALO and the relationship between the two time bases.
3. Range-to-altitude: `alt − range·cos(off-nadir)` from `look_vector_ecef`, compared with HALO `z` and the surface. Report; the regridding itself is task (d), not the reader.

*Task (a) run 2026-10-14: results in `docs/reports/task_a_readers.md`; new facts in `docs/findings_for_pk.md` until folded into v1.11.*

### 1.4 MetNav ICARTT (`metnav`)

`NURTURE-MetNav_G3_20260217_R0.ict`. ICARTT (FFI 1001), 1 Hz `Time_Start` (s), 24 variables incl. `Pitch_Angle` (up+), `Roll_Angle` (rt+), `True_Heading`, `Static_Air_Temp`, `Static_Pressure`. Fill −99999, LLOD flag −8888, ULOD flag −7777. Candidate 1 Hz attitude source for `[G21]`.

**Checks for task (a):** parse the header; variables, units, rate, fills; compare attitude with HALO `imu_pitch`/`imu_roll`.

*Task (a) run 2026-10-14: results in `docs/reports/task_a_readers.md`; new facts in `docs/findings_for_pk.md` until folded into v1.11.*

### 1.5 STM co-registration (task d; appendix §A.8(4), §A.11(b))

Segment 15:30–16:00 UTC (TPV transect). Time-match + CloudCube slant range → HALO 15 m altitude grid. Panels: (1) `532_total_attn_bsc`, (2) `532_dep`, (3) $Z_W$ gas-corrected only, (4) $\mathrm{DWR}_{W\text{-}G}$ with the C&S class-0 contour. Ka coverage check: fraction of ice pixels with W and G detected where Ka also detects, $Z_{\min}^{\mathrm{Ka}}(r) = Z_{1\,\mathrm{km}} + 20\log_{10}(r/1\,\mathrm{km})$, $r$ from `alt`.

*Run 2026-10-14: `docs/reports/task_d_coreg.md`.*

---

## 2. Scattering databases (`ssp_root`)

### 2.1 TAMUice2016 (`tamuice2016`) — lidar and visible/NIR ice source `[DECISION D18]`

- Local data: `Data_0.2_15.25/` only (0.2–15.25 µm, 396 wavelengths); the 16.4–99 µm set is absent and not needed.
- Nine habits × three roughness levels (`Rough000` smooth, `Rough003` moderate, `Rough050` severe).
- Folder names: `single_column`, `plate`, `hollow_column` (README says `HC`), `droxtal`, `HBR`, `SBR`, `8_columns`, `5_plates`, `10_plates`.
- `isca.dat`: 74 844 rows = 396 λ × 189 sizes, **size varies fastest**. Columns: λ (µm), $D_{\max}$ (µm), volume (µm³), projected area (µm²), $Q_{\text{ext}}$, $\omega_0$, $g$. Sizes 2–10 000 µm.
- Phase-matrix files: `P11.dat`, `P12.dat`, `P22.dat`, `P33.dat`, **`P43.dat`** (README lists P34), `P44.dat`. Each has one header line of 498 scattering angles (0–180°), then 74 844 rows in the same order as `isca.dat`. **P12…P44 are stored as ratios to P11.** $\iint P_{11}\,d\Omega = 4\pi$ — equivalent to `index.md` §4.2's $\int_0^\pi P_{11}\sin\theta\,d\theta = 2$.
- Size: about 600 MB of ASCII per P-file, about 97 GB total.
- Still unverified `[GAP G20]`: which native wavelengths bracket 0.355, 0.532, 0.935, 1.064 and 1.645 µm (wavelength interpolation is an $\mathbf{M}$ term); refractive-index source and temperature treatment.
- Known physics caveat: the Yang-lineage backscatter-peak approximation near 180° for $D \gtrsim 40$ µm `[OPEN O26]`.

**`legacy/parse_yang_db.ipynb`** — the logic is sound (it de-normalizes P12…P44 by multiplying by P11, and computes $C_{\text{ext}} = Q_{\text{ext}} \times$ area). Issues to fix when porting:
- Wavelength selection uses `np.isclose(λ, target, atol=5e-3)` with no check that exactly 189 rows match. If two native λ fall within ±5 nm, sizes are silently mixed. Assert 189 rows and report the native λ actually selected.
- It reads each full P-file into pandas per call. Cache once per (habit, roughness, λ set).
- Paths are hard-coded to an old machine.

### 2.2 ARTS SSDB v1.1 (`arts_ssdb`) — provisional radar/radiometer source `[OPEN O27]`, `[DECISION D20]`

- `SSD/TotallyRandom/` is the primary data: one netCDF4 file per size, per habit folder.
  - Ice/SingleCrystals/Pristine: IDs 2–7, 9–16, 24, 27, 31, 36.
  - Ice/Aggregates/Pristine: IDs 1, 8, 17–22, 26, 28, 30, 32, 34.
  - Ice/Aggregates/Rimed: IDs 23, 29, 33.
  - Liquid: one folder.
- `SSD_ARO/AzimuthallyRandom/`: **liquid only locally** (`LiquidSpheroidAREq2Ch_Id35`, `ChebyshevRainDrop_Id37`). There is no ice azimuthally-random data on disk.
- `StandardHabits/FullSet/`: **derived** user-ready data (17 habits incl. mixes; `.mat` and ARTS XML). Two extensions were applied:
  - pchip temperature interpolation to 190–270 K in 20 K steps;
  - for large particles, properties above the DDA limit $x = 10$ are **filled by copying the highest computed frequency** up to 886.4 GHz.
  
  ⇒ Do not use StandardHabits for G-band without checking coverage; use `SSD/` for the catalog.
- **Damaged files — do not read:** `IconHail.meta.xml` (0 B), `IconHail.xml` (exactly 20 KiB), `IconHail.xml.bin` (truncated), and `IconCloudIce (copy).mat` (exactly 50 MiB). Joe will re-download or delete them.
- Vendor interface: `DataInterfaces/Python/` (`utils.py`, `assp.py`, `sph.py`, …). Needs netCDF4 and numpy; scipy for interpolation; pyarts only for ARTS-format export (not needed); SHTns only for azimuthally random data (not needed). Use it read-only, or reimplement with netCDF4.
- Documentation in the store: `ReadMe.pdf` and `TR.pdf` (technical report).

**Catalog content for task (c)** — one row per habit (ID, name, class):
1. Size grid: number of sizes, $D_{\max}$ and $D_{\text{veq}}$ ranges.
2. $m(D)$ from per-file mass: fit $m = \alpha D_{\max}^{\beta}$ and compare with the ReadMe-tabulated $\alpha$, $\beta$.
3. $A(D)$ (projected or cross-sectional area), if present in the metadata.
4. Frequency and temperature grids actually present, **per size**.
5. **Largest $D_{\text{veq}}$ with real (non-extended) data at 94.1, 228 and 247.2 GHz.** The $x \le 10$ DDA limit means large particles lack the highest frequencies; G-band (~239 GHz) needs interpolation between 228 and 247.2 `[GAP G1]`.
6. Orientation: totally random vs azimuthally random.
7. Name matches with TAMUice2016 habits (plate, column, bullet rosette, 8-column aggregate). Flag them only: their geometric equivalence is unverified and is an `[O27]` question for a later session.

### 2.3 SIR-LUT (`sirlut`) — liquid spheres at lidar wavelengths (candidate source; no `[DECISION]` yet)

- Chemyakin et al. (2021, *Front. Remote Sens.* 2:711106) Lorenz–Mie LUT using scale invariance from a base wavelength.
- The local file `CloudLUT_100_100_0.350_700bins_1400CRI_ln2rKr_Twomey.mat` is 3.5 GB, MATLAB v7.3 (HDF5).
  - Group `CloudLUT`: `RealPartsCRI`, `ImagPartsCRI`, `RadiiGridBins`, `ScatteringAngles`, and `DataRecord` holding **object references** (1400 = 56 real × 25 imaginary).
  - Per the legacy notebook: 700 radii and 203 angles. Variables `Cext`, `Csca`, `Cbsc`, `C11`, `C12`, `C33`, `C34`.
- The vendor `read.me` describes an **aerosol** LUT: radii 0.001–100 µm, $m_r$ 1.29–1.65, $m_i$ 0–0.05, base 0.355 µm. The vendor example `sir_sca_usage.py` (Python 2) reads a **different binary file** (`AerosolLUT_...0.355_650bins_2325CRI...dat`) that is not present. The cloud-LUT ranges therefore come only from the legacy notebook — unverified.
- **Semantics to verify before porting:** the vendor code treats `Cext`, `C11`, … as **kernel coefficients**, combined with PSD-dependent quadrature weights (`coef1`/`coef2`/`coef3`, "ln2rKr Twomey") and a ×2 correction at the first bin near 355 nm. The legacy notebook instead integrates them as cross-sections with a trapezoid rule, and labels them alternately µm² and µm⁻¹.

**`legacy/load_sir_lut.ipynb`** — issues:
- The (real, imag) unpacking order of the 1400 references is assumed (`i // 25`), not verified.
- `get_liquid_properties` uses keys `C_bsc`/`C_ext`, but the dataset uses `Cbsc`/`Cext`.
- The example $r_{\text{eff}}$ sweep of 10–2000 µm exceeds plausible cloud-droplet sizes and possibly the LUT range.
- Paths are hard-coded.

Not part of tasks (a) or (c).

### 2.4 TAMUoic2019 (`tamuoic2019`) — oriented ice, NURTURE QC/$\mathbf{M}$ only `[DECISION D19]`

The folder is **empty**: neither data nor the author's read routines are present. Not current work.
