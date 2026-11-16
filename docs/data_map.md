# data_map.md — data sources for nurture-oe

**Distilled from** `appendix_A_nurture_datasets.md` (2026-11-11) and `index.md` **v1.15** (Project Knowledge), plus a directory inspection of the local data stores (2026-10-05). §1 re-synced to the v1.12 PK snapshot on 2026-10-20; §1.6 added at the v1.13 fold (2026-10-25). Re-stamped at the v1.14 fold (2026-11-11): no data fact changed; v1.14 added the lidar multiple-scattering model (`appendix_C_lidar_multiple_scattering.md`), whose test references live in `tests/lidar_ms/reference/` and are not a data store. Re-stamped at the v1.15 fold (2026-11-16): no data fact changed. Project Knowledge is authoritative; this file is a working summary. Register IDs (D/G/O) refer to `index.md` §11.

**PK snapshot.** `docs/pk_snapshot/` holds copies of `index.md` (v1.15), `appendix_A_nurture_datasets.md`, `appendix_B_stm_figures.md`, `appendix_C_lidar_multiple_scattering.md`, `handoff.md` and `team_followups.md`, edited in place at the v1.15 fold (2026-11-16). The originals are on OneDrive; never write there. Findings of 2026-10-14, 2026-10-20, 2026-10-25 and 2026-11-12 are folded. Later facts stay in `docs/findings_for_pk.md`, not here, until they are folded.

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
- **Window `[G23]`:** a profile at t averages the C&S sub-profiles in **[t − 5 s, t + 5 s)**. This is the working co-registration convention; a ±10 s sweep could not resolve the offset, and in 10.5 s steps one sub-profile lies in no window (v1.12, §A.8(5b)). Profile spacing is irregular (median 10.0 s, 9.04–23 s), so consecutive windows can overlap. `lat`/`lon` are positions at t.
- Undefined within 0.5–2.5 km below the aircraft (product-dependent; `532_ext` 2.5 km).

**Name and attribute traps**
- Variable names start with digits (`532_ext`). Access them with explicit strings.
- Attribute names contain trailing spaces (`"Description "`, `"Units "`). Several variables carry both `Units ` and `units`, sometimes disagreeing (`Angstrom_Spherical`: `"sr"` vs `"1"`). Use CF `units` and warn on disagreement.
- `_FillValue = NaN` everywhere; no `mask_low`, no fill-reason flag `[G23]`.
- `State/Pressure` is in **atm**. `State_Type` scalar: 1 balloon sonde, 2 MERRA-2, 3 LaRC O₃ sonde, 4 dropsonde, 5 HRRR; **= 2 (MERRA-2) on 02-17** (closed). 1064 products list dropsondes as state source, 532 products MERRA-2 `[G23]`.
- **`Nav_Data/imu_pitch`, `imu_roll`: do not use for incidence angle `[G21]`.** On 02-17 `imu_roll` is neither MetNav sampled at t nor a centered 10 s mean; record it as an unexplained timing/averaging offset (or a different IMU). `imu_pitch` − MetNav pitch has median +0.29°. Attitude comes from MetNav (1 Hz).
- `UserInput/DEM_altitude`: long_name GMTED2010, global attribute GLOBE (cosmetic).
- Header ↔ README conflicts: appendix §A.10.3 (16 items).

**Roles and uncertainty** (appendix §A.10.2, §A.10.4, §A.10.6; `[D21]`, `[D22]`, `[G22]`)

| Variable(s) | Res. | `_unc`? | Role |
|---|---|---|---|
| `532_total_attn_bsc` (`1064_…`) | 10 s, 15 m | none | $\ln\beta'$ rows in liquid / optically thick cloud `[D21]` |
| `532_dep` (`1064_dep`) | 10 s | screened only | $\delta_v$ rows with $\beta'$; `1064_dep` rows take their mask from the product, not `_unc` (`1064_dep_cloud_screened_unc` is finite on 10,971 pixels where its product is NaN) `[G22]` |
| `532_bsc`, `532_bsr` (+ screened) | 10 s | screened only | $\ln\beta_p$ rows where HSRL-valid (above cloud, thin ice) `[D21]`, switch `[O33]`; `bsr` = $\beta_p/\beta_m$ |
| `1064_bsc`, `1064_bsr` | 10 s | screened only | Fernald, normalized to 532 ⇒ not independent of 532 |
| `532_ext`, `532_Sa` | 60 s, 300 m running mean | yes | not rows alongside $\beta_p$ + $\beta'$; excluded `[O33]` |
| `1064_ext`, `532_extinction_from_backscatter` | | none | derived; exclude |
| `h2o_mmr_v` | 60 s, 315–705 m | yes | $\ln q_v$ rows out of cloud; in-liquid values not credible; the README's "stops 100 m above cloud top" does **not** hold on 02-17 (median lowest valid bin 495 m below `cloud_top_height_935`), so exclude in-cloud rows with **our own cloud mask** `[O32]`; in-ice credibility open |
| `Aerosol_ID` | 10 s | — | class 1 = Ice (9 classes in header) |
| `/State` | (time, z) | none | MERRA-2 $T$, $p$ for $F(\mathbf{x})$ and $k_{\text{gas}}$ |

02-17 fill fractions in cloud (§A.10.5(b)): unscreened 532 $\beta'$/$\delta$/$\beta_p$ 0.91 liquid, 0.80 ice; 1064 0.17 / 0.48; screened fields and their `_unc` 0.055 / 0.50 (filled on identical pixels, except `1064_dep_cloud_screened_unc`); `532_ext` 0.04 / 0.46; `h2o_mmr_v` 0.66 / 0.93 under this definition (0.97 of 2026-10-06 is reproduced only by other masks; §A.8(5c)). In ice, DIAL $\mathrm{RH}_w$ median 89 %, with 19.5 % > 105 % (§A.10.5(i)). Cloud pixels there were defined from ≥50 % of the twenty C&S sub-profiles (liquid = class 0, ice = classes 1–2), >2.5 km below the aircraft: 10,395 liquid, 62,567 ice.

Running means are kept on the full grid ⇒ adjacent rows correlated `[O29]`. Readers return them as-is.

**Checks for task (a)** — done 2026-10-14 (`docs/reports/task_a_readers.md`; folded into v1.11): header dump, `State_Type` and units, `imu_*` vs MetNav, NaN census, `h2o_mmr_v` truncation. **Still open:** complete diff against 01-24 (needs the 01-24 dump in the repo).

### 1.2 HALO Cloud & Surface L2 (`halo_cloudsurface`)

Appendix §A.2. 02-17 granule: **33,765 profiles at 0.5 s**, 845 bins identical to Aerosol/WV `z` (§A.2.4).

- **Only two profile datasets:** `DataProducts/CloudMask/Type` (profile × 845) and `DataProducts/WaterCloudProfile/ExtinctionProfile532` (profile × 501). Everything else is one value per profile or layer boundaries, stored as (N, 1) or (N, 10).
- `CloudMask/Type`: 0 Water dominant, 1 Ice/Water mix **or** Ice/HOI mix, 2 Ice dominant, 3 HOI, 4 Aerosol (probably incl. clear air: 64 % of non-NaN pixels), plus **undocumented class 5** (12,906 pixels on 02-17) `[G17]`, `[G23]`. PI says the mask is averaged to 10 s × 15 m, but it varies within 10 s windows `[G23]`.
- `WaterCloudProfile/Altitude`: 1.249 m spacing, $\Delta z \in [-374.7, +249.8]$ m relative to `Cloud Top Height` `[G13 closed]`. Profile units are m⁻¹ (attribute says m). Profiles stop at $\tau \approx 2$; they are Crosbie $\alpha^*$: first guess only, never $\mathbf{y}$.
- 1-D: `Cloud/{Cloud Top Height, Multilayer Cloud Top/Base Heights (N,10), IAB532, Sc532, LID532, WaterExtinction532}`, `Surface/{correctedReflectance, paraReflectance, surfaceDepol, surfaceMol, WindSpeed, altSurface}`, `Nav/{gps_lat, gps_lon, gps_alt, gps_time (hours)}`. Names contain spaces. Root datasets `000_Readme`, `000_revision_date`.
- Metadata conflicts awaiting the PI: appendix §A.3 (`LID532`, layer-array shape, `Sc532`/`IAB532` dependence).
- No uncertainties `[G12]`.

**Checks for task (a)** — done 2026-10-14 (folded into v1.11). Window convention: see §1.1 (working convention `[G23]`).

### 1.3 CloudCube L2 (`cloudcube`)

Appendix §A.4 (readme + JPL email 2026-10-08; headers of the three 02-17 granules inspected 2026-10-14, folded into v1.11).

- NetCDF4, CF-1.9, L2, JPL. Variables `zm_{ka,w,g}` (dBZ), `katten_{ka,w,g}` (dB km⁻¹), `s0_{ka,w,g}` (dB), `range` (m from aircraft — slant range), `time` (s since 00 UTC, interval mid-time), `lat`, `lon`, `alt`, `roll`, `pitch`, `heading`, `incidence_angle`, `look_vector_ecef`, `landmask` (0 ocean, 1 river, 2 ice, 3 land), `nxyz`.
- Dims `(time, range)` and `look_vector_ecef (time, nxyz)` (the readme says `(range, time)`, `(3, time)`). 751 gates 19.986 m apart from −999.3 m (r = 0 at index 50; r < 0 unused); rays every 0.10044 s. `_FillValue = NaN` for `zm`/`katten`/`s0`; nav has no `_FillValue` (9.97e36 → NaN); `landmask` −2147483647. Nav/attitude are MetNav interpolated (lag 0 s).
- **`alt` is above the WGS84 ellipsoid** (20–44 m below HALO MSL on 02-17) and equals MetNav MSL + $N_{\text{EGM96}}$ to ±8 cm. **Never use it for altitude `[DECISION D23]`:** gate altitude $z = \text{alt}_{\text{HALO}} - r\cos\theta$, with HALO `gps_alt` (MSL) interpolated to ray times and $\cos\theta = -\mathbf{l}_{\text{ECEF}}\cdot\hat{u}$ from `look_vector_ecef`.
- **`zm` is NOT gas-corrected; `katten` is ONE-WAY** specific attenuation (JPL). $F(\mathbf{x})$ carries the full two-way gas path; display-only correction $Z_{gc}(r) = z_m + 2\int_0^r k_{\text{atten}}\,dr'$.
- **$N = 162$ independent samples per gate** (JPL); high-SNR random error ≈0.33 dB. SNR definition behind the sensitivities unknown `[G1]`, `[G24]`.
- Sensitivity at 1 km is a **global attribute per granule** (`radar_sensitivity_{ka,w,g}_dbz`). W −32.64 and G −37.58 dBZ in all three 02-17 granules; **Ka per granule: −16.32 dBZ for 15:13–18:59, −19.26 for 19:00–19:26.** Read it from each file. The granule's own Ka value is the **admission threshold** of the optional Ka block, $Z_{\mathrm{Ka}} > Z_{1\,\mathrm{km}} + 20\log_{10}(r/1\,\mathrm{km})$ `[O15]`. At 19:00 the Ka range response also changed (≈1 km skirt around strong echoes; σ0 Ka − W −1.6 to −3.3 dB). Mask 19:00-granule Ka rows near strong echoes; calibration per granule (§A.8(5e)).
- **The documented 1 km sensitivity is not the processor's detection threshold:** delivered `zm` reaches 5–11 dB below $Z_{\min}(r)$ in every band (0.1th percentile; Ka −7.7 to −9.0, W −9.8 to −10.9, G −6.3 to −7.0 dB) `[G24]`. Don't use the attribute as the noise floor in the $\mathrm{SNR}^{-1}$ term `[O31]`.
- **No uncertainty, MDS/SNR or QC fields** `[G24]`. **No Doppler for the whole campaign** `[G16 closed]`.
- Scope `[D12 amended 2026-10-09]`: W+G baseline; Ka an optional ice/snow row block above the melting layer per `[O15]` (condition 3 met 2026-10-14: 4.0 % coverage on 15:30–16:00, near-surface snow only; condition 1, Lebsock's concurrence, pending). Readers load all three bands.
- `look_vector_ecef` is computed from attitude (body-fixed nadir beam assumed). Surface ranges in turns support body-fixed in roll; pitch, lever arm and boresights unknown `[G27]`. Off-nadir ≈2.4° on the 15:30–16:00 straight leg.

**Checks for task (a)** — done 2026-10-14 (folded into v1.11). Still open with JPL: SNR definition and processor threshold, noise-subtraction state of `zm`, meaning of NaN `[G24]`; range resolution, near-field limit `[G1]`.

### 1.4 MetNav ICARTT (`metnav`)

`NURTURE-MetNav_G3_20260217_R0.ict`. ICARTT (FFI 1001), 1 Hz `Time_Start` (s), 24 variables incl. `Pitch_Angle` (up+), `Roll_Angle` (rt+), `True_Heading`, `Static_Air_Temp`, `Static_Pressure`. Fill −99999, LLOD flag −8888, ULOD flag −7777. Candidate 1 Hz attitude source for `[G21]`.

**The attitude source** for incidence angle (`[G21]` interim rule, v1.11). `GPS_Altitude_MSL` equals HALO `gps_alt` in the median (+0.1 m), but drifts −6 to +4 m from it over long stretches (v1.12). There is no ellipsoidal altitude field. 1 Hz with no gaps on 02-17; no LLOD/ULOD values.

### 1.5 STM co-registration (task d; appendix §A.8(4), §A.11(b))

Done 2026-10-14 (`docs/reports/task_d_coreg.md`; folded into v1.11, appendix §A.8(4), §A.11(b)). Segment 15:30–16:00 UTC (TPV transect). Each HALO profile at t takes the C&S sub-profiles and CloudCube rays in [t − 5 s, t + 5 s) `[G23]`; gate altitudes on the HALO MSL datum `[D23]`; nearest gate per 15 m bin. Ka coverage: $Z_{\min}^{\mathrm{Ka}}(r) = Z_{1\,\mathrm{km}} + 20\log_{10}(r/1\,\mathrm{km})$ with $r$ the slant range and $Z_{1\,\mathrm{km}}$ from the granule.

- **Co-registered radar rows are censored means** over detecting rays: biased high where the detected fraction is < 1 `[O21]` (sized in follow-up §A.8(5d)).
- **About one gate in three is duplicated** into two 15 m rows (perfectly correlated errors); decimate to native gates before building $\mathbf{S}_\epsilon$ `[O29]`.
- Window means carry the same representativeness error as the lidar `[O36]`.

*Follow-ups run 2026-10-20 (`docs/reports/followup_*.md`, folded into v1.12): geoid check, window sweep, 0.5 s × 15 m file `coreg_20260217_1530-1600_0p5s.nc` (`scripts/build_coreg_0p5s.py`; superseded 2026-10-25 by `coreg_cloud_20260217_1530-1600_0p5s.nc`, native 0.5 s lidar), DIAL RH_i, censored mean, Ka step at 19:00, mounting.*

### 1.6 HALO 0.5 s cloud file (`nurture_example`/`Cloud_Plots/`, key `halo_cloud`)

*Inspected 2026-10-25; folded into v1.13 (Appendix A §A.10.7, §A.8(7)).* `20260217_cloud_data.mat`, 2.73 GB, MATLAB v7.3 (HDF5, h5py), delivered 2026-10-23; interim, a CF-compliant version is announced. The `.fig`/`.png` files in the folder are the PI's quick-looks (MATLAB v5 figures).

- One struct `data`; groups `nav`, `cloudphase`, `cloud`, `column` (+`OD`), `cprofile`, `sprofile`, `surface`, `reference`, `QA`. Datasets are (profile, …) float64, gzip.
- 33,765 profiles at 0.5 s, nav identical to C&S. `cloudphase` (845 bins, `alt` = HALO `z`): `bsc`, `dep`, `bsc1064`, `dep1064`, `cp`, `S0`, `cT2p`, `X1064m_syn`. `cprofile` 501 bins cloud-top-relative (= C&S `WaterCloudProfile/Altitude`); `sprofile` 301 bins ±187 m around the surface.
- **No attributes at all** (no units, descriptions, fills). Reader `io/halo_cloud.py` attaches only `*_inferred` attributes backed by checks (`docs/reports/cloud_0p5s.md`).
- Traps: `bsc` is m⁻¹ sr⁻¹ (Aerosol/WV uses km⁻¹ sr⁻¹); `dep` is particulate, unscreened (|δ| > 1 in clear air); `bsc1064` is not on the Aerosol/WV 1064 scale; xarray `Dataset.roll` shadows the `surface` variable `roll` (use `ds["roll"]`).

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

*Task (c) run 2026-10-20: `docs/reports/arts_catalog.md` + CSVs (folded into v1.12, index §5.2): Ka and W cover all sizes; 228/247.2 GHz stop at Dveq ≈ 4.0–4.5 mm for 8 habits; 34 frequencies on disk.*

**Catalog content for task (c)** — one row per habit (ID, name, class):
1. Size grid: number of sizes, $D_{\max}$ and $D_{\text{veq}}$ ranges.
2. $m(D)$ from per-file mass: fit $m = \alpha D_{\max}^{\beta}$ and compare with the ReadMe-tabulated $\alpha$, $\beta$.
3. $A(D)$ (projected or cross-sectional area), if present in the metadata.
4. Frequency and temperature grids actually present, **per size**.
5. **Largest $D_{\text{veq}}$ with real (non-extended) data at the Ka-band grid frequencies, 94.1, 228 and 247.2 GHz** (Ka added 2026-10-20; CloudCube's exact Ka frequency is unknown `[G1]`). The $x \le 10$ DDA limit means large particles lack the highest frequencies; G-band (~239 GHz) needs interpolation between 228 and 247.2 `[GAP G1]`.
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
