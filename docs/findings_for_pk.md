# findings_for_pk.md — findings to fold into Project Knowledge

Append-only. Template and rules: see `CLAUDE.md`. Joe and Fable fold these into `index.md` / `appendix_A_nurture_datasets.md` and mark each entry `FOLDED <date>`.

### 2026-10-14 — CloudCube Ka sensitivity attribute differs between 02-17 granules
**FOLDED 2026-10-14**
Type: conflict-with-PK
Refs: G1, G24, O15, O31  Source: global attr `radar_sensitivity_ka_dbz` in all three 02-17 CloudCube granules
Finding: The 15:13–15:20 and 15:21–18:59 granules declare Ka **−16.32 dBZ @ 1 km**; the 19:00–19:26 granule declares −19.26, the value in appendix A §A.4 and O15. W (−32.64) and G (−37.58) agree in all three. The delivered Ka floor moves with the attribute: the 0.1th percentile of $z_m - Z_{\min}(r)$ is −9.0 dB for both granules (6–9 km range), so the ~3 dB change looks real, not a typo.
Impact: $Z_{\min}^{Ka}(r)$ in the O15 coverage check (4.0 % vs 10.9 % on 15:30–16:00); the Ka random-error term of the interim radar $\mathbf{S}_\epsilon$ (O31); any Ka presence mask. The sensitivity must be read per granule. Ask JPL at the STM what changed at 19:00.

### 2026-10-14 — Delivered zm extends 5–11 dB below the stated Z_min(r) in every band
**FOLDED 2026-10-14**
Type: GAP-update
Refs: G24, G1, O21, O31  Source: `zm_{ka,w,g}` vs `radar_sensitivity_*_dbz`, granules 15:21–18:59 (15:30–16:00) and 19:00–19:26
Finding: The 0.1th percentile of $z_m - Z_{\min}(r)$ (from each granule's own attribute) is Ka −7.7 to −9.0, W −9.8 to −10.9 and G −6.3 to −7.0 dB at 6–11.5 km range. The offset is about constant with range, so the $20\log_{10} r$ form holds. So the documented "sensitivity @ 1 km" is not the detection threshold of the delivered data; the SNR definition behind it is still unknown.
Impact: The SNR term $(1+\mathrm{SNR}^{-1})$ of the interim radar $\sigma_{\text{rand}}$ (§A.4(iv)) cannot use the attribute as the noise floor without the SNR definition. It also affects the O21 non-detection threshold and the "threshold slightly above $Z_{\min}$" rule (index §2.3), and which definition the O15 test uses ("finite" 12.8 % vs "> Z_min" 4.0 %).

### 2026-10-14 — CloudCube `alt` is ellipsoidal; HALO `z`/`gps_alt` are MSL (20–44 m apart)
**FOLDED 2026-10-14**
Type: fact
Refs: G1, O36; appendix A §A.11(b)  Source: CloudCube `alt:long_name` ("above WGS84 ellipsoid"); HALO `Nav_Data/gps_alt`; MetNav `GPS_Altitude_MSL`
Finding: HALO `gps_alt` equals MetNav MSL altitude (median +0.1 m), so HALO `z` is MSL. HALO − CloudCube altitude is a median 38.1 m (1–99 %: 18.7–44.0 m) in the 15:21–18:59 granule and 20 m in the others: the geoid undulation along track. Gates referenced to HALO MSL altitude put the W surface peak a median −11 m (15:30–16:00) and +7 m (19:00–19:25) from C&S `altSurface`; with CloudCube `alt` the residuals are −40 m and −13 m.
Impact: Radar–lidar vertical registration (about two gates if ignored); DWR alignment with lidar layers and the liquid mask; the S2 retrieval-grid choice. Co-registration uses HALO MSL altitude (`coreg.cloudcube_to_halo_grid`, `aircraft_alt="halo"`). The antenna lever arm is still unknown.

### 2026-10-14 — CloudCube L2 header facts (first header inspection)
**FOLDED 2026-10-14**
Type: GAP-update
Refs: G1, G24  Source: `docs/reports/headers/NURTURE-CLOUDCUBE_*.txt`
Finding: Three 02-17 granules: 15:13:00–15:20:10, 15:21:58–18:59:37 (4.7 GB) and 19:00:00–19:26:25. Dims are `(time, range)` and `look_vector_ecef (time, nxyz)`, not `(range, time)`/`(3, time)` as in the readme. There are 751 gates 19.986 m apart from −999.3 m (r = 0 at index 50), and rays every 0.10044 s at the interval mid-time. Fills: `zm`/`katten`/`s0` use `_FillValue = NaN`; nav variables have no `_FillValue` (netCDF default 9.97e36); `landmask` uses −2147483647. `katten_ka`/`katten_w` are finite at 22–25 % of r < 0 gates. Arrays are contiguous and uncompressed. The off-nadir angle from `look_vector_ecef` matches `incidence_angle` to 0.1°. CloudCube nav/attitude are MetNav interpolated (lag 0 s).
Impact: Gate spacing for $\mathbf{K}$/$F$ range weighting (the range *resolution* is still unknown, G1); co-registration; reader design.

### 2026-10-14 — HALO 10 s profile averages the C&S sub-profiles in [t − 5 s, t + 5 s)
**FOLDED 2026-10-14**
Type: GAP-update
Refs: G23, O36  Source: HALO `time` vs C&S `Nav/gps_time`; `IAB532` vs column `532_total_attn_bsc`
Finding: HALO profile spacing is irregular (median 10.0 s, range 9.04–23 s), and `time:long_name` reads "instant time in UTC". Over 273 heterogeneous profiles, the centered window [t − 5, t + 5) gives the best agreement between windowed C&S `IAB532` and column $\beta'_{532}$ (correlation 0.52, vs 0.41 for [t, t + 10)). It holds 20 sub-profiles in >95 % of profiles. Alignment is resolved to about ±1 s. HALO `lat`/`lon` are positions at t (they match C&S and MetNav within 1 s), not window means. With this rule the §A.10.5 counts reproduce within 0.4 % (10,439 liquid / 62,700 ice vs 10,395 / 62,567); the 2026-10-06 window convention was not recorded.
Impact: The window rule for all co-registration and for $\mathbf{S}_{\text{repr}}$ from mask heterogeneity (O36). PI confirmation wanted.

### 2026-10-14 — HALO `imu_roll` is neither a sample at t nor a centered 10 s mean
**FOLDED 2026-10-14**
Type: GAP-update
Refs: G21, O28  Source: HALO `Nav_Data/imu_roll`, `imu_pitch` vs MetNav `Roll_Angle`, `Pitch_Angle` (1 Hz)
Finding: Over the 50 turning profiles (|roll| > 5°), RMS(HALO − MetNav roll) is 6.2° against MetNav at t and 5.3° against the centered 10 s mean. The best fit is a 10 s mean centered 6 s *after* t (2.3°). The clocks agree within 1 s (positions), so this is attitude-specific. `imu_pitch` − MetNav pitch has a median of +0.29°.
Impact: The HALO incidence angle for D19 flags and $\mathbf{M}$ sizing: do not take $\theta_{\text{inc}}$ from HALO `imu_*` until the PI confirms their definition. MetNav 1 Hz (window-averaged) is the defensible source.

### 2026-10-14 — `h2o_mmr_v` is not truncated 100 m above cloud top in most profiles
**FOLDED 2026-10-14**
Type: GAP-update
Refs: O32, D7  Source: `DataProducts/h2o_mmr_v`, `cloud_top_height_935`
Finding: Of the 1,329 profiles with `cloud_top_height_935`, only 25 % have their lowest valid DIAL bin ≥ 85 m above it. The median lowest valid bin is 495 m *below* `cloud_top_height_935` (5–95 %: −3,598 to +696 m). This is consistent with §A.10.5's in-liquid values.
Impact: $\ln q_v$ rows: the in-cloud exclusion (§A.10.6) must use our own cloud mask, not the product's documented truncation.

### 2026-10-14 — `1064_dep_cloud_screened_unc` filled where its product is NaN
**FOLDED 2026-10-14**
Type: conflict-with-PK
Refs: G22, G23  Source: `DataUncertainty/1064_dep_cloud_screened_unc` vs `DataProducts/1064_dep_cloud_screened`
Finding: Five of the six `_cloud_screened_unc` fields are finite on exactly the pixels of their screened product, as appendix A §A.10.4 states for all of them. `1064_dep_cloud_screened_unc` is the exception: it is finite on 10,971 pixels where `1064_dep_cloud_screened` is NaN (899,377 both finite).
Impact: Minor. Any 1064 $\delta$ row must take its mask from the product, not from `_unc`. Add to the HALO metadata list (§A.10.3).

### 2026-10-14 — State and variable inventory of the 02-17 Aerosol/WV granule
**FOLDED 2026-10-14**
Type: fact
Refs: G23; appendix A §A.8(1), §A.10.2  Source: `State/*`, full header dump
Finding: `State_Type` = 2 (MERRA-2). Units: Temperature K, Pressure **atm** (0.168–1.051), Number_Density m⁻³, Relative_Humidity %, U/V m/s. Only Nav/State/UserInput and root `time`/`z` lack dimension scales; `DataProducts`/`DataUncertainty` have `(time, z)` attached. Every variable named in §A.10.2 is present. 23 more variables are not named there (`*_bsc_Sa` (+`_unc`), `532_AOT_from_bsc*`, several `_unc`, nav, `DataProducts/Altitude` = "legacy … same as z"). `Units `/`units` differ semantically only for `Angstrom_Spherical`. A complete diff against 01-24 needs the 01-24 dump, which is not in the repo.
Impact: Closes the §A.8(1) `State_Type` item; $T$, $p$ (convert atm) for $F(\mathbf{x})$ and $k_{\text{gas}}$.

### 2026-10-14 — Ka coverage on the 02-17 15:30–16:00 TPV transect
**FOLDED 2026-10-14**
Type: fact
Refs: O15, D12, O21  Source: `docs/reports/task_d_coreg.md`; `ka_coverage_*_20260217_1530-1600.csv`
Finding: 18,944 ice pixels (C&S classes 1–2 in ≥ 50 % of sub-profiles, all T < 0 °C) have W and G detected. Ka also detects in **4.0 %** with the granule's −16.32 dBZ, 10.9 % with PK's −19.26, and 12.8 % with any finite `zm_ka`. All Ka detections are in the lowest about 2 km (10–12 km slant range: 14 / 36 / 42 %); coverage is 0 % at 6–8 km range. With 500 m clearance above the surface: 2.2 / 8.6 %.
Impact: O15 condition 3. Ka adds rows only in near-surface snow on this transect, not where the density/riming signal is wanted aloft; elsewhere Ka contributes only via the O21 non-detection penalty. A second segment (19:00–19:25, Ka −19.26) would test generality.

### 2026-10-14 — Gas path in DWR_W−G near the surface
**FOLDED 2026-10-14**
Type: fact
Refs: G24, O32; index §4.1  Source: `katten_w`, `katten_g` (15:30–16:00), display correction $2\int k\,dr$
Finding: The median two-way gaseous attenuation at 1 km altitude (r ≈ 11 km) is 1.3 dB at W and 5.9 dB at G, so about 4.7 dB of the near-surface DWR_W−G is gas alone. `katten_g` reaches 1.9 dB km⁻¹ (99th percentile).
Impact: STM slide framing; the $k_{\text{gas}}$ consistency check; DWR in $\mathbf{y}$ must be modeled with the full gas path in $F(\mathbf{x})$, not corrected beforehand.

### 2026-10-20 — Geoid check confirms D23; CloudCube `alt` = MetNav MSL + N(EGM96)
**FOLDED 2026-10-20**
Type: fact
Refs: D23, G27  Source: `docs/reports/followup_a_geoid.md`; CloudCube `alt`, MetNav `GPS_Altitude_MSL`, PROJ EGM2008/EGM96 grids
Finding: MetNav MSL − CloudCube `alt` equals −N(EGM96) to within ±0.08 m (1–99 %) in all three 02-17 granules. With EGM2008 the residual is −0.9 to +0.7 m. −N is 14–45 m along track. MetNav has no ellipsoidal altitude field. So CloudCube `alt` is consistent with MetNav MSL converted to the ellipsoid with EGM96; the conversion is inferred, not documented by JPL.
Impact: Closes the geoid-confirmation item of §A.8(5a). D23's premise (20–44 m difference = geoid undulation) holds.

### 2026-10-20 — HALO `gps_alt` departs from MetNav MSL by up to ~6 m for long stretches
**FOLDED 2026-10-20**
Type: GAP-update
Refs: D23, G27, G21  Source: HALO `Nav_Data/gps_alt` vs MetNav `GPS_Altitude_MSL` (followup_a)
Finding: HALO − MetNav altitude has median 0.0 m over the flight but 5–95 % −6.1 to +3.8 m (std 3.2 m). It is a slow drift, not noise or lag: profile-to-profile RMS 1.4 m; RMS smallest at lag 0; unchanged in level flight. Examples: −5.2 m median over 16:15–17:10, +1.2 m over 19:00–19:25. It correlates with N (r = 0.64), which hints at a different geoid realization in HALO's GPS (unconfirmed).
Impact: Radar–lidar registration under D23 is unaffected (both referenced to HALO `gps_alt`), but "HALO MSL" and MetNav MSL are not the same absolute datum to better than ≈0.4 bin. This matters if MetNav, a DEM or MERRA-2 geopotential heights are combined with HALO `z`. Ask J. Collins which GPS/geoid HALO uses and whether `z` hangs from `gps_alt`.

### 2026-10-20 — Window offset sweep does not resolve the 10 s window offset
**FOLDED 2026-10-20**
Type: GAP-update
Refs: G23, O29, O36  Source: `docs/reports/followup_b_window.md`; C&S `IAB532` vs column HALO β′532, 273 heterogeneous profiles
Finding: Correlation vs window-centre offset (−10…+10 s, 0.5 s steps) is broad. Pearson peaks at −2.5 s (fixed 155-profile set, 0.592 vs 0.580 at 0 s); Spearman peaks at −0.5 s (0.490 vs 0.486), half-drop width −4 to +2 s. The bootstrap 95 % interval of the peak covers −10 to +9 s (Spearman: 58 % of replicates within ±1 s). Data are consistent with [t − 5, t + 5) but the offset is not resolved to a few seconds; the 2026-10-14 "±1 s" was unsupported. Bookkeeping: 79 of 1,689 windows lack 20 sub-profiles (59 hold 19); 236 consecutive pairs are < 10 s apart, 183 actually share sub-profiles, 190 sub-profiles sit in two windows, 341 of 33,765 in none.
Impact: G23 stays as written in v1.11 (working convention). The time correlation from shared sub-profiles is weak (≤ 1 of 20 per pair). A 0.5 s × 15 m co-registered file now exists for 15:30–16:00 (display and diagnostics).

### 2026-10-20 — `h2o_mmr_v` ice fill fraction: 0.932 is PK's definition; 0.97 reproduced only by other masks
**FOLDED 2026-10-20**
Type: fact
Refs: O32, G2  Source: `docs/reports/followup_c_dial.md`
Finding: With C&S ice (classes 1–2 in ≥ 50 % of sub-profiles, > 2.5 km below aircraft) the fill is 0.928–0.943 for every window and variant tested. 0.97 is reproduced by ice restricted to windows where the lidar sees the surface in ≥ 1 sub-profile (0.969, n = 51,004; complement 0.774), or by `Aerosol_ID` = 1 as the ice mask (0.971, n = 45,400). The 2026-10-06 definition is not recorded, so which one was used is undetermined. Fill varies by hour (0.69 at 17 UTC, 0.98 at 15–16 UTC).
Impact: §A.10.5(b) ice entry: 0.93 under the stated definition. DIAL gaps in ice sit in columns the lidar does not penetrate.

### 2026-10-20 — In-ice DIAL RH_i: median 112 %, deep low-level snow 143 %; no verdict
**FOLDED 2026-10-20**
Type: GAP-update
Refs: O32, D7, G2  Source: `h2o_mmr_v`, `/State` T, p, `Relative_Humidity` (followup_c)
Finding: In 58,463 C&S ice pixels, RH_i (Murphy & Koop; reproduces `RHi_v` to 0.2 %) has 5/25/50/75/95 % = 83/101/112/125/166 %, with 55 % > 110 %. By depth below the layer top, medians are 105–112 % in the top 1 km and 143 % at 2–4 km (82 % > 110 %); the > 150 % mode sits mostly in 15:30–16:00, 0.5–3.4 km above ground at 247–262 K. MERRA-2 RH in the same pixels is ≤ 100 % (median 94 %). DIAL − MERRA-2 is +25 points in ice, −6 in clear air and +9 within 500 m above ice. Random uncertainty is small (median 2 %); a 1 K T error moves RH_i by ≈10 %.
Impact: The O32 in-ice sub-question remains open. Upper-ice values are near saturation; persistent 140–165 % in precipitating snow is hard to reconcile with ice physics without a T (≈3–4 K) or DIAL systematic. Dropsonde T for 02-17 (the 1064 products cite dropsondes) would separate them. Affects D7 Architecture A ice rows and B's in-cloud RH.

### 2026-10-20 — MERRA-2 `/State/Relative_Humidity` saturation reference undocumented; behaves as over ice
**FOLDED 2026-10-20**
Type: new-GAP-candidate
Refs: O32, G23  Source: HALO `/State/Relative_Humidity` (`Description` "Relative Humidity in percent")
Finding: Its maximum is ≈100 % at 230–270 K, where ice saturation expressed over water would be 69–92 %. That fits an ice (or ice/liquid blend) reference, not pure water; inferred only.
Impact: Any comparison of DIAL or retrieved q_v with MERRA-2 RH, and B's in-cloud RH fallback, needs the reference stated (PI or MERRA-2 documentation).

### 2026-10-20 — Censored-mean bias of co-registered radar rows sized
**FOLDED 2026-10-20**
Type: GAP-update
Refs: O21, O29, O36, O31  Source: `docs/reports/followup_d_censored.md` (15:30–16:00, 10 s × 15 m)
Finding: `*_all` in the cache is the detecting-ray mean without the f ≥ 0.5 threshold, not an all-ray mean. Partial pixels are 14 % (W), 16 % (G) and 25 % (Ka) of valid pixels (f < 0.9: 9/10/17 %). On real partial pixels the bias (always high) is bounded by [non-detects at the floor, non-detects at 0]: 1.7–2.6 dB at f = 0.5–0.6, 0.9–1.2 at 0.7–0.8, 0.5–0.7 at 0.8–0.9 dB, the same in every band. Synthetic censoring of fully detected echoes gives about half the lower bound (0.8–0.9 dB at f ≈ 0.55). For DWR_W−G the bias has either sign (W-censored high, G-censored low); per-pixel interval medians −0.2 to +0.9 dB at f_min ≈ 0.5.
Impact: The f ≥ 0.5 rule admits rows biased by up to ≈2.5 dB, far above the ≈0.03 dB random term. Rows with f ≥ 0.9 are ≲ 0.2 dB (95th percentile). A one-sided penalty cannot absorb the DWR bias. Input to S2's choice between restricting rows and a censored likelihood (O21).

### 2026-10-20 — Ka range response changed at 19:00; Rayleigh-ice DWR test not possible; σ0 Ka−W drops 1.6–3.3 dB
**FOLDED 2026-10-20**
Type: GAP-update
Refs: G24, O15, G1, G5, O30  Source: `docs/reports/followup_e_ka_step.md`; `zm_*`, `s0_*`, `katten_*` across 18:59:37 / 19:00:00
Finding: The Ka surface echo's contiguous extent below the surface peak (no atmosphere there, so a measure of range response) jumps from 360 m to 1,200–1,300 m at 19:00 (W and G unchanged within 40 m), with ≈1 km above the peak. Later in the granule (19:20–19:26) it reaches 1.6 km, with Ka "detecting" 0.5 M gates > 1.5 km above the surface where W detects 27 k. Ka detects no atmospheric echo in 18:42–18:59:37, so the Rayleigh-ice DWR_Ka−W test has no data before the boundary. Gas-corrected σ0_Ka − σ0_W steps −3.3 dB (±1 min; 90 % −3.7 to −2.7), −2.4 (±2), −2.0 (±3) and −1.6 (±6); the G − W control moves +0.5 to +1.5 dB. All columns contain C&S liquid. Surface-peak `zm_ka` falls 3.3–3.7 dB, but a wider range response lowers the peak by itself.
Impact: The 19:00 change is a Ka configuration or processing change, not only a gain change, and the calibration-vs-sensitivity dichotomy of §A.8(5e) does not apply as posed. The Ka calibration term must be per granule. Ka rows in the 19:00 granule within ≈1–1.5 km of strong echoes are contaminated by the skirt, so the 21 % Ka fill and the −19.26 coverage bracket are likely inflated; the 15:21 granule (360 m, like W) is unaffected. Ask JPL what changed (waveform, pulse compression, gain, processing) and how `s0` is computed, including gas correction.

### 2026-10-20 — `look_vector_ecef` is computed from attitude; surface range favours a body-fixed antenna
**FOLDED 2026-10-20**
Type: GAP-update
Refs: G27, D23  Source: `docs/reports/followup_f_mounting.md`; `look_vector_ecef`, `heading`/`pitch`/`roll`, W surface range, HALO `DEM_altitude`
Finding: In body axes, `look_vector_ecef` is exactly the body z-axis (tilts 0.000°, roll up to 37.8°). θ = arccos(cos pitch · cos roll) to 0.0045° RMS in turns, so the L2 geolocation assumes a nadir, body-fixed beam with zero offset, and the commissioned θ-vs-attitude test is circular. Independent test, 427 turning samples: θ_eff = arccos(h_AGL/r_sfc) tracks θ_b (correlation 0.87, median −0.9°; r_sfc/h 1.030 vs 1/cosθ_b 1.035); a stabilized beam would give θ_eff ≈ 0.
Impact: Body-fixed is favoured. Boresight offsets of a few tenths of a degree are not resolvable this way (DEM and terrain error), and any real offset is absent from the delivered look vector: ≈0.9 m in gate altitude and ≈21 m along track per 0.1° at 12 km. G27's lever arm, boresight and the HALO `tilt_angle` stay open.

### 2026-10-20 — ARTS SSDB catalog: Ka/W cover every size; G band truncated for 8 habits; ReadMe errata
**FOLDED 2026-10-20**
Type: fact
Refs: O27, D20, G1, O23  Source: `docs/reports/arts_catalog.md`, `arts_catalog.csv`, `arts_catalog_sizes.csv` (`SSD/TotallyRandom`, 35 habits, 1,649 size files)
Finding: Ka grid frequencies 31.3/31.5/35.6 GHz and 94.1 GHz have computed data for every size of every habit. Both 228 and 247.2 GHz (needed to interpolate CloudCube G) stop below the largest size for 8 habits: ICON/GEM graupel, hail and snow (Ids 29, 30, 32, 33, 34) at Dveq 3.97 mm (247.2) / 4.46 mm (228) of 5.00; 8-Column Aggregate at 4.36 of 5.00; Large Block Aggregate at 4.32 of 4.61; Spherical Graupel at 4.53 of 5.29. 34 frequencies on disk (the ReadMe's 10.65 GHz absent); all temperatures exist wherever a frequency exists; all ice is totally random (SSD_ARO holds liquid only). The DataSummary α, β describe the large-size end: all-size fits differ (e.g. Evans Snow β 2.07 vs 2.41, Sector Snowflake 1.79 vs 1.44). ReadMe vs disk: Liquid Sphere α 480 (ReadMe; the ice value) vs 523.6 on disk; Spherical Graupel 13/2.69 vs 16.05/2.716; Sector Snowflake Dmax/Dveq max 12.0/1.415 mm vs 10.24/1.31 mm; Short Column 34 vs 35 sizes.
Impact: S1 habit/database choice (O27) and S3 radar operator tables. Ka needs no special coverage handling. G-band rows for large graupel/hail/snow need a PSD truncation or extrapolation M term above Dveq ≈ 4 mm. Use per-size masses rather than the tabulated power law for small particles.

### 2026-10-20 — In 10.5 s HALO steps one C&S sub-profile lies in no [t − 5, t + 5) window
**FOLDED 2026-10-20**
Type: GAP-update
Refs: G23, O36  Source: HALO `time` vs C&S `Nav/gps_time` (whole flight); `scripts/build_coreg_0p5s.py`
Finding: Of 245 sub-profiles that fall in no window during ordinary 10.46–10.5 s HALO steps (another 96 sit in gaps > 11 s), 226 lie between two full 20-sub-profile windows and only 15 border a 19-count window. A 10.5 s step spans 21 sub-profiles at the 0.50138 s C&S cadence. Under the working rule, one sub-profile per such step enters no HALO average. Whether HALO skips it or uses a slightly wider window is unresolved (the offset sweep cannot resolve it).
Impact: 10 s y rows are unaffected (those sub-profiles and their radar rays simply carry no lidar row). The 0.5 s display file now holds the nearest profile within 0.5 s of a window edge, flagged `halo_hold` = 2. Ask J. Collins how sub-profiles are assigned to 10 s profiles.

### 2026-10-20 — Correction: in-ice DIAL RH_i must be judged against water saturation too
**FOLDED 2026-10-20**
Type: GAP-update
Refs: O32, D7, G2  Source: `followup_c_dial.json` (`RHw_DIAL_all_ice`, `by_depth_below_layer_top_m`), corrects the 2026-10-20 entry "In-ice DIAL RH_i: median 112 %…"
Finding: The same in-ice DIAL vapour pressure gives RH_w 5/25/50/75/95 % = 63/79/89/100/144 %. Water saturation corresponds to RH_i 113–146 % (median 125 %) at these temperatures, so the RH_i median of 112 % lies between ice and water saturation, as expected in mixed-phase cloud. The implausible part is narrower: 19.5 % of ice pixels exceed RH_w 105 %, concentrated 2–4 km below the layer top (median RH_w 127 %, 62 % > 110 %; low-level snow of 15:30–16:00).
Impact: Replaces "55 % of ice pixels > 110 % RH_i is hard to reconcile". The upper ≈2 km of ice is physically plausible; the deep low-level snow is not without a ≈3 K T bias or a DIAL in-cloud systematic. O32 stays open; dropsonde T/q_v is the discriminating test.

### 2026-10-25 — HALO 0.5 s cloud file delivered for 02-17: structure and coverage
**FOLDED 2026-10-25**
Type: fact
Refs: O36, D21, G23  Source: `Cloud_Plots/20260217_cloud_data.mat` (2.73 GB, MATLAB v7.3); `docs/reports/cloud_0p5s.md`
Finding: One struct `data` with groups `nav`, `cloudphase`, `cloud`, `column`, `cprofile`, `sprofile`, `surface`, `reference`, `QA`. The whole flight is at 0.5 s (33,765 profiles; nav identical to C&S, `gps_time` differs in 37 profiles by ≤ 0.072 s). `cloudphase` is on the 845-bin HALO `z` grid (identical). The file declares no units, descriptions or fill values; the PI announced a CF-compliant version (.mat or .h5).
Impact: A 0.5 s lidar source for co-registration and for retiring most of S_repr (O36). Reader: `nurture_oe.io.halo_cloud`.

### 2026-10-25 — 0.5 s file lacks attenuated backscatter and uncertainties relayed on 2026-10-09
**FOLDED 2026-10-25**
Type: conflict-with-PK
Refs: D21, D22, O36, G22  Source: `cloudphase/*`, `cprofile/*`, `sprofile/*` vs appendix A §A.8(3) / index O36 (contents relayed 2026-10-09)
Finding: Full-grid 0.5 s fields are β_p and δ_p (532, 1064) and the phase mask only. There is no full-grid β′ (total or per channel) and no uncertainty of any kind. Channel signals (`Ppara`, `Pperp`, `Pmol`, `AB`) exist only within ≈ −375/+250 m of cloud top (`cprofile`) and ±187 m of the surface (`sprofile`, with `SNR`). The HALO lead has verbally preferred β_p over β′; a request is pending (Joe, at/after the NURTURE meeting).
Impact: If this is all that arrives: D21 liquid/optically thick rows stay on 10 s β′ (Aerosol/WV), so a 0.5 s retrieval grid means mixed resolution (10 s rows as means of F over member columns); D22 drops to rung (iii) or the placeholder diagonal S_ε; O36 is retired only for HSRL-valid (above-cloud, thin-ice) β_p rows. Affects S2 (row construction) and S4.

### 2026-10-25 — Phase mask `cp` equals C&S CloudMask/Type; class 5 labelled "NLM"
**FOLDED 2026-10-25**
Type: GAP-update
Refs: G17, G23  Source: `cloudphase/cp` vs C&S `CloudMask/Type`; strings in the PI's `20260217_F1_cloudphase.fig` / `.png`
Finding: `cp` equals `Type` on all 19,730,512 finite pixels of the flight, with identical NaN pattern. The PI's colorbar labels, in order, are WATER, MIX, ICE, SPEC, AER., NLM, i.e. 0–5 with class 5 = "NLM" (meaning not given). So the 0.5 s mask and the C&S mask are one product; C&S `Type` is 0.5 s data on the HALO grid.
Impact: Narrows G17 (class 5 has a name; definition still needed) and the G23 mask-time-resolution question (the file is per 0.5 s profile; whether the classification itself uses 10 s averaging is still unknown).

### 2026-10-25 — 0.5 s β_p,532 in m⁻¹ sr⁻¹ and δ particulate (inferred); 1064 nm not the Aerosol/WV product
**FOLDED 2026-10-25**
Type: fact
Refs: D21, G22, G23  Source: 10 s window means ([t−5, t+5), 15:30–16:00) of `cloudphase/{bsc,dep,bsc1064,dep1064}` vs Aerosol/WV
Finding: `bsc` / `532_bsc` (km⁻¹ sr⁻¹) median 1.009 × 10⁻³ (IQR 0.949–1.089; n = 26,174), so `bsc` is β_p,532 in m⁻¹ sr⁻¹. β_p-weighted `dep` matches `532_aer_dep` (median difference ≤ 0.002 in each `532_bsr` bin) and differs from volume `532_dep` by 0.15 (bsr < 2) to 0.016 (bsr ≥ 10): `dep` is particulate δ_p. `bsc1064` is 0.55 × `1064_bsc` × 10⁻³ (5–95 % 0.10–0.78); the ice colour ratio is 0.57 vs 0.97 in Aerosol/WV, so 1064 calibration or units are unknown. `dep1064` matches neither 1064 product (0.006–0.04). All fields are unscreened (|δ_p| > 1 in 4.3 % of finite pixels).
Impact: Readers carry these as `*_inferred` attributes, never `units`. The 1.009 agreement is also weak support for the [t−5, t+5) window convention (no offset sweep run). 1064 rows from this file need the PI's calibration.

### 2026-10-25 — HALO lidar is processed as body-fixed; PK's 0.3–0.5 km radar–lidar offset assumed a nadir lidar
**FOLDED 2026-10-25**
Type: conflict-with-PK
Refs: G27, O36, D23, G21  Source: `surface/offset_angle`, `surface/pitch`, `surface/roll`; CloudCube `look_vector_ecef`; `docs/reports/check_halo_cloud.json`
Finding: `offset_angle` = cos(pitch)·cos(roll) to 1e-16 on all 33,744 finite profiles (751 with |roll| > 10°; acos median 2.75°). HALO processing therefore assumes a body-fixed lidar along the body z-axis, the same model as CloudCube's look vector (follow-up §A.8(5f)). Both beams lean forward with pitch together; radar gates lie 190 m (8 km) to 520 m (1 km altitude) ahead of the aircraft nadir, and so do the lidar's, if the model holds. The "≈ 0.3–0.5 km radar–lidar along-track offset" in O36/G27 assumed a nadir lidar.
Impact: The relative geometric offset is ≈ 0 apart from lever arm and boresights (still G27). The 0.5 s co-registration uses a body-fixed lidar (equivalent to time matching). Confirm lidar mounting with J. Collins.

### 2026-10-25 — Empirical radar–lidar lag ≈ −1.1 s at 0.5 s, constant with depth
**FOLDED 2026-10-25**
Type: new-GAP-candidate
Refs: G23, G27, O36  Source: Spearman(ln β_p, Z_W) over C&S-ice pixels vs radar time shift, 15:30–16:00, 1 km bands (`scripts/check_halo_cloud.py`, `figures/check_halo_cloud_lag.png`)
Finding: Peak lag (radar time − lidar time) −1.06, −0.93, −1.09, −1.19 s at 10.5, 9.5, 8.5, 7.5 km below the aircraft (peak ρ 0.91, 0.94, 0.85, 0.70; lower bands broad, ±0.25 s). A nadir lidar would give −2.13 to −1.53 s, scaling with depth; a body-fixed lidar 0. The constant ≈ 1.1 s (≈ 260 m, two 0.5 s profiles) is timing-like: a time-stamp or clock offset between CloudCube (interval mid-time) and the HALO 0.5 s stamps (centre assumed, G23). One segment, no bootstrap; not applied (`time_offset_s = 0`).
Impact: At 0.5 s the lag misregisters radar and lidar by about two columns; at 10 s it is ≈ 10 % of a window. Either an S_repr/M term or a correction once its cause is known. Ask J. Collins (0.5 s time-stamp convention) and JPL (time base).

### 2026-11-12 — lidar_ms test K5 (milestone M6) needs the Mie fixture that is created at M9
Type: conflict-with-PK
Refs: D26, D32, O40  Source: appendix C §C.4.2 (M6: `polarization.py`, K1–K6; M9: `tests/lidar_ms/mie_fixture.py`), §C.4.3 K5
Finding: K5 runs closure C2 on `halo_ice_over_liquid` "with Mie-fixture droplet tables", but the Mie fixture is a milestone-M9 file. Decided with Joe on 2026-11-12: K5 is a skipped test at M6 (`tests/lidar_ms/test_polarization.py`) and runs once M9 exists. Every other M6 test passes; K3 and K4 use synthetic tables (a droplet-shaped table with c_perp = 0, and a Yang-shaped table with constant c_perp/p).
Impact: The C2-under-ice prediction (δ_liq in the first liquid bin < 1/10 of C1's 0.33) is still untested. That bears on O40. Appendix C §C.4.2 should list K5 under M9.

### 2026-11-12 — Under `pvc`, A_s ≥ 0.5 for isotropic near-180° scattering: C1 flag bits 0 and 2 need synthetic tables
Type: OPEN-candidate
Refs: O41, D26  Source: `lidar_ms` M6, `tests/lidar_ms/test_polarization.py` (K6); parameter scan at HALO geometry and at a spaceborne fixture geometry (range 700 km, ρ_tr 50 µrad, ρ_fov 65 µrad)
Finding: Small-angle multiple scattering can at most halve the apparent particle extinction (equivalent-medium limit), so the accumulated single-scatter fraction of a liquid layer satisfies A_s ≥ 0.5 (isotropic) and is larger with the droplet factor (F ≤ 1). The minimum in the scan was 0.50, versus 0.63 for HALO liquid at 5–60 km⁻¹. Hence δ_acc ≤ 0.17, and C1's validity limit (δ_acc > 0.45, bit 0) cannot be reached with closed-form descriptors. A scan over two liquid layers (5–60 km⁻¹, a = 5–20 µm, gaps 30 m–1.8 km, FOV 1–4 mrad) gave no negative C1 increment (bit 2): after the first layer, forward-scattered light still in view dominates the attenuated unscattered beam (M ≈ 0.5–1.3 at the second top, local A_s ≈ 0.6 vs the threshold 2√A_s − A_s ≈ 0.96). K6 therefore raises bits 0 and 2 with synthetic near-180° tables: p rising to 8 off backscatter for bit 0, a backscatter spike with F → 0 for bit 2.
Impact: Within F(x), C1 only reaches δ_acc ≈ 0.11–0.17 in liquid without wide-angle light. Hu's δ_acc up to 0.45 belongs to footprints where wide-angle multiple scattering matters, and §C.2.6 does not route that light through C1. The O41 failure mode (negative C1 increment at a second layer's top) was not observed under pvc. The validation campaign (M12) should check both against `mc`.

### 2026-11-12 — lidar_ms M3–M5: test-construction choices that appendix C leaves open
Type: fact
Refs: D27  Source: `tests/lidar_ms/test_anisotropy.py` (N2), `src/lidar_ms/derivatives.py` (`check_jvp`), `tests/lidar_ms/test_derivatives.py` (J1)
Finding: (i) N2: the 301 table nodes are spaced quadratically in γ/Θ over 0–10, giving a 4.4e-5 maximum error (appendix C measured 6e-5). Uniform spacing over 0–10 gives 4.4e-4, which fails the 2e-4 tolerance. (ii) `check_jvp` reports the normwise mismatch max|jvp − fd| / max|fd| per profile, tangent and output. Element-wise ratios are dominated by central-difference roundoff (≈1e-16 |F| / step) in elements whose derivative is far below the block maximum. Measured on `pvc`: α 6.9e-10 (appendix C: 5e-10), a 5.5e-9, β 3.1e-10, δ 3.7e-10 (`halo_liquid_30`), and up to 2.5e-8 on `halo_ice_over_liquid`. (iii) J1 compares jvp's real part with `simulate` at 1e-14 relative to each output's maximum. The real part carries the complex step's h² term (≈1e-68 here), which no element-wise test passes where the true value is 0.
Impact: None on F(x). It records how the fixed tolerances were applied.

### 2026-11-12 — CLAUDE.md locates `schema_sketch.py` and `make_refs.py` in `tests/lidar_ms/reference/`
Type: fact
Refs: D32  Source: repository layout; `tests/lidar_ms/README.md`
Finding: The two scripts are in `tests/lidar_ms/`, as in the bundle README; only `pvc_port.py`, the `.in`/`.ref` files and the licence files are in `reference/`.
Impact: Cosmetic; update the CLAUDE.md "Current work" line at the next refresh.
