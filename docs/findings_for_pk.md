# findings_for_pk.md — findings to fold into Project Knowledge

Append-only. Template and rules: see `CLAUDE.md`. Joe and Fable fold these into `index.md` / `appendix_A_nurture_datasets.md` and mark each entry `FOLDED <date>`.

### 2026-10-14 — CloudCube Ka sensitivity attribute differs between 02-17 granules
Type: conflict-with-PK
Refs: G1, G24, O15, O31  Source: global attr `radar_sensitivity_ka_dbz` in all three 02-17 CloudCube granules
Finding: The 15:13–15:20 and 15:21–18:59 granules declare Ka **−16.32 dBZ @ 1 km**; the 19:00–19:26 granule declares −19.26, the value in appendix A §A.4 and O15. W (−32.64) and G (−37.58) agree in all three. The delivered Ka floor moves with the attribute: the 0.1th percentile of $z_m - Z_{\min}(r)$ is −9.0 dB for both granules (6–9 km range), so the ~3 dB change looks real, not a typo.
Impact: $Z_{\min}^{Ka}(r)$ in the O15 coverage check (4.0 % vs 10.9 % on 15:30–16:00); the Ka random-error term of the interim radar $\mathbf{S}_\epsilon$ (O31); any Ka presence mask. The sensitivity must be read per granule. Ask JPL at the STM what changed at 19:00.

### 2026-10-14 — Delivered zm extends 5–11 dB below the stated Z_min(r) in every band
Type: GAP-update
Refs: G24, G1, O21, O31  Source: `zm_{ka,w,g}` vs `radar_sensitivity_*_dbz`, granules 15:21–18:59 (15:30–16:00) and 19:00–19:26
Finding: The 0.1th percentile of $z_m - Z_{\min}(r)$ (from each granule's own attribute) is Ka −7.7 to −9.0, W −9.8 to −10.9 and G −6.3 to −7.0 dB at 6–11.5 km range. The offset is about constant with range, so the $20\log_{10} r$ form holds. So the documented "sensitivity @ 1 km" is not the detection threshold of the delivered data; the SNR definition behind it is still unknown.
Impact: The SNR term $(1+\mathrm{SNR}^{-1})$ of the interim radar $\sigma_{\text{rand}}$ (§A.4(iv)) cannot use the attribute as the noise floor without the SNR definition. It also affects the O21 non-detection threshold and the "threshold slightly above $Z_{\min}$" rule (index §2.3), and which definition the O15 test uses ("finite" 12.8 % vs "> Z_min" 4.0 %).

### 2026-10-14 — CloudCube `alt` is ellipsoidal; HALO `z`/`gps_alt` are MSL (20–44 m apart)
Type: fact
Refs: G1, O36; appendix A §A.11(b)  Source: CloudCube `alt:long_name` ("above WGS84 ellipsoid"); HALO `Nav_Data/gps_alt`; MetNav `GPS_Altitude_MSL`
Finding: HALO `gps_alt` equals MetNav MSL altitude (median +0.1 m), so HALO `z` is MSL. HALO − CloudCube altitude is a median 38.1 m (1–99 %: 18.7–44.0 m) in the 15:21–18:59 granule and 20 m in the others: the geoid undulation along track. Gates referenced to HALO MSL altitude put the W surface peak a median −11 m (15:30–16:00) and +7 m (19:00–19:25) from C&S `altSurface`; with CloudCube `alt` the residuals are −40 m and −13 m.
Impact: Radar–lidar vertical registration (about two gates if ignored); DWR alignment with lidar layers and the liquid mask; the S2 retrieval-grid choice. Co-registration uses HALO MSL altitude (`coreg.cloudcube_to_halo_grid`, `aircraft_alt="halo"`). The antenna lever arm is still unknown.

### 2026-10-14 — CloudCube L2 header facts (first header inspection)
Type: GAP-update
Refs: G1, G24  Source: `docs/reports/headers/NURTURE-CLOUDCUBE_*.txt`
Finding: Three 02-17 granules: 15:13:00–15:20:10, 15:21:58–18:59:37 (4.7 GB) and 19:00:00–19:26:25. Dims are `(time, range)` and `look_vector_ecef (time, nxyz)`, not `(range, time)`/`(3, time)` as in the readme. There are 751 gates 19.986 m apart from −999.3 m (r = 0 at index 50), and rays every 0.10044 s at the interval mid-time. Fills: `zm`/`katten`/`s0` use `_FillValue = NaN`; nav variables have no `_FillValue` (netCDF default 9.97e36); `landmask` uses −2147483647. `katten_ka`/`katten_w` are finite at 22–25 % of r < 0 gates. Arrays are contiguous and uncompressed. The off-nadir angle from `look_vector_ecef` matches `incidence_angle` to 0.1°. CloudCube nav/attitude are MetNav interpolated (lag 0 s).
Impact: Gate spacing for $\mathbf{K}$/$F$ range weighting (the range *resolution* is still unknown, G1); co-registration; reader design.

### 2026-10-14 — HALO 10 s profile averages the C&S sub-profiles in [t − 5 s, t + 5 s)
Type: GAP-update
Refs: G23, O36  Source: HALO `time` vs C&S `Nav/gps_time`; `IAB532` vs column `532_total_attn_bsc`
Finding: HALO profile spacing is irregular (median 10.0 s, range 9.04–23 s), and `time:long_name` reads "instant time in UTC". Over 273 heterogeneous profiles, the centered window [t − 5, t + 5) gives the best agreement between windowed C&S `IAB532` and column $\beta'_{532}$ (correlation 0.52, vs 0.41 for [t, t + 10)). It holds 20 sub-profiles in >95 % of profiles. Alignment is resolved to about ±1 s. HALO `lat`/`lon` are positions at t (they match C&S and MetNav within 1 s), not window means. With this rule the §A.10.5 counts reproduce within 0.4 % (10,439 liquid / 62,700 ice vs 10,395 / 62,567); the 2026-10-06 window convention was not recorded.
Impact: The window rule for all co-registration and for $\mathbf{S}_{\text{repr}}$ from mask heterogeneity (O36). PI confirmation wanted.

### 2026-10-14 — HALO `imu_roll` is neither a sample at t nor a centered 10 s mean
Type: GAP-update
Refs: G21, O28  Source: HALO `Nav_Data/imu_roll`, `imu_pitch` vs MetNav `Roll_Angle`, `Pitch_Angle` (1 Hz)
Finding: Over the 50 turning profiles (|roll| > 5°), RMS(HALO − MetNav roll) is 6.2° against MetNav at t and 5.3° against the centered 10 s mean. The best fit is a 10 s mean centered 6 s *after* t (2.3°). The clocks agree within 1 s (positions), so this is attitude-specific. `imu_pitch` − MetNav pitch has a median of +0.29°.
Impact: The HALO incidence angle for D19 flags and $\mathbf{M}$ sizing: do not take $\theta_{\text{inc}}$ from HALO `imu_*` until the PI confirms their definition. MetNav 1 Hz (window-averaged) is the defensible source.

### 2026-10-14 — `h2o_mmr_v` is not truncated 100 m above cloud top in most profiles
Type: GAP-update
Refs: O32, D7  Source: `DataProducts/h2o_mmr_v`, `cloud_top_height_935`
Finding: Of the 1,329 profiles with `cloud_top_height_935`, only 25 % have their lowest valid DIAL bin ≥ 85 m above it. The median lowest valid bin is 495 m *below* `cloud_top_height_935` (5–95 %: −3,598 to +696 m). This is consistent with §A.10.5's in-liquid values.
Impact: $\ln q_v$ rows: the in-cloud exclusion (§A.10.6) must use our own cloud mask, not the product's documented truncation.

### 2026-10-14 — `1064_dep_cloud_screened_unc` filled where its product is NaN
Type: conflict-with-PK
Refs: G22, G23  Source: `DataUncertainty/1064_dep_cloud_screened_unc` vs `DataProducts/1064_dep_cloud_screened`
Finding: Five of the six `_cloud_screened_unc` fields are finite on exactly the pixels of their screened product, as appendix A §A.10.4 states for all of them. `1064_dep_cloud_screened_unc` is the exception: it is finite on 10,971 pixels where `1064_dep_cloud_screened` is NaN (899,377 both finite).
Impact: Minor. Any 1064 $\delta$ row must take its mask from the product, not from `_unc`. Add to the HALO metadata list (§A.10.3).

### 2026-10-14 — State and variable inventory of the 02-17 Aerosol/WV granule
Type: fact
Refs: G23; appendix A §A.8(1), §A.10.2  Source: `State/*`, full header dump
Finding: `State_Type` = 2 (MERRA-2). Units: Temperature K, Pressure **atm** (0.168–1.051), Number_Density m⁻³, Relative_Humidity %, U/V m/s. Only Nav/State/UserInput and root `time`/`z` lack dimension scales; `DataProducts`/`DataUncertainty` have `(time, z)` attached. Every variable named in §A.10.2 is present. 23 more variables are not named there (`*_bsc_Sa` (+`_unc`), `532_AOT_from_bsc*`, several `_unc`, nav, `DataProducts/Altitude` = "legacy … same as z"). `Units `/`units` differ semantically only for `Angstrom_Spherical`. A complete diff against 01-24 needs the 01-24 dump, which is not in the repo.
Impact: Closes the §A.8(1) `State_Type` item; $T$, $p$ (convert atm) for $F(\mathbf{x})$ and $k_{\text{gas}}$.

### 2026-10-14 — Ka coverage on the 02-17 15:30–16:00 TPV transect
Type: fact
Refs: O15, D12, O21  Source: `docs/reports/task_d_coreg.md`; `ka_coverage_*_20260217_1530-1600.csv`
Finding: 18,944 ice pixels (C&S classes 1–2 in ≥ 50 % of sub-profiles, all T < 0 °C) have W and G detected. Ka also detects in **4.0 %** with the granule's −16.32 dBZ, 10.9 % with PK's −19.26, and 12.8 % with any finite `zm_ka`. All Ka detections are in the lowest about 2 km (10–12 km slant range: 14 / 36 / 42 %); coverage is 0 % at 6–8 km range. With 500 m clearance above the surface: 2.2 / 8.6 %.
Impact: O15 condition 3. Ka adds rows only in near-surface snow on this transect, not where the density/riming signal is wanted aloft; elsewhere Ka contributes only via the O21 non-detection penalty. A second segment (19:00–19:25, Ka −19.26) would test generality.

### 2026-10-14 — Gas path in DWR_W−G near the surface
Type: fact
Refs: G24, O32; index §4.1  Source: `katten_w`, `katten_g` (15:30–16:00), display correction $2\int k\,dr$
Finding: The median two-way gaseous attenuation at 1 km altitude (r ≈ 11 km) is 1.3 dB at W and 5.9 dB at G, so about 4.7 dB of the near-surface DWR_W−G is gas alone. `katten_g` reaches 1.9 dB km⁻¹ (99th percentile).
Impact: STM slide framing; the $k_{\text{gas}}$ consistency check; DWR in $\mathbf{y}$ must be modeled with the full gas path in $F(\mathbf{x})$, not corrected beforehand.
