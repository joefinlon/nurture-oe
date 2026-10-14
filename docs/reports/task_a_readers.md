# Task (a) — NURTURE readers and data-level checks (2026-02-17)

*2026-10-14. Code: `src/nurture_oe/io/`. Numbers: `scripts/run_task_a_checks.py` → `docs/reports/task_a_checks.json`. Header dumps: `docs/reports/headers/`. Findings logged in `docs/findings_for_pk.md`.*

## Readers

| Module | Function | Returns |
|---|---|---|
| `io/halo_aerowv.py` | `read_halo_aerowv(path, groups, variables, time_range)` | (time 1689, z 845) flat Dataset with verbatim names (`"532_ext"`). Nav/State/UserInput lack dimension scales and are mapped by length; the reader raises if a length is ambiguous. Warns on a `Units `/`units` conflict. |
| `io/halo_cloudsurface.py` | `read_halo_cloudsurface(path, variables, time_range)` | (profile 33,765) with dims `z` 845, `wcp_dz` 501 and `layer` 10. (N,1) arrays are squeezed. Names with spaces are kept. `000_Readme` goes into attrs. |
| `io/cloudcube.py` | `read_cloudcube(path, bands, time_range, profiles)` | (time, range 751), lazily sliced. Parses `radar_sensitivity_*_dbz` into `zmin_1km_dbz_*`, and warns when the value differs from PK. |
| `io/metnav.py` | `read_metnav(path)` | 1 Hz ICARTT. Missing values become NaN. LLOD/ULOD get boolean companions (none present on 02-17). |

All readers return native grids and units, NaN fills and every attribute verbatim. Each adds `source_file`/`hdf5_path` and a `datetime64` `time` coordinate, keeping the raw value as `time_native`. The ICARTT and C&S files do not store the flight date, so it comes from the data-date line and the file name respectively.

**HDF5 vs CF fills.** The h5py `fillvalue` (0.0 in the HALO files) is the HDF5 allocation fill, not a CF fill. The readers ignore it and use only `_FillValue`/`missing_value`. CloudCube nav variables have no `_FillValue`; the netCDF default (9.97e36) is converted to NaN.

## HALO Aerosol/WV

- **Grid.** 1,689 × 845. `z` runs from −300.0 to 12,351.2 m in 14.990 m steps, **ascending** (closes data_map check 7). Time spacing is irregular: median 10.0 s, range 9.04–23 s. `time:long_name` reads "instant time in UTC".
- **Variable diff.** All variables named in appendix A §A.10.2 are present. 23 variables are present but not named there: the `_bsc_Sa` family, `532_AOT_from_bsc*`, several `_unc` fields, Nav fields and `DataProducts/Altitude`. The latter is described as a "legacy altitude array … same as z". A complete diff against the 01-24 file needs the full 01-24 dump (`HALO_L2_description.txt`), which is not in the repo.
- **`State`.** `State_Type` = **2 (MERRA-2)**. Units: Temperature K, **Pressure atm** (0.168–1.051), Number_Density m⁻³, Relative_Humidity %, U/V m/s.
- **`Units ` vs `units`.** These are notation differences everywhere except `Angstrom_Spherical` (`sr` vs `1`). That one raises the only reader warning.
- **PK regressions (§A.10.5).**
  - Cloud pixels more than 2.5 km below the aircraft: 10,439 liquid and 62,700 ice, against PK's 10,395 / 62,567 (within 0.4 %). The window convention used on 2026-10-06 was not recorded; no variant reproduces PK exactly.
  - Fill fractions (liquid / ice): unscreened 532 0.909 / 0.804; 1064 `total_attn_bsc` 0.172 / 0.479; screened fields 0.055 / 0.504; `532_ext` 0.042 / 0.457; `Aerosol_ID` 0.013 / 0.345; `h2o_mmr_v` 0.665 / 0.932. All match PK.
- **`_unc` masks.** Five of the six `_cloud_screened_unc` fields are filled on exactly the pixels of their screened product. **`1064_dep_cloud_screened_unc` is the exception:** it is finite on 10,971 pixels where `1064_dep_cloud_screened` is NaN.
- **`Aerosol_ID` vs C&S.** Classes present: 1, 3, 4, 5, 6, 7, 9. P(C&S ice | `Aerosol_ID`=1) = 0.47. P(`Aerosol_ID`=1 | C&S ice) = 0.34.
- **NaN census by region** (finite fraction; full table in the JSON):

  | Region | `532_total_attn_bsc` | `532_dep`, `532_bsc` | `h2o_mmr_v` |
  |---|---|---|---|
  | Above the aircraft | 0 | 0 | 0 |
  | 0–0.5 km below | 0 | 0 | 0 |
  | 0.5–2.5 km below | 0.47 | 0.94 | 0.71 |
  | Below the DEM | 0 | 0 | 0.002 |

- **`h2o_mmr_v` truncation `[O32]`.** The README's "100 m above cloud top" holds in only 25 % of the 1,329 profiles that have `cloud_top_height_935`. The median lowest valid DIAL bin is **495 m below** `cloud_top_height_935` (5th–95th percentile: −3,598 to +696 m). In-cloud $q_v$ rows must be excluded by our own mask, not by the product's truncation.

## HALO Cloud & Surface

- 33,765 profiles at a median 0.501 s, with gaps up to 13.5 s. `z` is identical to Aerosol/WV.
- **`CloudMask/Type` class counts:**

  | Class | Pixels |
  |---|---|
  | 0 (water) | 228,635 |
  | 1 | 399,206 |
  | 2 (ice) | 873,004 |
  | 3 (HOI) | 28,122 |
  | 4 (aerosol) | 18,188,639 |
  | 5 (undocumented) | 12,906 |

- `WaterCloudProfile/Altitude`: 501 bins at 1.249 m, minimum −374.7 m (as PK).
- **10 s window convention** (needed for co-registration; `[G23]`). A HALO profile at time t averages the C&S sub-profiles in **[t − 5 s, t + 5 s)**. Over the 273 heterogeneous profiles, correlation with column $\beta'_{532}$ is 0.52, against 0.41 for [t, t + 10) and 0.42 for [t − 10, t). The window holds 20 sub-profiles in more than 95 % of profiles. The alignment is resolved to about ±1 s. This is empirical; the PI should confirm it.

## Attitude: HALO `Nav_Data` vs MetNav `[G21]`, `[O28]`

- MetNav is 1 Hz with no gaps (17,733 rows, 14:46:19–19:41:51). Its `Pitch_Angle` is up-positive and `Roll_Angle` right-positive. It has no LLOD/ULOD values.
- **Clocks agree.** HALO and C&S `lat` lead MetNav by 0.5–1 s. CloudCube `lat`, `roll` and `pitch` are MetNav interpolated (best lag 0 s, residual roll 0.003°).
- **`imu_roll` matches neither hypothesis.** RMS differences over the 50 turning profiles (|roll| > 5°):

  | HALO `imu_roll` compared with | RMS in turns |
  |---|---|
  | MetNav sampled at t | 6.2° |
  | MetNav 10 s mean centered at t | 5.3° |
  | MetNav 10 s mean centered at t + 6 s (best fit) | 2.3° |

  The clocks agree, so this is an attitude-specific timing or averaging convention (or a different IMU). Don't use HALO `imu_*` for incidence angle until the PI confirms what they are. MetNav, or CloudCube's copy of it, is the defensible 1 Hz source.
- **Pitch offset.** `imu_pitch` − MetNav pitch has a median of +0.29°: a possible IMU or mount offset.

## CloudCube (three granules)

| | 15:13:00–15:20:10 | 15:21:58–18:59:37 | 19:00:00–19:26:25 |
|---|---|---|---|
| rays | 4,283 | 130,027 | 15,787 |
| `radar_sensitivity_ka_dbz` | **−16.32** | **−16.32** | −19.26 (PK value) |
| HALO `gps_alt` − CC `alt` (median) | 19.9 m | 38.1 m (1–99 %: 18.7–44.0) | 20.0 m |

- **Dims.** `(time, range)` and `look_vector_ecef (time, nxyz)`. The readme says `(range, time)` and `(3, time)`.
- **Range grid.** 751 gates, 19.986 m apart, from −999.3 m; r = 0 is at index 50. Rays every 0.10044 s, with `time` at the interval mid-time.
- **Fills.** `zm`, `katten` and `s0` use `_FillValue = NaN`. Nav has no `_FillValue`. `landmask` uses the int fill −2147483647 (not seen). `katten_ka`/`katten_w` are finite at 22–25 % of r < 0 gates; `katten_g` at 0 %.
- **Vertical datum.** `alt:long_name` says "above WGS84 ellipsoid". HALO `gps_alt` equals MetNav `GPS_Altitude_MSL` (median +0.1 m), so HALO `z` is MSL. The 20–44 m difference is the geoid undulation along track.
- **Range → altitude.**
  - Off-nadir angle from `look_vector_ecef` agrees with `incidence_angle` to 0.1° (median).
  - W surface peak altitude − C&S `altSurface`, median: −11 m on 15:30–16:00 and +7 m on 19:00–19:25 using the HALO MSL reference. Using CloudCube's ellipsoidal `alt` instead gives −40 m and −13 m.
- **Detection floor vs stated sensitivity `[G24]`.** In every band and both inspected segments, finite `zm` reaches 5–11 dB below $Z_{\min}(r)$ computed from the granule's own attribute. That is the 0.1th percentile of $z_m - Z_{\min}(r)$:

  | Band | 6–9 km | 9–11.5 km |
  |---|---|---|
  | Ka | −9.0 (15:30), −8.7 (19:00) | −8.5, −7.7 |
  | W | −10.9, −10.6 | −9.9, −9.8 |
  | G | −7.0, −6.9 | −6.3, −6.4 |

  The offset is roughly constant with range, so the $20\log_{10} r$ form holds. The Ka offset is the same relative to −16.32 and to −19.26, so the Ka floor itself moved by about 3 dB between granules. The attribute change looks real.
- **Gaseous attenuation (`katten`, 15:30–16:00), dB km⁻¹, 1st / 50th / 99th percentile:**

  | Band | 1 % | 50 % | 99 % |
  |---|---|---|---|
  | Ka | 0.004 | 0.022 | 0.117 |
  | W | 0.004 | 0.023 | 0.32 |
  | G | 0.004 | 0.048 | 1.92 |

## Not done / open

- Complete Aerosol/WV diff against 01-24: needs the 01-24 dump.
- **Radar near-field.** Gates are used from r ≥ 0 (`min_range` parameter). The radar near-field limit is undocumented `[G1]`.
- **Mount offsets.** The lever arm between the radar antenna and the HALO reference is unknown. The −11 m / +7 m surface residuals are within ±1 gate.
