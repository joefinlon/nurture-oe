# HALO 0.5 s cloud file: reader, 0.5 s co-registration, review of the 2026-10-20 display file

2026-10-25 · Claude Code (nurture-oe) · base: `index.md` v1.12, appendix A (2026-10-20)

Scripts: `scripts/check_halo_cloud.py` → `docs/reports/check_halo_cloud.json`, `figures/check_halo_cloud_lag.png`, `headers/20260217_cloud_data.txt`; `scripts/build_coreg_0p5s.py` → `<cache_dir>/coreg_cloud_20260217_1530-1600_0p5s.nc`; `scripts/stm_fourpanel.py --0p5s` → `figures/stm_20260217_1530-1600_0p5s.png`.

## Summary

- **The file.** `Cloud_Plots/20260217_cloud_data.mat` (2.73 GB, MATLAB v7.3 = HDF5, delivered 2026-10-23) holds the whole flight at 0.5 s.
  - It has the same 33,765 profiles and nav as C&S and the same 845-bin `z` as Aerosol/WV.
  - It declares **no units, descriptions or fill values**.
- **What it contains:**
  - 0.5 s **particulate β_p,532 and δ_p,532** on the full grid. β_p is in m⁻¹ sr⁻¹ (inferred); δ_p is particulate (inferred).
  - 1064 nm counterparts.
  - The phase mask `cp`, identical to C&S `CloudMask/Type`.
- **What it does not contain:**
  - No **full-grid attenuated backscatter**.
  - No **uncertainty** of any kind.
  - Channel signals (∥, ⊥, molecular) and `AB` exist only within about ±300 m of cloud top and ±190 m of the surface.
- **Lidar pointing.** HALO's own `surface/offset_angle` equals cos(pitch)·cos(roll) on every profile. HALO therefore processes the lidar as body-fixed, like the CloudCube look vector. The radar–lidar relative geometric offset is then about 0, **not 0.3–0.5 km**: PK G27/O36 assumed a nadir lidar.
- **Timing.** An empirical radar–lidar lag of **−0.9 to −1.2 s** (radar earlier; ≈ 220–280 m) is constant with depth below the aircraft. That pattern is timing-like, not geometric. It is unexplained and not applied.
- **Co-registration.** The 0.5 s file is rebuilt from native 0.5 s lidar fields. The held-10 s lidar values of 2026-10-20 are removed, while the 10 s path is unchanged.

## 1. File inventory (`io/halo_cloud.py`)

`read_halo_cloud(path, group=…, variables=None, time_range=None)` reads one group plus `nav`, slicing lazily by time. Its conventions:
- Names are the path below the group with `/` → `_`.
- Values are native: no units conversion, no screening.
- No `units` attribute is set. Evidence-backed inferences go in `units_inferred` / `meaning_inferred` / `evidence` (table `INFERRED`).
- When the announced CF-compliant file declares attributes, those govern, and the reader warns on any disagreement with `INFERRED`.

| Group | Dims | Contents (names verbatim) | Status |
|---|---|---|---|
| `nav` | profile | `gps_time` (h), `gps_alt`, `gps_lat`, `gps_lon` | identical to C&S `Nav/*`; `gps_time` differs in 37 profiles by ≤ 0.072 s |
| `cloudphase` | profile × `z` (845) | `bsc`, `dep`, `bsc1064`, `dep1064`, `cp`, `S0`, `cT2p`, `X1064m_syn`; `alt` = HALO `z` | §2 |
| `cloud` | profile | `IAB` (≡ C&S `IAB532`), `IAB1064_uncal`, `LID`, `LID1064`, `MT2p`, `norm` | C&S-like column values |
| `column` | profile, `layer` (10), `od_coord` (200), `od_level` (5) | `cth`, `cths`, `cbh`, `cbhs`, `columnclass`, `Sc`, `Sc_thin`, `ext`, `cprofilewaterfraction`, `tempLevel{0C,Minus20C,Minus40C}{Low,High}Alt`, `OD_{ext,dep,OD_depth}` | undocumented |
| `cprofile` | profile × `wcp_dz` (501, 1.249 m; = C&S `WaterCloudProfile/Altitude`) | `AB`, `AB1064`, `IAB`, `IAB1064`, `LID`, `LID1064`, `Ppara`, `Pperp`, `Pmol`, `Ptotal`, `Ptotal1064`, `beta_m`, `moltrans2p`, `dep`, `dep1064`, `ext`, `ext_HSRL`, `ext_thin` | cloud-top-relative only |
| `sprofile` | profile × `sfc_dz` (301, ±187 m) | `Ppara`, `Pperp`, `Pmol`, `Ptotal` (+1064), `SNR` | surface-relative only |
| `surface` | profile | 32 fields incl. `altSurface`, `DEMalt`, `pitch`, `roll`, `offset_angle`, `windSpeed`, `reflectance*` | `offset_angle` §2 |
| `reference` | profile | `Altitude`, `CT2p`, `CT2p_1064`, `GR`, `N_ref`, `lambda`, `lambda_unc`, `ref_unc`, `weight` | undocumented (calibration reference?) |
| `QA` | profile; `block_*` 34 × 201 | saturation counts/maxima per channel; channel cross-correlations | not interpreted |

## 2. What the data establish (15:30–16:00 unless stated)

- **Grids.** `cloudphase/alt` equals HALO `z` exactly; `cprofile/alt` equals C&S `WaterCloudProfile/Altitude`.
- **Phase.** `cp` equals C&S `CloudMask/Type` on all 19,730,512 finite pixels of the flight, with the same NaN pattern. The colorbar strings of the PI's `20260217_F1_cloudphase.fig` are, in order, WATER, MIX, ICE, SPEC, AER., NLM. This labels class 5 "NLM" (`[G17]`); the meaning of NLM is unknown.
- **β_p,532 units.** 10 s window means ([t−5, t+5), fixed) of `bsc` divided by Aerosol/WV `532_bsc` (km⁻¹ sr⁻¹) give a median of **1.009 × 10⁻³** (IQR 0.949–1.089, n = 26,174). So `bsc` is β_p in m⁻¹ sr⁻¹. The PI's figure label is β_p (Mm⁻¹ sr⁻¹).
- **δ_532 is particulate.** β_p-weighted 10 s means of `dep` were compared with both Aerosol/WV depolarization products:

  | `532_bsr` bin | − `532_aer_dep` | − `532_dep` (volume) |
  |---|---|---|
  | < 2 | +0.002 | +0.152 |
  | 2–10 | +0.001 | +0.061 |
  | ≥ 10 | +0.001 | +0.016 |

  `dep` is unscreened: |δ| > 1 in 4.3 % of finite pixels (clear air).
- **1064 nm is not the Aerosol/WV 1064 product.** `bsc1064` window means are **0.55 ×** `1064_bsc` × 10⁻³ (median; 5–95 % 0.10–0.78). The ice colour ratio `bsc1064/bsc` is 0.57 (IQR 0.54–0.60), against 0.97 (0.93–1.03) for Aerosol/WV `1064_bsc/532_bsc`. Calibration and units are unknown. C&S-style `_uncal`/`_noncal` names elsewhere hint at an uncalibrated 1064 channel. `dep1064` differs from both `1064_dep` and `1064_aer_dep` by 0.006–0.04, so whether it is particulate is undetermined.
- **Undocumented fields.** `S0` (median 21, 1–99 % −402 to 441) and `cT2p` / `X1064m_syn` (≈ 3×10¹⁶ / 1.5×10¹⁷, signal-like magnitudes) are not interpreted.
- **Lidar pointing.** `surface/offset_angle` = cos(`pitch`)·cos(`roll`) to 1e-16 on all 33,744 finite profiles, including 751 with |roll| > 10°. acos gives 2.75° median (1–99 %: 1.1–20.2°). HALO's processing assumes a body-fixed, body-z lidar, which is the same model JPL uses for `look_vector_ecef` (follow-up §A.8(5f)).

## 3. Radar–lidar registration at 0.5 s

**Radar gate offsets from the aircraft nadir.** Along-track and cross-track offsets of the radar gates (median and 5–95 %), with 5 rays per pixel:

| Altitude | Along-track (ahead) | Cross-track | Time equivalent |
|---|---|---|---|
| 1 km | 520 m (489–613) | −20 m (−83 to +60) | ≈ 2.2 s |
| 2 km | 473 m | −18 m | |
| 4 km | 379 m | −15 m | |
| 6 km | 285 m | −11 m | |
| 8 km | 190 m | −7 m | ≈ 0.8 s |

**Lag diagnostic** (`figures/check_halo_cloud_lag.png`). The method:
- Ice pixels only (`cp` ∈ {1, 2}), correlating ln β_p with Z_W.
- The radar is re-aggregated by ray time with every ray shifted by the lag, so lags are not tied to the 0.5 s profile step.
- The peak is the vertex of a parabola through the best three lags.

| Band | Depth below aircraft | Observed peak (radar − lidar) | Peak ρ | Nadir lidar would give | Body-fixed lidar gives |
|---|---|---|---|---|---|
| 1–2 km | 10.5 km | −1.06 s | 0.907 | −2.13 s | 0 |
| 2–3 km | 9.5 km | −0.93 s | 0.942 | −1.93 s | 0 |
| 3–4 km | 8.5 km | −1.09 s | 0.848 | −1.73 s | 0 |
| 4–5 km | 7.5 km | −1.19 s | 0.696 | −1.53 s | 0 |

**Reading:**
- A geometric pointing difference grows with depth, but the observed lag does not. A constant ≈ −1.1 s is a timing signature, about 2 lidar profiles or ≈ 260 m. It could come from a clock or time-stamp convention offset between CloudCube (interval mid-time) and the HALO 0.5 s stamps (centre assumed, `[G23]`). It could also be a boresight difference that is constant in time rather than angle, which is unlikely.
- The nadir-lidar footprint correction would over-shift the radar by 0.4–1.1 s.
- The lower bands' curves are broad, so their peaks are uncertain to ≈ ±0.25 s. This is one segment, and there is no bootstrap.

**Code.** `coreg.cloudcube_to_halo_grid(..., assign="footprint")` assigns each (ray, z) gate to the window that contains the time the lidar footprint reaches it:

$$t_{\text{eff}} = t_{\text{ray}} + \frac{\text{along} - \text{along}_{\text{lidar}}}{v_g} + \text{time\_offset\_s}$$

- Options for `lidar_pointing`:
  - `"body_fixed"` (default): along_lidar = along;
  - `"nadir"`;
  - a fixed along-track tangent.
- Joe's choice for the builder: body-fixed lidar and `time_offset_s = 0`. On 02-17 this is identical to time matching. The 1.1 s lag is recorded, not applied.
- `assign="time"` (default) is unchanged. A regression against the committed `coreg.py` on 15:40–15:45 gives identical arrays for every variable.
- Footprint mode needs non-overlapping windows (`subprofile_windows`). It adds the diagnostics `along_track_offset_mean`, `cross_track_offset_mean` and `footprint_dt_mean`.

**OE note.** Cross-track offsets (±60–85 m at 1 km, 5–95 %) cannot be corrected on a curtain. With 0.5 s columns about 120 m long, they are a representativeness term at the scale of a column. Like the remaining 1.1 s, they belong in $\mathbf{S}_{\text{repr}}$ / $\mathbf{M}$ when S2 builds 0.5 s rows (`[O36]`).

## 4. The 0.5 s co-registered file (`scripts/build_coreg_0p5s.py`)

- **General CLI.** `--start/--end` select the segment. CloudCube granules are picked from their `_S…_E…` file-name stamps (`config.cloudcube_granules`). Each granule is aggregated on its own, and every profile is tied to one granule, with `cc_granule` and `zmin_1km_dbz_{ka,w,g}` per profile. This was tested across 18:59:37 / 19:00:00: profiles carry −16.32 then −19.26 dBZ, and the gap stays empty.
- **Contents.**
  - Lidar, native 0.5 s: `bsc532`, `dep532`, `bsc1064`, `dep1064`, `phase`, plus `liquid` / `ice` from `phase` and `surface_alt_lidar`. Names are mapped in one dict, `LIDAR_FIELDS`, so a renamed CF file touches one line.
  - Radar: as for 10 s (`zm_*`, `*_all`, `*_detect_frac`), with `n_rays` and the offset diagnostics as (time, z).
  - Display only: `zgc_*` and `dwr_wg_gc`.
- **Size.** 15:30–16:00 gives 3,589 profiles × 845, with 5 rays per pixel (4–6 at the edges). The build takes 7 s.
- **Figure.** `stm_20260217_1530-1600_0p5s.png` shows:
  - β_p,532 in km⁻¹ sr⁻¹, converted from m⁻¹ sr⁻¹ ×10³ for display only;
  - δ_p,532, shown where β_p ≥ 10⁻³ km⁻¹ sr⁻¹;
  - Z_W;
  - DWR with the 0.5 s liquid contour.
- **Superseded file.** `coreg_20260217_1530-1600_0p5s.nc` (held 10 s lidar, 2026-10-20) is replaced and can be deleted from `cache_dir`.

## 5. Item 3: the 2026-10-20 tweak

| Piece | Decision | Why |
|---|---|---|
| Held 10 s β′/δ_v/T at 0.5 s, `halo_hold`, `HOLD_EXTEND` | **removed** | Copying one 10 s profile into 20 rows makes perfectly correlated duplicates. A diagonal $\mathbf{S}_\epsilon$ would count them 20 times (`[O29]`, worse). It was display-only; native 0.5 s β_p/δ_p now exist. |
| `coreg.containing_window` (+ 2 tests) | **removed** | It is a one-parent map and cannot represent sub-profiles shared by overlapping 10 s windows. For a 10 s row on a 0.5 s retrieval grid ($y = \frac1n\sum_k F(\mathbf{x}_k)$) the member map is `coreg._bounds`, now documented as such. |
| `coreg.subprofile_windows` | **kept** | It builds the 0.5 s windows. |
| 10 s path (`run_task_d.py`, `halo_windows`, time mode) | **unchanged** | It remains the y-row prototype. Whether S2 retrieves on 10 s or 0.5 s columns is open (`[O36]`). |

**Caveat for S2.** Per `[D21]`, 0.5 s β_p rows are valid only where the HSRL inversion holds (above cloud, between layers, thin ice). Liquid and optically thick rows need β′ + δ. If 0.5 s β′ never arrives, a 0.5 s grid means mixed resolution (β′ at 10 s in liquid) and the averaging operator above.

## 6. Requests for the HALO lead (proposed for `team_followups.md`)

1. **0.5 s attenuated backscatter at 532 nm** (total, or ∥ and ⊥) on the 845-bin grid. In dense cloud, multiple scattering biases `532_bsc` (PI, 2026-10-06). D21 and F(x) simulate β′ directly, so β_p cannot replace it in liquid.
2. **Per-channel random uncertainty**, split into signal and background terms, or the total σ (D22 rungs i–ii). Failing that, NSF and background per channel (rung iii).
3. **CF metadata**: units, descriptions and fill values for every variable (announced).
4. **Definitions**: whether `dep` / `dep1064` are particulate; the calibration of `bsc1064` (0.55× Aerosol/WV); class 5 "NLM"; `S0`, `cT2p`, `X1064m_syn`, `MT2p`, `norm`; the `reference` group.
5. **Timing and pointing**: where the 0.5 s time stamp sits in its interval. Is the lidar body-fixed (`offset_angle` suggests so)? Is any boresight or tilt known? The ≈ 1.1 s lag against CloudCube also goes to JPL (their time stamp is the interval mid-time).
6. **Other flights**: when the CF version, and files for other flights, will be available.

## 7. Open after this session

- Whether the 1.1 s lag holds on other segments and flights, or with a bootstrap. This is not commissioned.
- `[G23]` window-offset sweep with the 0.5 s β_p: an unapproved §A.8(6) candidate. The reader makes it cheap.
- The new CF file: re-run `scripts/check_halo_cloud.py` and `pytest tests/test_halo_cloud.py`. The reader warns where declared units contradict `INFERRED`.
