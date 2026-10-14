# Fold map for `index.md` v1.11 / appendix A / handoff (from Claude Code, 2026-10-14)

**Purpose.** Lets a Project-attached Chat session apply this Code session's results to Project Knowledge in one batched wrap-up pass (`index.md` §0).
**Sources:** `docs/findings_for_pk.md` (eleven entries dated 2026-10-14), `docs/reports/task_a_readers.md`, `docs/reports/task_d_coreg.md`, `docs/reports/figures/stm_20260217_1530-1600.png`.
**Rules.** No register numbers were minted here. Where a new entry seems warranted, it is marked *candidate*; §11 assigns numbers. Line numbers refer to the 2026-10-14 snapshot in `docs/pk_snapshot/` and are hints only; the quoted phrases are the anchors.
**Not checked:** `appendix_B_stm_figures.md`, which is not in the snapshot. Its "RF17 context facts" should be checked against the findings below (RF17 = the 2026-02-17 flight).

After folding, mark each `findings_for_pk.md` entry `FOLDED <date>`. Then tell Claude Code to move the `CLAUDE.md` / `docs/data_map.md` stamps to v1.11.

---

## 1. `index.md`

| Location | Current text (anchor) | Change | Finding |
|---|---|---|---|
| §1.1 instrument table, Ka row (~L60) | "−19.26 dBZ at 1 km `[GAP G1]`" | Per-granule value: −16.32 for 15:13–18:59, −19.26 for 19:00–19:26 on 02-17 | Ka sensitivity |
| §9.5 code repository (~L413) | layout list | Add `coreg.py`, `checks/ka_coverage.py`, `scripts/`, the cache file, and the `nurture-oe` env (mamba, conda-forge, ultraplot). Tasks (a) and (d) done 2026-10-14 | — |
| §11.2 **G1** (~L476) | "sensitivities at 1 km documented (Ka −19.26, …)" | *Narrowed 2026-10-14:* gate spacing 19.986 m (751 gates from −999.3 m; r = 0 at index 50); rays every 0.10044 s (mid-time); Ka sensitivity differs by granule. Still open: range *resolution*, near-field limit, frequencies, pulse scheme | Header facts; Ka sensitivity |
| §11.2 **G21** (~L496) | "averaged vs sampled unstated" | *Narrowed:* HALO `imu_roll` is neither (RMS in turns 6.2° sampled, 5.3° centered 10 s mean; best fit is a 10 s mean centered at t + 6 s, 2.3°). Clocks agree within 1 s. `imu_pitch` − MetNav is +0.29°. Use MetNav for θ_inc until the PI confirms. Optional addition: on 02-17 the radar beam leans ≈47 m/km forward (pitch), so ≈0.56 km ahead at the surface; whether HALO tilts the same way is unknown (`tilt_angle`) | imu_roll |
| §11.2 **G23** (~L498) | "PI says averaged to 10 s × 15 m (double-checking)" | *Narrowed:* the 10 s profile window is empirically [t − 5, t + 5) (±1 s); `lat`/`lon` are instantaneous; spacing irregular, 9–23 s. `1064_dep_cloud_screened_unc` is filled on 10,971 pixels where its product is NaN | Window; `_unc` mask |
| §11.2 **G24** (~L499) | "the SNR definition behind the 1 km sensitivities …" | *Narrowed:* delivered `zm` reaches 5–11 dB below the stated $Z_{\min}(r)$ in every band, with the $r^2$ form intact. The attribute is not the delivered detection threshold, so the SNR definition is now blocking. Fill convention: `_FillValue = NaN` for `zm`/`katten`/`s0` | zm floor; header facts |
| §11.3 **O15** (~L521) | formula "$-19.26 + 20\log_{10}(r/1\,\mathrm{km})$" and "(Claude Code task d)" | Constant becomes the per-granule $Z_{1\,\mathrm{km}}$. Record the condition-3 result: 15:30–16:00, 18,944 ice pixels with W+G; Ka also detects in **4.0 %** (file value −16.32) / 10.9 % (−19.26) / 12.8 % (finite), all in the lowest ≈2 km; 0 % at 6–8 km range. Decision on whether condition 3 is met is Joe's | Ka coverage |
| §11.3 **O21** | non-detection | Note: on 02-17 Ka contributes almost only through this penalty; the threshold definition depends on G24 | Ka coverage; zm floor |
| §11.3 **O31** (~L537) | interim radar $\mathbf{S}_\epsilon$ | The random term's noise floor cannot be the attribute (see G24); the Ka value is per granule | zm floor; Ka sensitivity |
| §11.3 **O32** (~L538) | "profiles stop 100 m above cloud top" | **Contradicted on 02-17:** only 25 % of 1,329 profiles truncate ≥ 85 m above `cloud_top_height_935`; the median lowest valid bin is 495 m *below* it. In-cloud exclusion must use our own mask | h2o truncation |
| §11.3 **O36** | 10 s averaging | Window rule [t − 5, t + 5) (see G23) | Window |
| §11.3 **O29** | correlated rows | Co-registered radar: 20 m gates on 15 m bins ⇒ adjacent rows share gates; decimate | task (d) method |
| *candidate* §11.1 or §4.1 note | — | Vertical datum: HALO `z` and `gps_alt` are MSL; CloudCube `alt` is WGS84-ellipsoidal (20–44 m apart on 02-17). Co-registration references gates to HALO MSL altitude. Settled by data, so possibly a short DECISION ("radar gates referenced to HALO MSL altitude") | Datum |
| Change log (end) | — | v1.11 line summarizing the above | — |

## 2. `appendix_A_nurture_datasets.md`

| Location | Change | Finding |
|---|---|---|
| §A.1 status table, CloudCube row (~L23) | "no file header inspected yet" → headers inspected 2026-10-14 (three granules: 15:13–15:20, 15:21–18:59, 19:00–19:26) | Header facts |
| §A.4 intro (~L125) | Same; add granule list and the 0.10044 s / 19.986 m / 751-gate facts | Header facts |
| §A.4 "Documented sensitivities" (~L129) | Ka is per granule (−16.32 / −19.26); W and G constant. Delivered floor 5–11 dB below $Z_{\min}(r)$ | Ka sensitivity; zm floor |
| §A.4.1 table: `zm`, `katten` dims (~L135–136) | `(range, time)` → **`(time, range)`** | Header facts |
| §A.4.1 `time` row (~L139) | interval length = 0.10044 s | Header facts |
| §A.4.1 nav row (~L140) | `alt` is **above the WGS84 ellipsoid**; nav and attitude are MetNav interpolated (lag 0) | Datum |
| §A.4.1 `look_vector_ecef` (~L142) | `(3, time)` → **`(time, nxyz)`**; unit length except 98 rays (≤4e-4) | Header facts |
| §A.4(i) Ka formula (~L152) and the "−3.7 dBZ at 6 km" sentence | Per-granule constant; add the coverage result | Ka coverage |
| §A.4(iv) random term | Note that the attribute is not the delivered floor | zm floor |
| §A.8(1) | `State_Type` = 2 (MERRA-2) done; imu result (G21) | State; imu_roll |
| §A.8(2) | "file header" item done; remaining: SNR definition, calibration, absorption model, noise subtraction, range resolution, near-field | Header facts |
| §A.8(4) STM co-registration | Mark done 2026-10-14; point to the figure and report; record the method in one paragraph and the Ka result; slide-framing number: two-way gas at 1 km altitude W 1.3 dB vs G 5.9 dB ⇒ ≈4.7 dB of near-surface DWR_W−G is gas | Ka coverage; gas path |
| §A.10.1 / §A.10.5(f) | Window rule [t − 5, t + 5); spacing irregular 9–23 s; `time:long_name` "instant time"; PK counts reproduce within 0.4 % (10,439 / 62,700) | Window |
| §A.10.2 `h2o_mmr_v` row (~L273) | "Extends … to 100 m above cloud top" → contradicted on 02-17 (see O32) | h2o truncation |
| §A.10.2 `/Nav_Data` row (~L278) | imu timing result | imu_roll |
| §A.10.2 (new rows or a note) | 23 variables present but not named (`*_bsc_Sa`, `532_AOT_from_bsc*`, several `_unc`, nav, `DataProducts/Altitude` = legacy copy of `z`); State units (Pressure atm) | Inventory |
| §A.10.3 list | Add item: `1064_dep_cloud_screened_unc` mask exception | `_unc` mask |
| §A.10.4 "exactly the pixels" (~L313) | "five of six; `1064_dep_cloud_screened_unc` is the exception" | `_unc` mask |
| §A.11(b) co-registration (~L404) | Add the datum step (HALO MSL vs CloudCube ellipsoid) and the method actually used. Surface check: W surface peak within −11 m (15:30–16:00) and +7 m (19:00–19:25) of C&S `altSurface` | Datum; task (d) |
| §A.11(c) fill semantics (~L406) | CloudCube convention now known: NaN `_FillValue` for profiles; nav has none (netCDF default) | Header facts |
| Footer (~L430) | "Next: Claude Code task (d)" → done; next items below | — |

## 3. `handoff.md` (replaced at wrap-up)

Remove "Claude Code task (d)" and "Claude Code task (a)" from NEXT SESSION NEEDS. Task (c), the ARTS catalog, stays pending. Suggested new items:

- **Joe, at the STM / with JPL:**
  - Why Ka sensitivity changed at 19:00 (−16.32 → −19.26).
  - The SNR definition behind the 1 km sensitivities (now blocking G24/O21/O31).
  - Range resolution and the near-field limit.
- **Joe, with the HALO team (J. Collins meeting):**
  - Definition and timing of `imu_pitch`/`imu_roll`.
  - The 10 s window convention.
  - The `1064_dep_cloud_screened_unc` mask.
  - The `h2o_mmr_v` truncation rule.
- **Decision for Joe:** does a 4 % (file threshold) Ka coverage on 15:30–16:00 meet O15 condition 3? A second segment, 19:00–19:25 with Ka −19.26 and higher Ka fill, can be run with the same code.
- **NEXT SESSION READS:** §A.4, §A.8(4), §A.11(b); `docs/reports/task_d_coreg.md` in the repo.

## 4. Findings → entries (cross-reference)

| `findings_for_pk.md` entry (2026-10-14) | Type | Primary targets |
|---|---|---|
| Ka sensitivity differs between granules | conflict-with-PK | §1.1, G1, O15, O31, §A.4 |
| zm extends 5–11 dB below Z_min(r) | GAP-update | G24, O21, O31, §A.4(iv) |
| `alt` ellipsoidal vs HALO MSL | fact | §A.4.1, §A.11(b), *candidate* decision |
| CloudCube header facts | GAP-update | G1, G24, §A.1, §A.4, §A.11(c) |
| 10 s window [t − 5, t + 5) | GAP-update | G23, O36, §A.10.1, §A.10.5(f) |
| `imu_roll` timing | GAP-update | G21, O28, §A.10.2 |
| `h2o_mmr_v` truncation | GAP-update | O32, §A.10.2 |
| `1064_dep_cloud_screened_unc` | conflict-with-PK | G22/G23, §A.10.3, §A.10.4 |
| State / inventory | fact | §A.8(1), §A.10.2 |
| Ka coverage result | fact | O15, D12, §A.8(4) |
| Gas path in DWR | fact | §A.8(4), §4.1 |
