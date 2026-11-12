# CLAUDE.md — nurture-oe

Python code for the NURTURE → FALCON optimal-estimation (OE) cloud retrieval framework: data readers, scattering-database tools, the lidar multiple-scattering package `lidar_ms`, and (later) the other forward models and the solver.

**Distilled from Project Knowledge** (claude.ai Project "OE Retrieval Frameworks: NURTURE → FALCON"): `index.md` **v1.14**, `appendix_A_nurture_datasets.md` (2026-11-11) and `appendix_C_lidar_multiple_scattering.md` (2026-11-11). Project Knowledge is authoritative. If this file, `docs/data_map.md` or the code disagree with it, stop and flag the conflict to Joe. Do not resolve it yourself.

Joe is a capable scientist who is still learning OE: explain OE-relevant choices briefly as you make them.

## Read before working

- `docs/data_map.md` — every data source: paths, formats, structure, known metadata traps, open register items. **Required before touching any reader.**
- `config/paths.yaml` — all data locations. Never hard-code a path.
- `docs/pk_snapshot/` — read-only copies of the PK files (`index.md` v1.14, `appendix_A_nurture_datasets.md`, `appendix_B_stm_figures.md`, `appendix_C_lidar_multiple_scattering.md`, `handoff.md`, `team_followups.md`; refreshed 2026-11-11 at the v1.14 fold; appendix B is unchanged since v1.13). **Required before any `lidar_ms` work:** appendix C §C.2–C.4. The originals are on OneDrive: never write there. Ask Joe to refresh the snapshot when PK moves on.

Environment: `mamba env create -f environment.yml` (conda-forge only), then `pip install -e .` in the `nurture-oe` env. Plots use ultraplot.

## Hard rules

1. **Never invent instrument or database specifications.** If a needed fact is missing, expose it as a parameter, tag it `[GAP]` in a comment or docstring, and record it in `docs/findings_for_pk.md`.
2. **Data stores are read-only.** Never write, rename or delete anything under the paths in `config/paths.yaml`. Caches go to `cache_dir`.
3. **Header governs README.** Where a product's file metadata and its README disagree, the reader follows the file. It also emits a `warnings.warn` naming the conflict. Never silently "fix" units, names or attributes.
4. **Readers do not interpret.** They return native grids, native units and native fill values (as `NaN` / masked), with all attributes attached. They do no regridding, no co-registration, no unit conversion beyond what the file declares, and no merging of products. Those steps belong to separate, explicit functions.
5. **`legacy/` is reference only.** Do not import from it or edit it. Port deliberately: each legacy notebook has known issues listed in `docs/data_map.md`. `legacy/multiscatter-1.2.11/` is Hogan's C source (Apache-2.0), the reference for the `lidar_ms` port; its departures from our formulation are listed in appendix C §C.2.9.
6. **No register numbers.** `[DECISION]`, `[GAP]` and `[OPEN]` numbering (D#, G#, O#) is owned by `index.md` §11. Cite existing IDs, but never create new ones. Record new findings in `docs/findings_for_pk.md` (template below); Joe and Fable fold them into Project Knowledge. Appendix C has its own local labels — departures Dep1–Dep11, gates Gate1–Gate4, assumptions Asm1–Asm9, milestones M0–M13, tests ST0–ST4, A1–A9, N, R, J, K, X, P. They are not register numbers: cite them as written and never shorten Dep to D or Gate to G.
7. **Ask before adding a dependency** beyond numpy, scipy, xarray, h5py, netCDF4 and pandas. For `lidar_ms` that means NumPy and SciPy only — no Numba, JAX or GPU `[DECISION D28]`.
8. **`lidar_ms` never imports `nurture_oe`** `[DECISION D29]`. Its tolerances, reference files (`tests/lidar_ms/reference/`) and equations are fixed: if a test cannot be met, stop and record it in `docs/findings_for_pk.md` (appendix C §C.4.0).

## Conventions

- Notation (match in names, docstrings and plots):
  - state $\mathbf{x}$, observations $\mathbf{y}$, forward model $F(\mathbf{x})$, Jacobian $\mathbf{K}$;
  - prior $\mathbf{x}_a$, $\mathbf{S}_a$; observation + forward-model error $\mathbf{S}_\epsilon = \mathbf{O} + \mathbf{M}$;
  - posterior $\hat{\mathbf{S}}$, averaging kernel $\mathbf{A}$.
  
  Size-spanning quantities ($Z_e$, $\beta$, $\alpha$, $q_v$, $N_0'$) enter $\mathbf{x}$ and $\mathbf{y}$ as natural logs.
- Epistemic tags in comments and docs: `[DECISION Dn]` settled, `[GAP Gn]` missing fact, `[OPEN On]` unresolved design question.
- Python with NumPy-style docstrings and minimal comments. Prefer small, runnable functions over frameworks.
- Package layout: one distribution, two import packages. `src/nurture_oe/` — `io/` (instrument readers), `ssp/` (scattering databases), later `forward/lidar_adapter.py`. `src/lidar_ms/` — lidar multiple-scattering model (file tree in appendix C §C.3.1). `tests/` (with `tests/lidar_ms/` and its frozen `reference/`), `notebooks/` (exploration only; nothing imports from it).
- Large files: open lazily (h5py, xarray `chunks=`, memory-mapped or line-indexed access). Never load a whole TAMUice2016 `P*.dat` (~600 MB of ASCII each) more than once. Convert the needed subset to netCDF in `cache_dir`, and record the source file, wavelength(s) and code version in its attributes.
- Tests run against the example flight (`nurture_example` in `config/paths.yaml`). Mark tests that read >100 MB as `@pytest.mark.slow`.

## Data facts that bite (v1.14; data facts unchanged since v1.13)

- **Heights `[DECISION D23]`:** everything is on the HALO MSL datum (HALO `z`, `gps_alt`; `/State` on the same grid). CloudCube gate altitude is $z = \text{alt}_{\text{HALO}} - r\cos\theta$, with θ from `look_vector_ecef` and HALO altitude interpolated to ray times. Never use CloudCube `alt` (WGS84 ellipsoid) for altitude.
- **Ka sensitivity is per granule:** read `radar_sensitivity_ka_dbz` from each file (−16.32 dBZ for 15:13–18:59 on 02-17, −19.26 for 19:00–19:26). It is also the Ka admission threshold for the optional Ka block `[OPEN O15]`, condition 3 met.
- **The documented 1 km sensitivity is not the processor's detection threshold.** Delivered `zm` reaches 5–11 dB below $Z_{\min}(r)$ in every band `[GAP G24]`.
- **HALO 10 s window `[t − 5 s, t + 5 s)`** is a working convention; a ±10 s sweep could not resolve its offset `[GAP G23]`. In 10.5 s steps one sub-profile lies in no window.
- **Attitude:** never take incidence angle from HALO `imu_*`; use MetNav `[GAP G21]`. `imu_roll` has an unexplained timing/averaging offset; don't describe it as a shifted 10 s mean.
- **`h2o_mmr_v`:** the README's 100 m truncation above cloud top does not hold. Exclude in-cloud rows with our own cloud mask `[OPEN O32]`.
- **`1064_dep` rows** take their mask from the product, not from `_unc` `[GAP G22]`.
- **Co-registered radar rows** are censored means over detecting rays (bias ≥ 1.7 dB at detected fraction ≈ 0.5, ≲ 0.2 dB at ≥ 0.9) `[OPEN O21]`, and on the 15 m grid about one gate in three is duplicated `[OPEN O29]`. `*_all` is the unthresholded detecting-ray mean.
- **Ka, 19:00 granule:** the Ka range response changed (≈1 km skirt around strong echoes); mask Ka rows near strong echoes; Ka calibration is per granule `[GAP G24]`.
- **`look_vector_ecef`** is computed from attitude (body-fixed nadir beam assumed), so it can't test mounting `[GAP G27]`.
- **HALO 0.5 s cloud file** (`io/halo_cloud.py`): interim, no metadata. `bsc` is β_p,532 in m⁻¹ sr⁻¹ and `dep` is δ_p (both inferred, unscreened); `bsc1064` is not on the Aerosol/WV scale; `cp` ≡ C&S `CloudMask/Type`. No ATB and no uncertainty `[DECISION D22]`, `[OPEN O36]`.
- **Both beams are airframe-fixed** (HALO `offset_angle` = cos pitch · cos roll), so 0.5 s co-registration uses a body-fixed lidar (= ray-time matching). An empirical ≈1.1 s radar−lidar lag (radar earlier) is unexplained and not applied `[GAP G27]`.
- **The 0.5 s × 15 m co-registered file** is for display and diagnostics; 0.5 s β_p rows are D21-valid only where HSRL-valid.
- **HALO geometry for `lidar_ms` `[GAP G25]`:** `Instrument.from_full_angles(532e-9, 0.8e-3, 1.0e-3)` (divergence 0.8 mrad full angle at 1/e², field of view 1 mrad full angle). There is no F-Lidar preset: never invent one `[GAP G6]`.

## Current work (updated 2026-11-11)

**Commissioned: the `lidar_ms` build** — appendix C §C.4, milestones M0–M13, one milestone per commit, in order `[DECISION D24–D32]`. Inputs are in place: `tests/lidar_ms/reference/` (15 inputs, 50 reference outputs of `multiscatter` 1.2.11, `pvc_port.py`, `schema_sketch.py`, `make_refs.py`) and `legacy/multiscatter-1.2.11/`. Returns to Project Knowledge: `docs/reports/lidar_ms_validation.md`, `docs/reports/lidar_ms_benchmark.md`, and `docs/findings_for_pk.md` entries for `[OPEN O35]` and `[OPEN O40]`–`[OPEN O46]` (appendix C §C.4.7). Not in this build: `src/nurture_oe/forward/lidar_adapter.py` (waits for the S1 LUT schema) and any comparison with HALO data.

Earlier tasks (a), (c), (d), the follow-ups of Appendix A §A.8(5) and the HALO 0.5 s cloud file (reader + 0.5 s co-registration, §A.8(7)) are done and folded (reports in `docs/reports/`). Candidates awaiting approval are listed in Appendix A §A.8(6): dropsonde check of DIAL, radar pitch test, Ka skirt on other flights, rebuild of the 10 s cache, window test with the 0.5 s HALO file, the 1.1 s lag on other segments. When the CF-compliant cloud file arrives, re-run `scripts/check_halo_cloud.py` and `tests/test_halo_cloud.py`.

**Not current:** TAMUoic2019 (`[DECISION D19]`, data not yet downloaded); the radar and passive forward models, the HALO lidar adapter, the solver, $\mathbf{S}_\epsilon$ construction (execution sessions S1–S4 in Claude Chat first).

## `docs/findings_for_pk.md` template

Append entries; never edit past ones.

```
### <YYYY-MM-DD> — <short title>
Type: GAP-update | new-GAP-candidate | OPEN-candidate | conflict-with-PK | fact
Refs: <existing IDs, e.g. G22, O29>  Source: <file / variable / attribute>
Finding: <one or two sentences; numbers with units>
Impact: <which y rows, S_eps, presence masks, or sessions (S1–S10) it affects>
```
