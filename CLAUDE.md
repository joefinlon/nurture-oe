# CLAUDE.md — nurture-oe

Python code for the NURTURE → FALCON optimal-estimation (OE) cloud retrieval framework: data readers, scattering-database tools, and (later) forward models and the solver.

**Distilled from Project Knowledge** (claude.ai Project "OE Retrieval Frameworks: NURTURE → FALCON"): `index.md` **v1.10** and `appendix_A_nurture_datasets.md` (2026-10-09). Project Knowledge is authoritative. If this file, `docs/data_map.md` or the code disagree with it, stop and flag the conflict to Joe. Do not resolve it yourself.

Joe is a capable scientist who is still learning OE: explain OE-relevant choices briefly as you make them.

## Read before working

- `docs/data_map.md` — every data source: paths, formats, structure, known metadata traps, open register items. **Required before touching any reader.**
- `config/paths.yaml` — all data locations. Never hard-code a path.
- `docs/pk_snapshot/` — read-only copies of the PK files (`index-1-10.md`, `appendix_A_nurture_datasets.md`, `handoff.md`; copied 2026-10-14). The originals are on OneDrive: never write there. Ask Joe to refresh the snapshot when PK moves on.

Environment: `mamba env create -f environment.yml` (conda-forge only), then `pip install -e .` in the `nurture-oe` env. Plots use ultraplot.

## Hard rules

1. **Never invent instrument or database specifications.** If a needed fact is missing, expose it as a parameter, tag it `[GAP]` in a comment or docstring, and record it in `docs/findings_for_pk.md`.
2. **Data stores are read-only.** Never write, rename or delete anything under the paths in `config/paths.yaml`. Caches go to `cache_dir`.
3. **Header governs README.** Where a product's file metadata and its README disagree, the reader follows the file. It also emits a `warnings.warn` naming the conflict. Never silently "fix" units, names or attributes.
4. **Readers do not interpret.** They return native grids, native units and native fill values (as `NaN` / masked), with all attributes attached. They do no regridding, no co-registration, no unit conversion beyond what the file declares, and no merging of products. Those steps belong to separate, explicit functions.
5. **`legacy/` is reference only.** Do not import from it or edit it. Port deliberately: each legacy notebook has known issues listed in `docs/data_map.md`.
6. **No register numbers.** `[DECISION]`, `[GAP]` and `[OPEN]` numbering (D#, G#, O#) is owned by `index.md` §11. Cite existing IDs, but never create new ones. Record new findings in `docs/findings_for_pk.md` (template below); Joe and Fable fold them into Project Knowledge.
7. **Ask before adding a dependency** beyond numpy, scipy, xarray, h5py, netCDF4 and pandas.

## Conventions

- Notation (match in names, docstrings and plots):
  - state $\mathbf{x}$, observations $\mathbf{y}$, forward model $F(\mathbf{x})$, Jacobian $\mathbf{K}$;
  - prior $\mathbf{x}_a$, $\mathbf{S}_a$; observation + forward-model error $\mathbf{S}_\epsilon = \mathbf{O} + \mathbf{M}$;
  - posterior $\hat{\mathbf{S}}$, averaging kernel $\mathbf{A}$.
  
  Size-spanning quantities ($Z_e$, $\beta$, $\alpha$, $q_v$, $N_0'$) enter $\mathbf{x}$ and $\mathbf{y}$ as natural logs.
- Epistemic tags in comments and docs: `[DECISION Dn]` settled, `[GAP Gn]` missing fact, `[OPEN On]` unresolved design question.
- Python with NumPy-style docstrings and minimal comments. Prefer small, runnable functions over frameworks.
- Package layout: `src/nurture_oe/io/` (instrument readers), `src/nurture_oe/ssp/` (scattering databases), `tests/`, `notebooks/` (exploration only; nothing imports from it).
- Large files: open lazily (h5py, xarray `chunks=`, memory-mapped or line-indexed access). Never load a whole TAMUice2016 `P*.dat` (~600 MB of ASCII each) more than once. Convert the needed subset to netCDF in `cache_dir`, and record the source file, wavelength(s) and code version in its attributes.
- Tests run against the example flight (`nurture_example` in `config/paths.yaml`). Mark tests that read >100 MB as `@pytest.mark.slow`.

## Current work (updated 2026-10-14)

Tasks (a), (c) and (d) are independent of each other, except that (d) uses the (a) readers. Each task's deliverables are code, tests, and a short report under `docs/reports/`, with every unresolved item added to `docs/findings_for_pk.md`.

**(a) NURTURE readers + data-level checks.** Modules: `io/halo_aerowv.py`, `io/halo_cloudsurface.py`, `io/cloudcube.py`, `io/metnav.py`.

1. Before writing any reader, dump each example file's full header (groups, variables, dims, dtypes, attributes, fill values) to `docs/reports/headers/`. Diff the HALO Aerosol/WV header against the 2026-01-24 dump summarized in `docs/data_map.md`; the variable set may differ between flights.
2. Then run the checks listed per product in `docs/data_map.md` ("Checks for task (a)").

**(c) ARTS SSDB catalog** (groundwork for `[OPEN O27]`). Module: `ssp/arts_catalog.py`. Output: `docs/reports/arts_catalog.md` plus a machine-readable CSV. Per-habit content is listed in `docs/data_map.md` §ARTS. Use the vendor Python interface in the data store (`DataInterfaces/Python`) read-only, or plain netCDF4.

**(d) STM co-registration + Ka coverage check** (appendix §A.8(4), §A.11(b); handoff 2026-10-09). Segment 2026-02-17 15:30–16:00 UTC (CloudCube granule `…S20260217a152158_E20260217a185937…`). Time-match and map CloudCube slant range onto the HALO 15 m altitude grid. Produce the four-panel figure: `532_total_attn_bsc`, `532_dep`, $Z_W$ (gas-corrected only), and $\mathrm{DWR}_{W\text{-}G}$ with the C&S class-0 contour. Ka coverage = fraction of ice pixels with W and G detected where Ka also detects, with $Z_{\min}^{\mathrm{Ka}}(r)$ taken from each granule's attribute. Modules: `coreg.py`, `checks/ka_coverage.py`. Report: `docs/reports/task_d_coreg.md`.

**Not current:** TAMUoic2019 (`[DECISION D19]`, data not yet downloaded); forward models, solver, $\mathbf{S}_\epsilon$ construction.

## `docs/findings_for_pk.md` template

Append entries; never edit past ones.

```
### <YYYY-MM-DD> — <short title>
Type: GAP-update | new-GAP-candidate | OPEN-candidate | conflict-with-PK | fact
Refs: <existing IDs, e.g. G22, O29>  Source: <file / variable / attribute>
Finding: <one or two sentences; numbers with units>
Impact: <which y rows, S_eps, presence masks, or sessions (S1–S10) it affects>
```
