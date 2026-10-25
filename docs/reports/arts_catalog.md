# Task (c) — ARTS microwave SSDB v1.1 catalog (incl. Ka band)

*2026-10-20. Groundwork for `[OPEN O27]` (index §5.2, §10 S1); `[DECISION D20]`, `[GAP G1]`. Code: `src/nurture_oe/ssp/arts_catalog.py`; driver `scripts/run_arts_catalog.py` (≈6 min, read-only). Outputs: `docs/reports/arts_catalog.csv` (one row per habit, every column below and more) and `docs/reports/arts_catalog_sizes.csv` (one row per size file: shape data, frequencies and temperatures present, flags for the target frequencies).*

## Scope and method

- **Source.** `SSD/TotallyRandom/` only, 35 habit folders, 1,649 size files, read with plain netCDF4. Each size file has one group per (frequency, temperature).
- **What a missing group means.** The size was not computed at that frequency (DDA size-parameter limit). `SSD/` holds no extrapolated values. Copying the highest computed frequency upward is done only in `StandardHabits/`, which is not used here (data_map §2.2). The damaged IconHail / IconCloudIce files are in `StandardHabits/`, not `SSD/`.
- **Target frequencies.**
  - **Ka:** the grid offers 31.3, 31.5 and 35.6 GHz. CloudCube's exact Ka frequency is unknown `[G1]`, so all three are reported.
  - **W:** 94.1 GHz.
  - **G:** 228 and 247.2 GHz bracket CloudCube's G band (exact value unknown `[G1]`, ≈239 GHz), so interpolation needs both.
- **Grids on disk.** 34 frequencies, 1–886.4 GHz. **The ReadMe's 10.65 GHz is not on disk**, matching the per-habit `DataSummary` (`N_FREQ=34`). Temperatures: 190/230/270 K for ice, 230–310 K in 20 K steps for liquid.
  - Wherever a (size, frequency) exists, all temperatures exist. No size has a partial temperature set.
- **Orientation.** All 35 habits are `totally_random` (attribute `orient_type` and folder). `SSD_ARO/AzimuthallyRandom` holds liquid only (Id 35, 37); **no azimuthally oriented ice is on disk**.
- **$m(D)$.** The data are `ShapeData/mass` against `diameter_max` (minimum circumscribing-sphere diameter).
  - Shown: the `DataSummary` α, β ($m = \alpha D_{\max}^\beta$, SI) and our least-squares fit over all sizes.
  - The CSV also has a fit over the upper half of log $D_{\max}$, because the ReadMe says its α, β "are calculated for the upper end of d_max", and the ReadMe-table values.
- **$A(D)$.** `ShapeData/diameter_area_eq_aerodynamical` is populated for every size. The CSV gives $A = \pi D_{a}^2/4$ and a power-law fit $A = \gamma D_{\max}^\sigma$; the table shows σ. Its exact definition (projected area for which orientation?) is not stated in the file. Treat it as provisional.

## Catalog

Sizes in mm. "Largest Dveq" is the largest volume-equivalent diameter with computed data at that frequency, at all temperatures. **Bold** in the G column means G coverage stops below the largest size.

| Id | Habit (ReadMe name) | Class | Sizes | $D_{\max}$ | $D_{\text{veq}}$ | α / β (DataSummary) | β fit (all sizes) | σ of $A(D)$ | Largest Dveq, Ka 35.6 | Largest Dveq, W 94.1 | Largest Dveq with both 228 & 247.2 | TAMUice2016 name match (flag only) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Evans Snow Aggregate | Aggregates/Pristine | 35 | 0.03–11.8 | 0.05–2.51 | 0.2296 / 2.408 | 2.07 | 1.75 | 2.51 | 2.51 | 2.51 |  |
| 2 | 5-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–2.23 | 0.3965 / 2.423 | 2.49 | 1.78 | 2.23 | 2.23 | 2.23 | SBR (bullet rosette) |
| 3 | Sector Snowflake | SingleCrystals/Pristine | 34 | 0.02–10.2 | 0.02–1.31 | 0.0008177 / 1.444 | 1.79 | 1.73 | 1.31 | 1.31 | 1.31 |  |
| 4 | Perpendicular 3-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–2.14 | 0.4271 / 2.464 | 2.53 | 1.77 | 2.14 | 2.14 | 2.14 | SBR (bullet rosette) |
| 5 | Flat 3-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–1.88 | 0.2386 / 2.423 | 2.48 | 1.77 | 1.88 | 1.88 | 1.88 | SBR (bullet rosette) |
| 6 | 6-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–2.37 | 0.4824 / 2.425 | 2.49 | 1.76 | 2.37 | 2.37 | 2.37 | SBR (bullet rosette) |
| 7 | Column type 1 | SingleCrystals/Pristine | 45 | 0.01–10.0 | 0.01–1.82 | 0.03743 / 2.049 | 2.24 | 1.64 | 1.82 | 1.82 | 1.82 | single_column |
| 8 | 8-Column Aggregate | Aggregates/Pristine | 39 | 0.02–9.71 | 0.01–5.00 | 65.45 / 3.000 | 3.00 | 2.00 | 5.00 | 5.00 | **4.36** | 8_columns |
| 9 | Plate type 1 | SingleCrystals/Pristine | 45 | 0.01–10.0 | 0.01–2.60 | 0.7564 / 2.477 | 2.49 | 2.02 | 2.60 | 2.60 | 2.60 | plate |
| 10 | Perpendicular 4-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–2.07 | 0.3184 / 2.423 | 2.49 | 1.75 | 2.07 | 2.07 | 2.07 | SBR (bullet rosette) |
| 11 | Flat 4-Bullet Rosette | SingleCrystals/Pristine | 45 | 0.02–10.0 | 0.01–2.07 | 0.3184 / 2.423 | 2.49 | 1.78 | 2.07 | 2.07 | 2.07 | SBR (bullet rosette) |
| 12 | Block Column | SingleCrystals/Pristine | 35 | 0.01–2.63 | 0.01–2.00 | 210.5 / 3.000 | 3.00 | 1.97 | 2.00 | 2.00 | 2.00 | single_column |
| 13 | Short Column | SingleCrystals/Pristine | 35 | 0.02–3.30 | 0.01–2.00 | 106.5 / 3.000 | 3.00 | 2.00 | 2.00 | 2.00 | 2.00 | single_column |
| 14 | Long Column | SingleCrystals/Pristine | 35 | 0.02–4.83 | 0.01–2.00 | 33.98 / 3.000 | 3.00 | 1.98 | 2.00 | 2.00 | 2.00 | single_column |
| 15 | Thick Plate | SingleCrystals/Pristine | 35 | 0.02–3.25 | 0.01–2.00 | 112.3 / 3.000 | 3.00 | 2.00 | 2.00 | 2.00 | 2.00 | plate |
| 16 | Thin Plate | SingleCrystals/Pristine | 35 | 0.03–5.06 | 0.01–2.00 | 29.66 / 3.000 | 3.00 | 2.05 | 2.00 | 2.00 | 2.00 | plate |
| 17 | Small Column Aggregate | Aggregates/Pristine | 35 | 0.11–3.86 | 0.04–0.74 | 0.1438 / 2.449 | 2.45 | 2.10 | 0.74 | 0.74 | 0.74 |  |
| 18 | Large Column Aggregate | Aggregates/Pristine | 35 | 0.37–20.0 | 0.13–3.02 | 0.2479 / 2.428 | 2.43 | 2.05 | 3.02 | 3.02 | 3.02 |  |
| 19 | Small Plate Aggregate | Aggregates/Pristine | 35 | 0.10–7.05 | 0.05–1.38 | 0.08097 / 2.253 | 2.24 | 1.83 | 1.38 | 1.38 | 1.38 | 5_plates / 10_plates (plate aggregates) |
| 20 | Large Plate Aggregate | Aggregates/Pristine | 34 | 0.35–22.9 | 0.20–4.56 | 0.2152 / 2.262 | 2.26 | 1.85 | 4.56 | 4.56 | 4.56 | 5_plates / 10_plates (plate aggregates) |
| 21 | Small Block Aggregate | Aggregates/Pristine | 35 | 0.10–7.33 | 0.07–1.67 | 0.213 / 2.332 | 2.29 | 1.93 | 1.67 | 1.67 | 1.67 |  |
| 22 | Large Block Aggregate | Aggregates/Pristine | 35 | 0.35–21.9 | 0.25–4.61 | 0.3224 / 2.253 | 2.25 | 1.85 | 4.61 | 4.61 | **4.32** |  |
| 23 | Spherical Graupel | Aggregates/Rimed | 30 | 0.62–9.74 | 0.45–5.29 | 16.05 / 2.716 | 2.72 | 2.02 | 5.29 | 5.29 | **4.53** |  |
| 24 | Ice Sphere (Mie) | SingleCrystals/Pristine | 200 | 0.00–50.0 | 0.00–50.0 | 480 / 3.000 | 3.00 | 2.00 | 50.00 | 50.00 | 50.00 |  |
| 25 | Liquid Sphere (Mie) | SingleCrystals/Pristine | 200 | 0.00–50.0 | 0.00–50.0 | 523.6 / 3.000 | 3.00 | 2.00 | 50.00 | 50.00 | 50.00 |  |
| 26 | Tyynela Dendrite Aggregate | Aggregates/Pristine | 35 | 0.60–20.8 | 0.23–3.33 | 0.1 / 2.248 | 2.25 | 1.80 | 3.33 | 3.33 | 3.33 |  |
| 27 | ICON Cloud Ice | SingleCrystals/Pristine | 45 | 0.01–10.0 | 0.01–2.93 | 1.59 / 2.560 | 2.56 | 2.03 | 2.93 | 2.93 | 2.93 |  |
| 28 | ICON Snow | Aggregates/Pristine | 35 | 0.12–20.0 | 0.09–3.22 | 0.03159 / 1.951 | 1.96 | 1.74 | 3.22 | 3.22 | 3.22 |  |
| 29 | ICON Graupel | Aggregates/Rimed | 35 | 0.17–6.66 | 0.09–5.00 | 396.4 / 3.135 | 3.14 | 2.14 | 5.00 | 5.00 | **3.97** |  |
| 30 | ICON Hail | Aggregates/Pristine | 35 | 0.12–5.35 | 0.09–5.00 | 382.1 / 2.994 | 3.01 | 2.01 | 5.00 | 5.00 | **3.97** |  |
| 31 | GEM Cloud Ice | SingleCrystals/Pristine | 45 | 0.01–3.09 | 0.01–3.00 | 440 / 3.000 | 3.00 | 2.00 | 3.00 | 3.00 | 3.00 |  |
| 32 | GEM Snow | Aggregates/Pristine | 35 | 0.17–10.5 | 0.09–5.00 | 25.07 / 2.863 | 2.82 | 2.16 | 5.00 | 5.00 | **3.97** |  |
| 33 | GEM Graupel | Aggregates/Rimed | 35 | 0.12–6.60 | 0.09–5.00 | 173.6 / 2.965 | 2.95 | 2.03 | 5.00 | 5.00 | **3.97** |  |
| 34 | GEM Hail | Aggregates/Pristine | 35 | 0.12–5.03 | 0.09–5.00 | 528.2 / 3.017 | 3.05 | 2.01 | 5.00 | 5.00 | **3.97** |  |
| 36 | Brown-Francis 6-Column Rosette | SingleCrystals/Pristine | 32 | 0.17–10.0 | 0.14–1.83 | 0.01854 / 1.900 | 1.90 | 1.51 | 1.83 | 1.83 | 1.83 |  |

Ka at 31.3 and 31.5 GHz has exactly the same coverage as 35.6 GHz for every habit (CSV columns `dveq_max_um_at_31.3GHz`, `..._31.5GHz`).

## Findings

1. **Ka and W cover the full size grid of every habit.**
   - The DDA limit bites only above ≈200 GHz.
   - Ka adds no coverage constraint beyond W.
2. **G band is truncated for 8 habits.**

   | Habits | Ids | Largest Dveq (mm), 228 / 247.2 GHz | Of size grid max (mm) |
   |---|---|---|---|
   | ICON graupel, ICON hail, ICON snow, GEM graupel, GEM hail | 29, 30, 32, 33, 34 | 4.46 / 3.97 | 5.00 |
   | 8-Column Aggregate | 8 | 4.36 at both | 5.00 |
   | Large Block Aggregate | 22 | 4.32 | 4.61 |
   | Spherical Graupel | 23 | 4.53 | 5.29 |

   - For these, a G-band forward model interpolating between 228 and 247.2 GHz has no computed data above Dveq ≈ 4.0–4.5 mm. Extrapolation or a PSD truncation would be an $\mathbf{M}$ term `[G1]`.
   - All other habits are complete at 228/247.2 GHz. The size grids of the snow aggregates end at Dveq 2.5–5.0 mm anyway, so for large snow the PSD tail beyond the grid is the bigger issue for all bands.
3. **$m(D)$: the DataSummary α, β are not a global fit.**
   - Where β < 3 for crystals, our all-size fit differs strongly. Examples: Evans Snow β 2.07 (fit) vs 2.41 (summary), fit RMS 0.24 in log10 $m$; Sector Snowflake 1.79 vs 1.44; Column type 1 2.24 vs 2.05; bullet rosettes ≈2.49 vs 2.42.
   - The upper-half fit comes closer (e.g. Sector Snowflake 1.44, Column 2.03), consistent with the ReadMe note.
   - **So the tabulated α, β describe large particles, and $m(D)$ for small particles departs from them.** Use the per-size masses, not the power law, when building the $F(\mathbf{x})$ tables.
   - Habits with β = 3 exactly (plates, columns of fixed aspect ratio, 8-column aggregate, spheres, GEM cloud ice) are constant-density shapes.
4. **ReadMe vs disk (the file governs):**
   - **Liquid Sphere α:** DataSummary 523.6 (= $\rho_w\pi/6$) vs ReadMe 480, which is the ice-sphere value. A ReadMe transcription error, most likely.
   - **Spherical Graupel α/β:** 16.05 / 2.716 (DataSummary and our fit) vs ReadMe 13 / 2.69.
   - **Sector Snowflake size range:** $D_{\max}$ to 10.24 mm and Dveq to 1.31 mm on disk vs ReadMe 12.0 / 1.415 mm (34 sizes in both).
   - **Short Column:** 35 sizes on disk vs ReadMe 34.
   - **Large Block Aggregate α/β:** 0.322 / 2.253 vs ReadMe 0.35 / 2.27.
   - **Rounding only:** other ReadMe α differences up to ≈20 % come from rounding to 1–2 significant figures (e.g. Flat 3-Bullet 0.239 vs "0.2").
   - **Status:** `DataSummary` `STATUS=WORKING` for 33 habits and `COMPLETED` only for the two Mie spheres. What "working" implies is not documented.
5. **Name matches with TAMUice2016 (flag only; geometric equivalence unverified, `[O27]`):**
   - plate ↔ Plate type 1, Thin Plate, Thick Plate;
   - single_column ↔ Column type 1, Block/Short/Long Column;
   - SBR (bullet rosette) ↔ 5- and 6-Bullet Rosettes and the 3-/4-bullet variants;
   - 8_columns ↔ 8-Column Aggregate;
   - 5_plates / 10_plates ↔ Small/Large Plate Aggregate.

   No ARTS counterpart for TAMU hollow_column, droxtal or HBR (hollow bullet rosette).

## Not done here

- Habit choice and the optical ↔ microwave consistency comparison (`[O27]`, S1).
- Per-size values of $\sigma_{bsc}$, $\sigma_{ext}$ or $\omega_0$ at the target frequencies. The data are there (`SingleScatteringData`), but extraction belongs with the S1/S3 table builder.
- `StandardHabits/` coverage (derived, extended data).
