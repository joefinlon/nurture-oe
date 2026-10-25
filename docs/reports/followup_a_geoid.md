# Follow-up (a) — geoid check of D23 (2026-02-17)

*2026-10-20. Appendix A §A.8(5a); `[DECISION D23]`, `[GAP G27]`. Script: `scripts/followup_a_geoid.py` → `docs/reports/followup_a_geoid.json`, figure `figures/followup_a_geoid.png`. Geoid module: `src/nurture_oe/geoid.py`.*

## Method

- **Geoid models.** $N$ along track from the PROJ grids `us_nga_egm08_25.tif` (EGM2008, 2.5′) and `us_nga_egm96_15.tif` (EGM96, 15′), from the conda-forge `proj-data` package. Network grid fetching is off.
- **Expected relation.** MSL = ellipsoidal − $N$, so MSL − CloudCube `alt` should equal $-N$.
- **Residual.** $\text{MSL} - \text{alt}_{CC} + N$, with CloudCube `lat`/`lon`/`alt` interpolated to the reference times.
- **Two MSL references:**
  - HALO `Nav_Data/gps_alt` at the 10 s profile times;
  - MetNav `GPS_Altitude_MSL` at 1 Hz.
- **MetNav fields.** MetNav has **no ellipsoidal altitude**: its only altitude fields are `GPS_Altitude_MSL` (m) and `Pressure_Altitude` (ft). The "does CloudCube `alt` match it" test therefore becomes: does `alt` equal MetNav MSL + $N$?

## Results

$-N$ along track is 17–18 m (15:13–15:20), 20–45 m (15:21–18:59) and 14–22 m (19:00–19:26). EGM2008 and EGM96 differ by −0.9 to +0.8 m (1–99 %).

**Residual MSL − `alt` + $N$ (m), 1st / 50th / 99th percentile:**

| Granule | Reference | n | EGM2008 | EGM96 |
|---|---|---|---|---|
| 15:13–15:20 | MetNav 1 Hz | 430 | −0.52 / −0.38 / −0.19 | **−0.03 / −0.02 / 0.00** |
| 15:21–18:59 | MetNav 1 Hz | 13,059 | −0.88 / −0.04 / +0.27 | **−0.03 / 0.00 / +0.05** |
| 19:00–19:26 | MetNav 1 Hz | 1,585 | −0.18 / +0.12 / +0.72 | **−0.08 / −0.01 / +0.06** |
| 15:13–15:20 | HALO 10 s | 43 | +1.4 / +2.0 / +2.8 | +1.6 / +2.3 / +3.0 |
| 15:21–18:59 | HALO 10 s | 1,302 | −7.0 / −1.0 / +4.5 | −6.8 / −0.8 / +4.3 |
| 19:00–19:26 | HALO 10 s | 150 | −3.7 / +1.3 / +4.8 | −4.2 / +1.3 / +4.1 |

The results are the same in level flight (|dz/dt| < 0.5 m s⁻¹).

1. **D23's premise holds.** The 20–44 m HALO − CloudCube difference is the geoid undulation. CloudCube `alt` equals MetNav `GPS_Altitude_MSL` + $N_{\text{EGM96}}$ to within ±8 cm (1–99 %) in all three granules. EGM2008 leaves a 0.2–0.9 m model residual. So CloudCube `alt` is consistent with MetNav MSL converted to the ellipsoid with EGM96, which also implies MetNav MSL is EGM96-referenced. JPL has not stated the conversion; this is inferred from a cm-level match, not documented.
2. **HALO `gps_alt` is not MetNav MSL to better than a few metres.**
   - HALO − MetNav: median 0.0 m over the flight, but 5–95 % −6.1 to +3.8 m (std 3.2 m). It does not shrink with a time lag (RMS is smallest at lag 0) or in level flight.
   - It is mostly a slow drift, not noise: profile-to-profile RMS change is 1.4 m. Examples: 16:15–17:10 median −5.2 m; 15:30–16:00 −0.2 m; 19:00–19:25 +1.2 m.
   - It correlates with $N$ (r = 0.64), which hints that the two altitudes use different geoid realizations. That is a guess; HALO's GPS and geoid are not documented.
   - The 2026-10-14 statement "`gps_alt` equals MetNav MSL (median +0.1 m)" is right about the median and silent about this spread.

## What this means for co-registration

- **D23 is unaffected in principle.** Gate altitudes reference HALO `gps_alt`, the same altitude the HALO `z` grid presumably hangs from. Radar–lidar registration therefore stays internally consistent whatever HALO's geoid is.
- **"HALO MSL" is not the same absolute datum as MetNav MSL.** The two differ by up to about 6 m (0.4 bin) over long stretches. This matters if MetNav altitude, a DEM or MERRA-2 geopotential heights are ever mixed with HALO `z`.
- **Not tested.** Whether HALO `z` is built from `gps_alt` or from another altitude; the HALO GPS and geoid model (ask J. Collins); the antenna lever arm `[G27]`.
