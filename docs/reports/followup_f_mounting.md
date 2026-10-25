# Follow-up (f) — CloudCube antenna mounting: body-fixed or stabilized? (2026-02-17)

*2026-10-20. Appendix A §A.8(5f); `[GAP G27]`, `[DECISION D23]`. Script: `scripts/followup_f_mounting.py` → `followup_f_mounting.json`, figure `figures/followup_f_mounting.png`.*

**Data.** All three granules, every 10th ray (≈1 Hz): 15,011 samples, 447 of them in turns (|roll| > 5°, up to 37.8°). Attitude is CloudCube's `heading`/`pitch`/`roll`, which is MetNav interpolated at lag 0 (task (a)).

## 1. The commissioned test is circular

- **Fit.** Off-nadir θ from `look_vector_ecef` against $\theta_b = \arccos(\cos\theta_{\text{pitch}}\cos\phi_{\text{roll}})$ gives θ = 0.000° + 1.0000 θ_b. The residual RMS is 0.0008° overall and 0.0045° in turns.
- **Body frame.** Rotated into the body frame (x forward, y right, z down), the look vector is **exactly the body z-axis**: fore/aft and lateral tilt are 0.000° (1–99 %) at every roll angle.
- **What that means.** `look_vector_ecef` is computed from the navigation attitude, assuming an exactly nadir, body-fixed beam with no boresight offset. It is not an independent measurement of where the antenna points. The comparison reproduces JPL's processing assumption, so it cannot tell body-fixed from stabilized, and it cannot measure a mounting offset.
- **Other fields.** `incidence_angle` agrees with this θ to −0.19…+0.16° (1–99 %), so it is presumably derived the same way.

## 2. Independent test: surface-echo range

The W surface peak (> 10 dBZ, beyond 3 km range) gives the measured slant range $r_{\text{sfc}}$. The height above ground is $h = $ HALO `gps_alt` − `UserInput/DEM_altitude` (at nadir, > 5 km).
- Body-fixed predicts $r_{\text{sfc}} = h/\cos\theta_b$.
- Stabilized (vertical) predicts $r_{\text{sfc}} \approx h$.
- The effective angle is $\theta_{\text{eff}} = \arccos(h/r_{\text{sfc}})$.

| | n | θ_eff − θ_b (°), 5 / 50 / 95 % | $r_{\text{sfc}}/h$, 5 / 50 / 95 % | $1/\cos\theta_b$, 5 / 50 / 95 % |
|---|---|---|---|---|
| Turns (|roll| > 5°) | 427 | −9.9 / **−0.9** / +5.6 | 0.999 / 1.030 / 1.158 | 1.006 / 1.035 / 1.142 |
| Level (|roll| < 1°) | 14,265 | −1.8 / +1.5 / +2.7 | — | — |

- In turns, θ_eff follows θ_b (correlation 0.87; median θ_eff 13.9° vs θ_b 14.8°), and $r_{\text{sfc}}/h$ matches $1/\cos\theta_b$. A stabilized antenna would give θ_eff ≈ 0° and $r_{\text{sfc}}/h$ ≈ 1 throughout.
- **The data favour a body-fixed antenna.**
- The ±10° scatter in turns comes from using the DEM at nadir for a footprint up to ≈7 km to the side, plus terrain within the beam. A few rays give θ_eff = 0 (surface higher than the nadir DEM).

## Limits

- **No offset resolved.** This test separates body-fixed from stabilized; it does not measure a boresight offset. The level-flight median of +1.5° and the 5–95 % spread of 4.5° are dominated by DEM and terrain error (≈ 20–40 m in $h$ at 12 km gives degree-level θ_eff error at small angles). A boresight or lever-arm offset of a few tenths of a degree `[G27]` is below what this can see.
- **Geolocation consequence.** Because the L2 geolocation assumes zero offset, any real boresight offset is not in `look_vector_ecef`. Gate altitudes would carry an error $\approx r\,\sin\theta\,\delta\theta$: about 0.9 m at 12 km range for δθ = 0.1° at θ = 2.4°, negligible against 15 m bins. The along-track ground offset is ≈ $r\,\delta\theta$, about 21 m per 0.1°.
- **Lidar side.** The HALO lidar pointing (`tilt_angle`, `[G21]`) remains unknown, so the radar–lidar relative boresight is still open.
