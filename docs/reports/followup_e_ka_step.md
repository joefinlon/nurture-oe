# Follow-up (e) — Ka calibration vs sensitivity across 18:59:37 / 19:00:00 (2026-02-17)

*2026-10-20. Appendix A §A.8(5e); `[GAP G24]`, `[OPEN O15]`, `[GAP G1]`. Script: `scripts/followup_e_ka_step.py` → `followup_e_ka_step.json`, figure `figures/followup_e_ka_step.png`.*

**Method.** Native rays and gates within ±1, ±2, ±3 and ±6 min of the boundary. The last ray of the 15:21 granule is at 18:59:37.9 and the first of the 19:00 granule at 19:00:00.1. The aircraft is at 12.06 km over land (`landmask` 3) at incidence 1.88–1.89°.
- Step = median(after) − median(before). 90 % intervals come from a block bootstrap (10 s blocks of rays, 2,000 resamples).
- Display gas correction with `katten`, as in task (d).

## 1. DWR$_{\mathrm{Ka}-W}$ in Rayleigh-regime ice: not testable here

**Selection:**
- Ka, W and G all detected; |DWR$_{W-G}$| ≤ 1 dB; $T$ < −5 °C (HALO `/State`).
- Ka above a *common* threshold on both sides: $Z_{\min}(r)$ with −16.32 dBZ (the less sensitive granule's value) + 3 dB.
- ≥ 500 m or ≥ 1.5 km above the surface echo, with or without $Z_W$ < −5 dBZ.

**Result:**
- **Before the boundary, zero gates pass in every variant and window.** Within 18:42–18:59:37, Ka detects no atmospheric echo at all, only the surface. The cloud is weak (column-max $Z_W$ median −8 dBZ), against a Ka threshold of ≈ +3 dBZ at 10 km range.
- The nearest Ka-detected ice before the boundary (18:36–18:42) lies within ≈1 km of the surface.
- After the boundary, gates ≥ 500 m above the surface pass (407–4,989) but none ≥ 1.5 km. They sit inside the range skirt described next. Their DWR$_{\mathrm{Ka}-W}$ of +14 to +27 dB at $Z_W$ < −5 dBZ is not physical for small ice.
- **|DWR$_{W-G}$| ≤ 1 dB alone does not select Rayleigh-regime ice.** Large aggregates also give W–G near 0 while Ka–W is large (+4 to +15 dB in an 18:36–19:15 trial).

## 2. The Ka range response changed at 19:00

Below the surface peak there is no atmosphere, so the contiguous detected extent below the peak measures each band's range response:

| Extent below the surface peak (m), median | Before | After | Step (90 %) |
|---|---|---|---|
| **Ka** (±6 min) | 360 | 1,299 | **+939 (+919 to +979)** |
| Ka (±1 min) | 360 | 1,199 | +839 (+799 to +859) |
| W (±6 min) | 400 | 360 | −40 (−60 to −20) |
| G (±6 min) | 180 | 160 | −20 |

- Above the peak, Ka goes from 320 m to 660–960 m.
- W and G "above" extents include real snow, so they are not comparable.
- From 19:00 the Ka surface echo carries a skirt of ≈1 km above and ≈1.3 km below the peak. W and G show nothing similar.
- **So at 19:00 the Ka waveform, range processing or receiver changed, not only a gain constant.** One candidate is pulse compression with range sidelobes, which would also lower the noise floor. That is a hypothesis to put to JPL, not an inference this test supports.
- Later in the granule the skirt grows: at 19:20–19:26, Ka extends 1.6 km below the surface peak (W 0.54 km), and Ka "detects" 496,711 gates more than 1.5 km above the surface where W detects 26,851 and G 33,453 (`19:20-19:26 check` in the JSON). W and G are not attenuated there: their surface peaks (34.3 and 16.8 dBZ) match 18:54 (34.8 and 16.7), so W missing those echoes is not the explanation. **Ka "detections" in the 19:00 granule near strong echoes are contaminated.**

## 3. Surface σ0 and surface-peak `zm`

| Span | σ0$_{Ka}$ − σ0$_W$ (gas-corr.) | σ0$_G$ − σ0$_W$ (gas-corr., control) | Peak `zm_ka` − `zm_w` |
|---|---|---|---|
| ±1 min | **−3.30** (−3.65 to −2.69) | +0.67 (+0.40 to +0.87) | −4.14 (−4.51 to −3.50) |
| ±2 min | −2.37 (−2.97 to −1.56) | +0.52 (+0.14 to +0.85) | −3.11 |
| ±3 min | −1.98 (−2.52 to −1.45) | +0.71 (+0.40 to +0.99) | −2.79 |
| ±6 min | −1.62 (−1.95 to −1.32) | +1.52 (+1.20 to +1.87) | −2.27 |

Raw (no gas correction) σ0 steps are within 0.1 dB of these for Ka − W. For G − W they are +0.25 to +0.86 dB.

- **Liquid.** None of the columns within ±6 min is liquid-free (C&S class 0 is present in every column), so liquid attenuation is not excluded. The G − W control moves by +0.5 to +1.5 dB, so the path or the surface also changed across the boundary. With more liquid before, part of the Ka − W step would be attenuation, with the same sign. Sizing that needs the liquid absorption at CloudCube's exact frequencies `[G1]`.
- **Window dependence.** The Ka − W step shrinks as the window widens: within each side, block medians drift by 1–3 dB (figure). The sharpest contrast, −3.3 dB within ±1 min, is the cleanest of these estimates. The ±6 min −1.6 dB mixes in scene change.
- **Surface-peak `zm`.** The peak falls 3.3–3.7 dB at Ka but only 0.3–0.8 dB at W and G. A broader range response spreads the surface return over more gates and lowers the peak with no calibration change at all, so the peak step does not measure calibration.
- **σ0 gas correction.** Whether `s0_*` is already gas-corrected is not documented. Raw σ0$_G$ − σ0$_W$ ≈ +0.1 to +0.9 dB; with $2\int k_{\text{atten}}\,dr$ added it becomes +5.8 to +7.3 dB. This does not affect the Ka − W step (Ka and W gas paths are similar), but matters for any absolute σ0 use `[G5]`, `[O30]`.

## Reading

- **The calibration-vs-sensitivity dichotomy does not fit.** The Ka *range response* changed at 19:00. That is a configuration or processing change, which can change sensitivity (consistent with the 2.94 dB attribute change and the 2.1–2.7 dB floor shift of v1.11) and the effective calibration together.
- **There is a Ka-specific drop in surface σ0 relative to W** of 1.6–3.3 dB, depending on span. Its sign is "Ka lower after", while the noise floor moved ≈3 dB *lower* (more sensitive). The ±1 min value (−3.3 dB, 90 % −3.7 to −2.7) is the size of the triple-frequency signal.
- **Not separable here:** how much of that drop is a calibration offset, liquid attenuation, or σ0 processing that depends on the range response. The Rayleigh-ice DWR test that would separate them has no data before the boundary.
- **For O15:** the Ka calibration term must at least be treated as per granule, and Ka rows in the 19:00 granule within ≈1–1.5 km of strong echoes (surface, heavy snow) are suspect. The "19:00 Ka fill 21 %" and the −19.26 dBZ coverage bracket are likely inflated by the skirt. The 15:30–16:00 coverage result is unaffected: in the 15:21 granule the Ka extent below the surface is 360 m, like W.
- **Questions for JPL:**
  1. What changed in the Ka chain at 19:00: waveform, pulse compression, receiver gain, or processing?
  2. How is `s0` computed (range integration, gas correction)?
  3. Is the 15:13–18:59 Ka calibration the same as the 19:00 one?
