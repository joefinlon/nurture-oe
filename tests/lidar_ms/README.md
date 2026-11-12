# lidar_ms handoff bundle

Produced 2026-11-11 in the design session for Appendix C (lidar multiple scattering); reissued the same day at the `index.md` v1.14 fold. Only this README differs from the first issue (departure label D2 renamed Dep2); every other file is byte-identical.
Unpack `reference/` into `tests/lidar_ms/reference/` of the `nurture-oe` repository. Use is defined in Appendix C §C.4.

## Contents

| Path | What it is |
|---|---|
| `reference/*.in` | 15 inputs in `multiscatter` ASCII format: six HALO-geometry cases (`halo_*`) written by `make_refs.py`, nine copied unchanged from `multiscatter-1.2.11/examples/` |
| `reference/*__<tag>.ref` | 50 outputs of `multiscatter` 1.2.11 (gcc, double precision). The first comment line of each file records the options used |
| `reference/LICENSE`, `reference/NOTICE` | `multiscatter` licence (Apache-2.0) and notice; they must stay with the nine shipped inputs |
| `pvc_port.py` | Small NumPy re-implementation of the PVC algorithm in `compat` and native modes. Matches the C code to 6e-9 on `halo_liquid_30`. Run: `python pvc_port.py <multiscatter executable>` |
| `schema_sketch.py` | Runnable sketch of the `Instrument`, `RangeGrid`, `OpticalProfile` dataclasses |
| `make_refs.py` | Regenerates every `.ref` file. Needs the C executable; link it with `-lm` after `-lmultiscatter` |

## Tags

| Tag | Options |
|---|---|
| `single` | `-algorithms single none` |
| `pvc` | `-algorithms fast none` |
| `pvc_hsrl` | `-algorithms fast none -hsrl` |
| `h06` | `-algorithms original none` |
| `explicit2` | `-algorithms explicit none -explicit 2` |
| `tdts_only` | `-algorithms none lidar` |
| `tdts_only_nolobe` | `-algorithms none radar` |
| `pvc_tdts` | `-algorithms fast lidar` |
| `single_tdts_nolobe` | `-algorithms single radar` |

## Known properties

- `i3rc_lidar_exp5__pvc*.ref` and `i3rc_lidar_exp6__pvc*.ref` contain NaN from the first clear bin below the cloud onward. That is the C code's behaviour, not a transcription error.
- With `-hsrl`, 1.2.11 applies the droplet near-backscatter factor to the air return and prints a warning. The `pvc_hsrl` references therefore show departure Dep2 of Appendix C §C.2.9.
- HALO geometry in the `halo_*` inputs: divergence 0.8 mrad (full angle, 1/e^2), field of view 1 mrad (full angle), instrument at 10.8 km. Cloud properties are test values, not retrievals.
