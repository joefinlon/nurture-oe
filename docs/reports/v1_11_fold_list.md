# Fold list for `index.md` v1.11 / appendix A (from Claude Code, 2026-10-14)

Source entries: `docs/findings_for_pk.md`, all dated 2026-10-14. No register numbers were minted. After folding, mark each entry `FOLDED <date>` and move the `CLAUDE.md` / `docs/data_map.md` stamps to v1.11.

## Substantial (change a register entry or a PK statement)

| # | Finding | Proposed fold target |
|---|---|---|
| 1 | Ka sensitivity is **−16.32** dBZ @ 1 km in the 15:13 and 15:21 granules (−19.26 only at 19:00); the floor shift is real | G1, G24, O15 (the formula's constant must be per granule), §A.4 "Documented sensitivities", §A.4(iv) |
| 2 | Delivered `zm` reaches 5–11 dB below the stated $Z_{\min}(r)$ in every band | G24 (SNR definition now demonstrably needed), O21, O31, §A.4(iv) random term |
| 3 | CloudCube `alt` is ellipsoidal; HALO is MSL (20–44 m) | §A.4.1 `alt` row, §A.11(b) co-registration (datum step), possibly an S2 note |
| 4 | Ka coverage on 15:30–16:00: 4.0 % (file) / 10.9 % (PK) / 12.8 % (finite), all in the lowest about 2 km | O15 condition 3 result; D12 "type/riming limited to gates where Ka detects" in practice |
| 5 | HALO 10 s window = [t − 5, t + 5) (empirical); `lat`/`lon` are instantaneous | G23 (narrow), O36, §A.10.1 / §A.10.5(f) |
| 6 | HALO `imu_roll` timing mismatch with MetNav; pitch offset +0.29° | G21 (narrow: attitude source and definition), O28 |
| 7 | `h2o_mmr_v` not truncated 100 m above cloud top in 75 % of profiles | O32, §A.10.2 `h2o_mmr_v` row |

## Facts (update tables / text only)

- CloudCube header facts: dims order, gate spacing 19.986 m, r = 0 at index 50, 0.10044 s rays, fills, three granules. Updates §A.4 ("no file header inspected yet" → inspected), §A.8(2).
- `State_Type` = 2, State units, variable inventory: closes the §A.8(1) item.
- `1064_dep_cloud_screened_unc` mask exception: §A.10.4 "exactly the pixels" sentence, §A.10.3 list.
- Gas path in DWR_W−G (1.3 vs 5.9 dB two-way at 1 km altitude): §A.8(4) slide framing.
