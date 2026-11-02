# Vortex tile goldens (steps 8-9; written by nxr; for Rheome and the Cortical Flow Editor)

Regenerate: `python tests/goldens/vortex_tiles_goldens.py` (deterministic, seeded, synthetic only). `tests/test_vortex_tiles_goldens.py` regenerates in memory, compares
with this folder, replays the shuffles from the stored permutations and asserts the planted answer. About 2.5 MB. File format, units and reading: those of
`tests/goldens/vortex/README.md` (`manifest.json`: per array `file`, `dtype`, `shape`, `unit`, `note`, `source` = nxr function and line; little-endian, C order; 0-based; NaN = missing).
Conventions a port must reproduce: `docs/plans/2026-11-02-vortex-port-notes.md`, "Steps 8-9".

**Scenario.** Two icospheres (642 vertices each, radius 60 mm, centres y = +/-90 mm, 0.6 mm jitter): `V_mm`, `F`, `hemi`, `labels` [9 levels x 1284] = `surface_ladder`. `support` (index%11 != 5) and `normals` are the stand-ins for `Cortex.support/normals`.
`steps_A`, `steps_B`: [6030 x 9] core steps (sense, from xyz, displacement xyz, track, frame), planted eastward on the upper half of each sphere, random elsewhere; run C (random everywhere) is stored only as sums.

| group | arrays |
|---|---|
| step 8 | `labels`, `nearest_support_A`, `cross_sums_A_l*`, `cross_n_A_l*`, `cross_key_A_l6/8`, `cross_vec_A_l6/8` (levels 6, 7, 8) |
| step 9 R | `stage_R_A`, `stage_shuf_A`, `stage_crossings_A` [3 levels x 4 stages], `burst_R_{A,B,C}`, `sums_{A,B,C}_l*`, `n_{A,B,C}_l*`, `stage_perms`, `stage_perm_lengths` |
| place variant | `dc_places`, `dc_result` (R, chance, n), `dc_perms`, `dc_perm_lengths`; input = the first 4000 rows of `steps_A` |
| repeat, tile dot, BH | `repeat_AB_l*`, `repeat_AC_l*`, `tile_dot_AB_l*`, `tile_dot_AC_l*`, `tile_dot_null_l*`, `tile_p_AB_l*`, `tile_q_AB_l*`, `bh_p`, `bh_q` |

**Randomness.** nxr calls only `rng.permutation(n)`; all permutations consumed are stored in call order (concatenated, with lengths) and replayed; a port needs no generator.
