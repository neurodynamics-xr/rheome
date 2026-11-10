# Vortex pair goldens (step 10a; written by nxr; for Rheome and the Cortical Flow Editor)

Regenerate: `python tests/goldens/vortex_pair_goldens.py` (deterministic, seeded, synthetic only). `tests/test_vortex_pair_goldens.py` regenerates in memory, compares with this folder, replays the stored Wishart draws and asserts the planted
answers. About 2.1 MB. Format as `tests/goldens/vortex/README.md` (`manifest.json`: per array `file`, `dtype`, `shape`, `unit`, `note`, `source` = nxr function and line; little-endian, C order; **complex values are two arrays, `_real` and `_imag`**).
Conventions: `docs/plans/2026-11-02-vortex-port-notes.md`, "Step 10a".

**Scenarios** (prefix `rot_`, `cw_`, `ctl_`): 8 sensors, 300 Hz, 20 s, a 10 Hz rhythm (Morlet 10 Hz / FWHM 4 Hz, analytic signal at 100 Hz by `V.analytic`) with slow phase (sd 1 rad) and amplitude (sd 40 %) drift, white noise sd 3.
`rot`: p cos + q sin (q a quarter cycle AFTER p); `cw`: p cos - q sin; `ctl`: one standing pattern (0.6 p + 0.8 q) cos, zero lag. p, q orthonormal sensor patterns. Samples 6-7 s excluded (`ok`); held-out blocks = 5 s even / odd.
Two operators: `operator_I` (identity, 8 x 8) and `operator_L` (a 6 x 8 Gaussian matrix); arrays are tagged `_I_` / `_L_`.

**Per scenario and operator.** `S_all/_even/_odd` (cov over ok / even / odd samples), `fa, fb, u` (the pair fitted on `S_even`; `_raw` = numpy's gauge, unsuffixed = canonical gauge), `proj`, `form` (gauge-free), `share` (on `S_all`), `s3` (on `S_odd`),
nulls (`null_S_all/_odd`, `share_null`, `s3_null`, `p_share`, `p_s3`), and `share_pop` (nxr's `share` on the planted population covariance). The Wishart draws consumed are `draws_all` / `draws_odd` (replay order in the manifest note).
