# Vortex signal goldens (steps 5-7; written by nxr; for Rheome and the Cortical Flow Editor)

Regenerate: `python tests/goldens/vortex_signal_goldens.py` (deterministic, seeded, synthetic only). `tests/test_vortex_signal_goldens.py` regenerates in memory,
compares with this folder and asserts the planted answer. About 2.4 MB. Reading, units and file format are those of `tests/goldens/vortex/README.md`
(`manifest.json`: per array `file`, `dtype`, `shape`, `unit`, `note`, `source` = nxr function and line; little-endian, C order; indices 0-based). Conventions a port must reproduce:
`docs/plans/2026-11-02-vortex-port-notes.md`, "Steps 5-7".

**Scenario.** 6 sensors, 300 Hz, 16 s: three planted 10 Hz bursts (half-rise/half-fall 2.5-4.0, 6.5-7.7, 10.0-12.0 s; tau 0.12-0.25 s) times a per-channel gain and phase,
a 0.15 s transient at 13.8 s, white noise sd 2, and frames 3060-3239 (8.5-9.0 s) flagged not good. Band: Morlet centre 10 Hz, FWHM 4 Hz; 360 frames/s = 36 per cycle.

**The null's randomness.** nxr draws ONE `rng.standard_normal((C, m))` (C order) in `lagfree_real`. `null_white` is that draw (`default_rng(777)`, PCG64); a port cannot
reproduce numpy's generator, so it reads `null_white` and applies the deterministic transform (`null_amp_spectrum`, `null_sqrt_cov`, `null_filtered_white`, `null_row_sd`
-> `null_lagfree` -> `null_envelope_matched`). For its own runs a port may use any generator: the null is a Gaussian surrogate, not a reproducible quantity.
