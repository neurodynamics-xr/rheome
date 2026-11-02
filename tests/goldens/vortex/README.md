# Vortex goldens (written by nxr; read by Rheome and the Cortical Flow Editor)

Regenerate: `python tests/goldens/vortex_goldens.py` (deterministic, seeded; no participant data). `tests/test_vortex_goldens.py` regenerates in memory
and fails if nxr's numbers move. Plan: `docs/plans/2026-11-02-vortex-ports.md`. Total about 4.1 MB.

**Reading.** `manifest.json` lists every array (`file`, `dtype` int32 or float64, `shape`, `unit`, `note`), the parameters, the scalars and the nxr commit.
Each `.bin` is raw little-endian, row-major (C order), no header: MATLAB `fread(f,prod(shape),'float64')` then `reshape(fliplr(shape))` and permute
(or `reshape` to the reversed shape and transpose); TypeScript `new Float64Array(buf)` / `new Int32Array(buf)` with index `i*cols + j`. Indices are 0-based.
Positions are in mm here (`V_mm`); nxr keeps its native surfaces in metres (x1000 to mm). **Curl and divergence are per METRE** (nxr's `FaceOperators` run on the metre surface): `curl_vertex`, `div_vertex`, `curl_heat10`, `curl_hat_level*` and the `threshold` are in flow units / m, and a port reproduces them only with positions `V_mm/1000` (with `V_mm` they come out 1000x smaller). `Z_curl`/`Z_div` are mode coordinates (mass in mm^2): flow units / m x mm. Frames are the columns of the "sensor" matrix: C = T = 12.

**The scenario.** Icosphere (3 subdivisions, 642 vertices), radius 60 mm. Two Gaussian rigid-rotation vortices (width 15 mm), one counter-clockwise
(+1, positive curl), one clockwise (-1), each moving 4 mm per frame along a great circle for 12 frames (fs 12 Hz, f0 2 Hz: a cycle = 6 frames, half = 3),
plus 1 % seeded noise. The flow is the "kernel" `K [3V x C]` of nxr, one column per frame, so `Hc @ eye(T)` is the curl movie.

| file | what | nxr source |
|---|---|---|
| `V_mm`, `F`, `normals`, `mass_mm2` (= `area_mm2`) | mesh (F counter-clockwise from outside); vertex normals (outward); a third of each face's area to its vertices | `Cortex.__init__` vortex.py:95-103, `sense.vertex_areas` sense.py:16 |
| `adjacency_indptr/indices` | symmetric CSR vertex adjacency | `support_adjacency` vortex.py:243 |
| `P` [K x V], `lam_per_mm2` [K], `mode_norm_s` | Laplace-Beltrami modes (cotan, mass-orthonormal, K = V = 642), eigenvalues in 1/mm^2, mode norms | built by the script; used as `Cortex.P/lam/s` vortex.py:94 |
| `planted_paths_mm` [2 x T x 3], `planted_sense` | the answer: centres (mm) and senses | script |
| `flow` [V x 3 x T] | the planted current plus noise | script |
| `curl_vertex`, `div_vertex` [V x T] | face-based tangent curl / divergence, area-weighted to vertices (positions cast to float32 first) | `FaceOperators.div_curl` _port/faces.py:46 |
| `Z_curl`, `Z_div`, `Z_norm` [K x T] | the readings in mode coordinates | `Cortex.readings` vortex.py:119, `coords` :105 |
| `mexhat_radii_mm`, `mexhat_gain_curl` [6 x K], `mexhat_curl_r25/r50_modes` [K x T] | Mexican-hat gain t*lam*e^(1-t*lam), t = (r/2)^2, and the hatted curl in modes | `gain`, `hat` vortex.py:125-132 |
| `curl_heat10` [V x T] | curl smoothed by the heat kernel e^(-lam*sigma^2/2), sigma = 10 mm | `core_operator` vortex.py:237 |
| `curl_hat_level6/7/8` [V x T] | curl through the level hat lam*s2*e^(-lam*s2/2), peak at the Weyl centre of the level, normalised to peak 1 | `hat_operator` vortex_cores.py:51 |
| `heat10_threshold_hist_bins` | [4096, lo, hi] of the threshold histogram | `_core_threshold` vortex.py:326 |
| `patch_frame/sense/centre_mm/area_mm2` | every core: frame, sense, weighted centre (mm), area (mm^2); ordered frame then sense (+1 first) | `_core_patches` vortex.py:343 |
| `steps` [S x 9] | per track step: sense, from xyz, displacement xyz, track id, frame | `core_stats(steps=)` vortex.py:361 |
| `dens_ccw/dens_cw` [V], `life_hist`, `speed_hist` | occupancy (fraction of frames above +thr / below -thr); 20 ms and 20 mm/s histograms | `core_stats` vortex.py:361 |
| manifest `scalars` | threshold, cores_per_frame, area_mm2, n_tracks, life_ms, life_cycles, speed_mm_s, net_speed_mm_s, disp_mm, straightness, reversal, persist_half, recur_cycle, jaccard_cycle, jaccard_half, frames | `core_stats` vortex.py:361 |

**Comparing.** `P`, `Z_*` and `mexhat_*_modes` depend on the sign and basis of degenerate modes: a port that builds its own modes compares through the vertex
outputs (`curl_heat10`, `curl_hat_level*`), a port that reads `P` can compare everything. Use rtol 1e-6 on floats (nxr runs the core operator in float32,
so `curl_heat10` and thresholds agree to about 1e-6 relative); patches, tracks and counts must agree exactly.
