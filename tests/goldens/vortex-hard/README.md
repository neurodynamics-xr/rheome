# Hard vortex goldens (written by nxr; for Rheome and the Cortical Flow Editor)

Cases the first set (`tests/goldens/vortex/`) does not exercise. Regenerate: `python tests/goldens/vortex_hard_goldens.py` (deterministic, no random numbers, no participant data).
`tests/test_vortex_hard_goldens.py` regenerates in memory, compares with this folder and asserts each scenario's known answer. About 1 MB in total.

**Format** as in the first set: per scenario a folder with `manifest.json` and little-endian row-major `.bin` files (float64 or int32, no header; same reading recipe as `../vortex/README.md`).
`manifest.json` here lists the scenarios; each scenario's own `manifest.json` has `arrays`, `scalars`, `scenario` (what is planted), `scenario_parameters` (vertex indices) and `mesh_ref`.

**The mesh is the first set's, by reference, not copied** (`mesh_ref` = `../../vortex`): `V_mm`, `F`, `area_mm2`, `adjacency_indptr/indices` (642 vertices, radius 60 mm). Read them from there.

**The input is the movie itself.** There are no modes and no flow here: each scenario's `curl_heat10` [V x T] (float32 values stored as float64; per METRE, like the first set) is the heat-smoothed curl handed straight to
`core_stats(Hc, eye(T), frames_used, f0=2, A, area_mm2, V_mm, fs=12, chunk=3000, steps=)`, with the cores put on chosen vertices (one vertex, or a vertex and its neighbour ring) at value +-1.
fs 12, f0 2: a cycle = 6 frames, half = 3. Arrays: `curl_heat10`, `frames_used` (the columns given to core_stats), `patch_*` (`_core_patches`; `patch_frame` is the POSITION in `frames_used`), `steps` (frame = the frame NUMBER),
`dens_ccw/cw`, `life_hist`, `speed_hist`, scalars. `area_used_mm2` appears only where core_stats got a different area than the mesh's.

**Threshold.** With zero background, fewer than 2 % of the entries are nonzero, so the 98th percentile falls in the zero bin and `_core_threshold` returns the histogram floor lo = hi x 1e-6 = 4e-6 (hi = 4 x max): every
nonzero entry is a core. A port must reproduce that too. (If the plateau is 2 % or more of the entries and has ONE value, the threshold is the grid point just above that value and nothing is a core: not used here.)
Only `threshold_edge` has a real percentile.

| scenario | what is planted | nxr's answer |
|---|---|---|
| `competition_single` | frame 0: two +1 cores A, B, 14.1 mm apart; frames 1-5: one core C nearer A than B (both within 15 mm) | C continues A's track; B's track closes after 1 frame; n_tracks 1, 5 steps |
| `competition_greedy` | frame 0: A, B; frames 1-5: c1 < c2 (vertex index), both nearer B, c2 nearer B than c1, c2 within 15 mm of A | greedy in patch order: c1 takes B although c2 is nearer; c2 gets A; 2 tracks, 10 steps |
| `competition_parallel` | two +1 cores per frame, 10 frames, 14-17.5 mm apart, one vertex per frame in parallel (the cross distance is never nearer than a core's own step) | 2 tracks of 10 frames (life 833 ms), 18 steps |
| `tie` | frame 0: +1 cores at mirror vertices a, b (y -> -y); frames 1-5: one core at c on the plane y = 0: distances to a and b EXACTLY equal (bit for bit, checked) | stable argsort: the first inserted track (lower vertex index, listed first) wins; the other closes; 1 track |
| `gap_zero_frame` | a +1 core moving 5 + 6 frames around an empty frame 5; a lone +1 core (frame 2); a -1 core over frames 8-9 (2 frames) | 2 tracks (5 and 6 frames; life 417 and 500 ms); the 2-frame track is 1 step in `steps` but not in life/path/disp; the lone core has no step |
| `gap_missing_frames` | the same, with frames 5, 6 left OUT of `frames_used` (their columns hold a -1 core that must not be seen) | 2 tracks of 5 frames (the jump 4 -> 7 closes them), `frames` = 10, the excluded core absent |
| `recurrence` | 18 frames; a +1 vortex (vertex + ring) stays in place, shifting one edge (8 mm) every 6 frames; a transient +1 vortex far away in frames 0-5 | recur_cycle 0.667 (12 of 18 old cores), jaccard_cycle 0.315, persist_half 0.857, reversal 0 |
| `reversal` | 12 frames; one place flips +1/-1 every 3 frames (ring cores); a steady +1 core far away | reversal 0.5, persist_half 0.5, 5 tracks (four of 3 frames, one of 12) |
| `threshold_edge` | 22 ring cores (154 entries >= 1.5, max 2.0) and one entry e = 0.036820836 in frame 11 on a bin edge; 7549 zeros + e = 98.0 % of 7704 | nxr (float32 log-binning) puts e in bin 2501: threshold 0.0369453; a float64 port gets bin 2500: 0.0368208 (one grid step, 0.34 %, lower); e is below both, no core; 22 / 12 cores per frame |

Chosen vertices are in each `scenario_parameters`. `competition_*`, `gap_*`, `recurrence`, `reversal` and `tie` use the mesh's own `area_mm2`.
