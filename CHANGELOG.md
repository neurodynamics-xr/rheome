# Changelog

## Unreleased

- **The editor's spectral chain, sensors to flow, in MATLAB.** `rheome.eigen.hodge_package` (the modes as vertex potentials, merged by mu, single), `rheome.inverse.hodge_weights` (the kernel carried onto the bridge and fused into the modes: B [2K x C]), `rheome.inverse.hodge_synthesis` (the current from the weights), `rheome.spectral.root` (a chunk's one-sided spectra, the editor's load), `rheome.spectral.column` (the real and analytic column at a time, under a gain) and `rheome.spectral.meyerband` (the band gain). Checked: B rebuilt from a real subject's kernel against the editor's stored B (whole matrix 4e-6, a typical row 3e-8: float32 storage; the worst rows at nearly coincident eigenvalues), the root and columns against the editor's own functions run on the same chunk (`tests/editor_spectral.mts`; equal to float32 rounding), and the synthesis read back to its weights (1e-15). Tests: `tSpectralRoot`, `tGlueHodgeParity/fusedWeightsReadBackTheirSynthesis`, `tNxrStoreParity`. The catalogue is now read in-process through MATLAB's Python (`sqlite3`, read-only), not the command-line tool.
- **The editor's spectral-chain inputs are read from an nxr store.** `rheome.io.read.zarr` reads any rectangular part of an uncompressed zarr v3 array, single-chunk, a regular grid or `sharding_indexed` shards (only the inner chunks a range touches, by the shard index; an absent chunk is the fill value; a compressed array is refused). On top of it, with the editor's own lookups: `nxrrecording` (any span of samples, fs and t0 from the time Line, channel names, bad flags), `nxrkernel` (the kernel, its channels as recording rows, the recording it applies to, and its stored fused weights B with the eigenvalues and glued mesh they were made with) and `nxrhodge` (a stored basis). `tZarrRead` holds the reader to arrays written by zarr-python (synthetic fixtures, `tests/fixtures/zarr/make.py`); `tNxrStoreParity` reads a recording chunk, the kernel and B against zarr-python on a real store and rebuilds B's basis (glue, flips, solve), its eigenvalues against those B records. `nxrsurface` now reads through `rheome.io.read.zarr`.
- **`rheome.io.read.nxrsurface` reads a subject's cortex from an nxr store** (the editor's dataset folder: `dataset.sqlite` and `<subject>.nxr.zarr`): the primary surface (or `Surface=`), its faces and its medial wall (the Desikan-Killiany `unassigned` vertices), by nxr-convert's own lookups. The catalogue is read with the `sqlite3` command-line tool (MATLAB's `sqlite` needs the Database Toolbox) and each array from its chunk file; only uncompressed single-chunk zarr v3 arrays are read, any other layout is refused. `tNxrStoreParity` runs the whole chain (read, glue, flips, solve of 200 modes) against nxr-convert's on every subject of `$RHEOME_NXR_STORE`, read-only; nothing from the store is kept in the repository.
- **The editor's surface preparation and Hodge-Dirac solve ported** (nxr-cortical-flow, via nxr-convert `eigen/`): `rheome.geom.glue` (the hemispheres glued by anatomy: each medial wall cut as a disk, a bridge of rings grown and stitched), `rheome.geom.meshkey` (the editor's mesh fingerprint) and `rheome.eigen.hodge_dirac` (the exact and coexact ladders, one lambda cut). Parity in `tGlueHodgeParity`: on the editor's golden mesh the same faces, kept vertices, bridge points (to 1e-15), mesh key, flips and pencil (K to 1e-12, M to 1e-9), and nxr-convert's solve (lambda to 1e-9, modes to 1e-6); on the editor's `cortex_pial_low` (38 678 glued faces) the glue is bit-identical to nxr-convert's, the flips the same, no edge clamped, and lambda to 1e-12.
- **`hodge_dirac` no longer yields a non-positive mass (request 8e5ae937).** `intrinsic_delaunay` now mollifies the edge lengths as cortical-flow does (`opt.mollify`, 1e-6 x mean edge; 0 on a valid mesh, `S.mollifyEps`), and `hodge_dirac` floors the 1-form mass star1 at `star1Floor` = 1e-6 (`H.star1Clamped` counts the edges): a cocircular edge has cotan weight 0 (or within the flip tolerance -1e-6 below it) and is never flipped, which made M_e <= 0 ('mass not positive'). Unchanged on a generic mesh (clamped 0, mollify 0; the cortical-flow exports still match to 1e-9).
- **Cortical-flow differential operators and eigen-solves ported** (for the nsp eigenmode pipeline): `rheome.operators.spin_dirac` (extrinsic, Hoffmann-Ye, Ye), `vertex_to_face`, `intrinsic_delaunay` (flip queue, mixed-Voronoi mass, signposts), `hodge_dirac` (d+delta on 0/1/2-forms), `connection_laplacian_idt` (vector Laplacian on the intrinsic Delaunay triangulation) and `rheome.eigen.spin_dirac` (Ye quotient solve). Parity with the nxr-cortical-flow exports: matrices to 1e-9 or better on bunny-closed and on a 39 380-face cortex with 10 678 flips; spectra (128 Hodge-Dirac, intrinsic D, three spin Diracs) match the Python reference. Tests: `tSpinDirac`, `tDiracExportParity`, `tHodgeDirac`, `tSpectraParity`, `tConnectionLaplacianIdt`. Not ported: mollification, tufted cover, `intrinsicGauge` (Rheome's `operators.gauge` is a different method).
- **`rheome.snapshot` reads the lab's figure tokens from the library.** Order: `Tokens=`, `$NXR_FIGURE_TOKENS`, `$HOME/workspace/library/styles/tokens/figure-tokens.json` if present, else the viewer's vendored `viewer/figure/figure-tokens.json` (the fallback when neither is set); a named file that is missing is an error with the id `rheome:viz:tokens` (was `rheome:snapshot:tokens`); the viewer is served the file used, and its sha256 and path are written to every `*.marks.json` as `tokens_sha256` and `tokens_path` (new `rheome.viz.tokenspath`).

### Changed behaviour (check results computed before this)

- **`rheome.geom.geodesic` measures surface distance correctly.** It was wrong on cortical meshes. On fsaverage6
  pial it saturated near 46 mm (true max ~240 mm) and was 20-30 % off within 40 mm. On an MS1 cortex
  (omega_sub0002, left) it was ±35 % (p5/p95) at 5-40 mm and 20-28 % short beyond 100 mm. On fsaverage, which has
  zero-area triangles, its Cholesky failed. It now follows geometry-central's robust heat method: intrinsic
  edge lengths, mollified, flipped to the intrinsic Delaunay triangulation, lumped mass, and a scale-invariant
  Poisson pin. On fsaverage5/6/fsaverage pial it agrees with geometry-central to <0.3 % at 5-40 mm and to
  <=1.1 mm everywhere. Both are 1-6 % above the exact edge-flip geodesic. Vertices on another connected
  component are `Inf`. Results that planted along it change: `flow.seedvortex` / `vortexatom` envelopes
  (MS1 G2 composition, rotation, diracangles) and `flow.movingvortex` paths (MS1 G11). In one placement the old
  path had a nominal 126 mm, but its ends were truly 178 mm apart, so it planted 0.141 m/s instead of 0.100.
  `tests/tGeodesic`; `tMovingVortex` no longer encodes the saturated 161 mm hemisphere.

### Added

- **Rheome's visualization package is the Cortical Flow viewer (`rheome.viz`; Diellor, 2026-10-14).** Every cortical
  visual is drawn by the viewer: `rheome.view` (interactive) and `rheome.snapshot` (print) are the primary API, and
  `rheome.viz` gives the viewer what it draws: `rheome.viz.tiles` (a `rheome.geom.tree` level as a labelling: tiles by
  level), `rheome.viz.points` (critical points as point sets, one a type, in the tokens' flow-source / flow-vortex) and
  `rheome.viz.marks` (a render's marks sidecar `<render>.marks.json`, nxr.marks@1, for figure_check E1).
  `rheome.io.write.bundle` writes `'Labellings'` and `'Points'` (the viewer's bundle keys); `rheome.view` takes
  `Labellings=`, `Points=`, `Labelling=` and its `.screenshot` writes the marks sidecar beside the PNG; `rheome.snapshot`
  takes `Labellings=`, `Points=` and keeps each render's marks sidecar. `rheome.show`'s cortex renders are the
  fallback (no WebGPU), tests and quick looks; the figures guide now draws with the viewer. The rasteriser-to-viewer
  capability table and the viewer's gaps are in `rheome.viz` (`help rheome.viz`).
- **`rheome.snapshot` draws selections, geodesic paths and shaped insets, and returns each render's camera**
  (figure system revision 3, Diellor's board comments 24-26).
  - `Selection`: a tile (`Vertices`) or a geodesic disk (`Centre`, `RadiusMM` along the folded surface). It is
    drawn as a blended fill in its axis colour at the tokens' alpha, or with `Dashed` as an extent. Its centre is
    the tokens' marker, one size for every centre.
  - `Paths`: geodesics. `From`/`To` are traced by `rheome.geom.geodesicpath`; `Points` with `How` are a path the
    caller traced. Each is ink at `cortex.paths` width with an arrowhead, and is hidden behind the surface and
    off the shown hemisphere.
  - `Inset.Shape`: `tile`, `disk` (cut by the region's mask) or `square`. With no `Vertices`, the inset takes the
    shot's first selection.
  - `out.cameras`: per render `R`, `s`, `off`, `crop`, `origin_mm` and `mm_per_px`, so a caller can place vector
    marks on the page.
- **`rheome.geom.geodesicpath`**: a geodesic between two vertices, traced across the faces by steepest descent of
  the heat distance. On fsaverage6 its median length is within 0.5 % of the exact geodesic, and the board's
  r3 pair comes out at 216.5 mm (exact 215.9; the old trace was 290.5).

- **`rheome.show` follows figure-system revision 1 (figure tokens v2.0.0).** One colormap per quantity: every
  colourbar names a token quantity (divergence, curl, potential, stream, magnitude, amplitude, speed, power, phase,
  code ...) and takes its map; the retired `flow-signed-field` / `'signed'` kind is gone. `surface`, `gif` and
  `filmstrip` take `'Quantity'` (unnamed signed data take the divergence map, the tokens' alias of the retired kind;
  filmstrip's colourbar text moves to `'Label'`), `resolution` takes `'FieldQuantity'`; ridgestrip shows |U| as
  amplitude, eventsensors' footprint and the scalograms as power, activity's |ω| as amplitude. The colormaps'
  descriptive names (tokens colormapName, never `PuOr_r`) are on each figure for its legend:
  `getappdata(h, 'colormapNames')`. Every figure takes `'Theme'` (`''` = Pergammon, 180 mm; `'neuroimage'` for MS1,
  190 mm) for its page and type family. Colourbars and keys sit in the panels' key bands (`nxr_figure 'Keys'`,
  `nxr_colorbar(..., 'F', F)`, `nxr_key`) instead of inside the axes or in a free cell; specparam's gains are labelled
  at their curves. Events are `nxr_event_trace` (labelled trace by default; eventsensors `'EventMode', 'span'` and
  `'EventLabel'`). Symbols and unit exponents are LaTeX maths (any text holding `$...$` is set by the LaTeX
  interpreter; units through `texunit`), never Unicode glyphs. `tests/tShowTokens` runs every figure in both themes
  against `figure_check.py --theme`, and also refuses the retired kind, colourbars without a token quantity, software
  colormap names and Unicode maths glyphs in `+show`.
- **`rheome.view` draws LAYERS: a field, a form, a representation.** The painted colours and the arrows are both
  fields, of two types (Diellor, 2026-10-14), so `Layers=` names what is drawn: one field, one form of it (a real
  field's value `f`; a vector field's `abs`, `nx`, `div`, `curl`, or the vector `x`, `xt`) and one representation
  (`colour`, `isolines`, `halfMax`; `arrows`), with its look, in draw order. A cell of `{field, form, as}` or a struct
  array with the viewer's layer keys; each is checked before anything is written (a vector field is never painted
  itself). It is sent as the Cortical Flow viewer's `layers` setting (viewer 0.1.0+57f83b55 and later). `Colormap`,
  `CLim` and `Threshold` set the colour layer's look, `ArrowDensity` and `ArrowScale` the arrows layers'; `out.set`
  acts on the layers last sent. The old viewer keys (`scalar`, `vector`, ...) are no longer sent.
- **Hodge potentials in the bundle (`rheome.io.write.bundle` `'Potentials'`, `rheome.view` `Potentials=`).** The
  scalar potential and the stream function of a vector field need a Poisson solve the viewer does not do: they are
  solved with `rheome.differential.helmholtz` and written as real fields `<v>_phi`, `<v>_psi` tagged
  `"derivedFrom": "<v>"` and `"form": "phi" | "psi"` (additive to `rheome.bundle@1`).

### Deprecated

- `rheome.view` `Field` and `Arrows`: rewritten into `Layers` with a warning (`rheome:view:deprecated`); removed in
  the next release.

- **`rheome.show` figures read the lab's figure tokens.** Every figure function (eventsensors, activity, critical,
  filmstrip, vectorstrip, ridgestrip, ridges, frame, jointscalogram, jointspectrum, specparam, resolution) and the
  panel functions (surface, vectors, spectrum2d) take their type sizes, line weights, trace/event colours and colormaps
  from `library/styles/tokens/figure-tokens.json` (through `library/styles/matlab`, `$NXR_STYLES` to move it), not
  literals. The figures are journal-sized (190 mm wide, NeuroImage) with panel letters and headings on the house
  grid (`nxr_figure`). Colourbars carry 'quantity (unit)' with round limits, and the default maps are the house
  maps (flow-signed-field, flow-magnitude, atlas-power). A new `'Export'` option writes PDF (fonts embedded), SVG and
  PNG (`nxr_save`); `'Title'` becomes the figure's name, since the legend carries the title. New options say the
  data's unit (`'Unit'`; filmstrip also `'Quantity'`). `tests/tShowTokens` renders each figure and requires its PDF
  to pass `library/styles/schematics/figure_check.py`, and no RGB literal or hand-picked map may remain in `+show`.

- **Print-ready cortical panels (`rheome.snapshot`), part of the figure system.** Every cortical panel of a figure
  is drawn the same way from a shot list (field, quantity, unit, reference, time or strip, surface, view). The
  Cortical Flow viewer draws the views in headless Chrome (Node driver, the viewer served read-only). The layout
  follows the figure tokens: colormap and limits by quantity (signed: PuOr_r, +-p95, symmetric; magnitude: cividis
  from 0), the colourbar in the viewer's own colours ("quantity (unit re reference)"), the four standard views (L
  lateral, L medial, R medial, R lateral), and strips with t above each frame and the interval stated once. The
  scale bar is geodesic: a mesh-edge path on the folded surface cut at exactly its length, drawn where it lies,
  with the camera read off each render and checked. Output is SVG (editable), PDF (fonts embedded, the views at
  full resolution; MATLAB's own PDF export resamples them to ~150 dpi) and PNG at the journal's dpi. The render
  scale is raised per shot so the views print at >= 500 dpi. `out.legend` gives the legend's sentence. The
  viewer-folder lookup moved from `rheome.view` to `+rheome/private/viewerFolder.m` (shared, unchanged). The
  viewer has no BrBG_r or twilight_shifted yet, so coupling and phase quantities need `Colormap=`.
- **Event boundaries read by their posterior (`rheome.detect.risefall` `Model="posterior"`; default of
  `measure_events` and of the §8 A rule in `groupevents`).** Diagnosis on planted episodes through MEG (truth
  known): the 10 dB failure is not the fit's start, bounds, margin, smoothing or loss, and the joint footprint-tile
  envelope does not help. The planted episode beats with the ongoing alpha of its band
  (|p + n|^2 = p^2 + n^2 + 2pn cos(dphi)); with that cross term removed the same fit passes. The least-squares
  point trades the resulting timing error into tau. The new model profiles the power double sigmoid over (tOn,
  tOff, tauR, tauF), takes the marginal posterior medians of tOn and tOff (tempered by the residuals'
  autocorrelation), and refits tau with the times held. Held-out plants (4 participants, 10 dB): every criterion
  is better than "quadrature"; timing is near one 165 ms tile and tau x0.73-0.98. `measure_eventplants` reads
  "posterior" alongside linear and quadrature, can fit the plant's footprint tiles (`Footprint`), and returns the
  plants' envelopes (third output) for estimator work.

- **The atlas's tiles in MATLAB (MS1 tree unification).** `rheome.geom.surfaceladder` ports nsp's atlas ladder
  (`nxr_convert/atlas/ladder.py`: equal-mass cuts across the geodesic diameter) step for step, floating-point order
  and sort emulation included (`rheome.utils.pysort` = CPython 3.12 `list.sort`). It is identical to ladder.py on
  the ten VR1 PREVENT-AD hemispheres at all 12 levels. `rheome.geom.tree` gains `Method="Ladder"` (re-build) and
  `Codes=` (read the stored codes, the only way to match an atlas store built before today's ladder.py: 3 of the
  10 VR1 hemispheres differ from it from level 7 on). The default stays Fiedler bisection; callers opt in.
- **Ambient-field and critical-point trajectory panels (MS1 1.4, 2.3)**, no instrument, on the sphere and the
  subject's left hemisphere. `rheome.scale.measure_ambient` (analysis `ambient`) plants a tilted field and reads it
  ambient vs tangential: the mean-curvature coupling 2H(J.n) against its closed form (H = 1/R on the sphere, cotan H
  on the cortex, also at > 60 mm). `rheome.scale.measure_cptrack` (analysis `cptrack`) plants a drifting rotor and a
  source that waxes and wanes, tracks them with `rheome.detect.track`, and reads core path, speed, curl and div
  against the closed form -2/s^2 at the core. Both write `<analysis>.csv` and figure data `<analysis>.mat`. Sphere:
  coupling 1.7 %, rotor path 1.4 mm median, speed 0.029 of 0.030 m/s, curl 5.8 %, source on/off within 1 frame
  (tests/tScaleAmbient).
- **fMRI on the same tiles as MEG (MS1 code 5/7).** `rheome.io.read.bold` reads nsp bst-fmri's overlay (the
  resting BOLD on the leadfield's cortex, with each run's fMRIPrep confounds); `rheome.scale.locate` finds it and
  `rheome.scale.importsubject` caches it as `fmri.mat`. `rheome.scale.measure_multimodal` (analysis `multimodal`)
  adds the field `fmri_infraslow`: per run percent signal change, motion/CSF/white-matter confounds and a trend
  regressed out, 0.01-0.1 Hz band-pass, and its power through the same graph-wavelet bank on the same tiles as the
  MEG band powers (%^2). Limitation stated: the BOLD point spread through the chain (fsaverage5 sampling,
  interpolation, haemodynamics) is not measured.
- **Cohort reducers on the hub (MS1 4.4, 5.5/5.6)**, `reduce/cohort.py` (Python 3.9, numpy + scipy; reads nsp's
  promoted `<sub>.tar.gz` in place). `gauge`: the cohort tile map of `measure_atlas` "gaugetensor" per dataset and
  template tile -- median and a distribution-free 95 % CI (order statistics) of anisotropy, normal share and the
  trace-normalised tensor, the median tensor's axis, axis agreement, polar tiles' axis withheld. `multimodal`:
  per subject of `measure_multimodal`, relative band power per tile, the whole-cortex relative power and Spearman rho
  across tiles between each relative band and each non-MEG field; then cohort median and CI per quantity (own-tree
  tiles are not the same place across people, so tiles are never pooled). Descriptive only.
- `rheome.scale.measure_multimodal` writes relative band power (`meg_<band>_rel`, the band's share of the four
  bands' stored sums, exact at every level and member) and includes it in `multimodal_assoc`.

- **Event boundaries and event statistics (MS1 SPEC1 §8)**, generic over an event table (one row per event:
  hemisphere, tiles, first and last time), so the tracker's detections and atlas-cell events are read alike.
  `rheome.detect.eventtable` fits `rheome.detect.risefall` to each event's own area-weighted envelope (tOn, tOff,
  tauR, tauF, build-up / plateau / decay, r2, height, extent). `rheome.detect.eventlocked` gives Morlet
  time-frequency power around anchors through a spatial filter. `rheome.scale.measure_events` (analysis `events`)
  reads the tracker's detections: their table, peak frequency and cycles, onset- and offset-locked TF against each
  event's matched random onsets, inter-event intervals, and the event-train and envelope autocorrelation against the
  mode-shift surrogate and a Poisson train. `rheome.scale.measure_eventplants` (analysis `eventplants`) plants
  episodes of known timing in the participant's own rest and reads the boundaries back. `rheome.scale.groupevents`
  gives per-cohort medians with CIs, split-half ICC, the §8 A acceptance and the group TF maps.
  `risefall` takes `TauMax` and flags fits at the bound (`atBound`). `measure_grouptrack` returns its detections
  (fifth output).

- **`rheome.view`, optional** (needs MATLAB R2026a+ or the Cortical Flow app; Rheome works without it on R2023b+):
  a result on the cortex in the 3D viewer from one call, `rheome.view(surf, result, Name=Value)` or
  `rheome.view(result)` with `result.surface`. It opens the Cortical Flow viewer's `viewer.html` (nxr-cortical-flow
  `editor/dist-viewer`) in a uihtml window, which is reused. Options are the viewer's settings (field, arrows,
  colormap, colour limits, threshold band, surface, hemisphere, camera, arrow density and scale, time, play and
  speed), checked before anything is written. `out.set(...)` changes the open view and `out.screenshot(file)` writes
  a PNG. A setting the viewer refuses is a warning. Without a viewer or WebGPU it draws the chosen field and time with
  `rheome.show`. `rheome.io.write.bundle` writes the result (rheome.bundle@1, with further surfaces and the
  settings). Guide: *How to view a result in 3D*.

- **Rheome MCP server** (`mcp/`): Rheome as a product, used through its public API rather than its source.
  It has tools for the scale-up analyses (`run`, `reduce`), the resting-state pattern atlas (`modes_discover`,
  `modes_reduce`), their tables (`readout`), jobs and versions. Every call records the Rheome version and its
  inputs (provenance), and MATLAB runs one at a time. It is not part of the `.mltbx`. See `mcp/README.md`.

## 1.1.0 — 2026-10-11

### Changed behaviour (check results computed with 1.0.x)

- **`rheome.flow.phasegradient`: `.velocity` now points along the propagation** (#18). In 1.0.x it pointed
  against it; wavefronts move at −ω∇φ/|∇φ|². Directions computed from `.velocity` with 1.0.x are reversed.
- **Vector optical flow retired** (#18): `rheome.dynamics.opticalflow_vector` fails on a normal current (what
  an inverse returns). It stays callable but warns on every call; use `rheome.dynamics.opticalflow_scalar` on
  the envelope.
- `rheome.flow` arrow smoothness is read after Levi-Civita transport, and `bandlimit` no longer returns a NaN
  `finestSigmaMM` (#23).

### Added

- **Ports from the predecessor toolbox** (parts 1–4): gauge flow per sphere patch in the group gauge and its
  group reduction (#15); `rheome.geom.spiralorder`, the pole-to-pole spiral order (#16); bad-channel rules v3,
  a FreeSurfer surface reader and connectome operator tests (#17); the phase-gradient sign fix above (#18).
- **Measures behind the multiscale-geometry paper**, per participant and as group reductions, in
  `rheome.scale` (#2–#12, #20, #24): band-limited cortical vector fields and tracked events on the sensor data;
  graph-wavelet tile envelopes (`rheome.scale.coefficients`); Helmholtz-band and noise floors; catalogue
  detectors against phase-randomised nulls; planted patterns, movers, moving 1/f changes and moving vortices
  read through each participant's own MEG; sleep slow oscillations as a positive control; the alpha
  orientation tensor per shared tile; atlas correspondence and connectome roll-up; MEG, PET and fibres on one
  wavelet bank and one dyadic tile ladder; parcel and tile equivalent-disc diameters and hemisphere areas.
- **Direction readouts and their nulls** (#26–#28, #33–#37): `upstreamshare` (share of tile-to-tile travel
  that runs up a cortical gradient, tiles ≥ 10 cm); a lagged-source null that refits generators with no
  travel through the participant's own head; slow-oscillation direction plants at realistic extents.
  **The resting lagged-source null (`laggedrest`) did not pass its calibration and is not validated** (#37);
  do not read direction claims from it.
- **`rheome.modes`**, a resting-state pattern atlas: large-scale complex modes discovered, back-fitted and kept
  only against nulls of the instrument, anatomy, background and second order, with a calibration on every
  head (#30–#32).

### Fixed

- An unreachable `+data` cache says why (drive not mounted, or link dangling), and the data tests skip with
  that reason (#21).
- MEGIN (Elekta/MEGIN) channels are selected (`'MEG MAG'` / `'MEG GRAD'`) (#29).
- Atlas events no longer crash at tile depth 1 (#14); a 120 s empty room yields its aperiodic null windows,
  and one too short drops that arm instead of aborting (#13); group tracking windows lie inside the clean
  span (#25); a span with no clean Welch segment is skipped, not an error (#36).

### Documentation

- **Help text scrubbed**: public help names no participant, private script or internal note, and a guard
  test (`tests/tPublicHelp`) keeps it so (#22). Earlier releases (1.0.0, 1.0.1) are left as published.
- Help-text gaps closed (H1 lines, defaults, Example blocks) (#19); the function reference is regenerated
  from the 1.1.0 help text.

## 1.0.1 — 2026-10-09

Citation release: Zenodo DOI in the README, `CITATION.cff` and `.zenodo.json`. No code change.

## 1.0.0 — 2026-10-08

First public release.
