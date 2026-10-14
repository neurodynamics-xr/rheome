# Changelog

## Unreleased

### Added

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
