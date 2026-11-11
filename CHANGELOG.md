# Changelog

## Unreleased

- **Filter kernels: units and definitions made unambiguous (Diellor, 2026-11-10: "Fix the ambiguities."; the same fix is in nxr-neuro's filterbanks).**
  - **Frequency is hertz, signed, everywhere.** `jtv.bank` and `jointfilterbank` pass `f` in Hz (not `omega = 2*pi*f`) to every time factor and joint kernel; `jointfilterbank.band(f1, f2)`'s handle takes Hz; `jointfilterbank.speedkernel(c, d)` takes `f` in Hz and forms `omega = 2*pi*f` itself, so `c` stays a speed in length/s (`d` is the angular width, rad/s). `forward.simulate` no longer hands `2*pi*f` to `resonator` / `gabor` / `travwave` (which read hertz). `jointfilterbank.Omega` stays as a read-only convenience. Every kernel's help states its units.
  - **`filters.signedfrequencies(f, half)`** (new) replaces the private `i_signed` copies in `gabor`, `travwave`, `resonator`, `stmatern`, `bandpass`: `'full'` (default) reads an FFT grid 0 ... fs-df as signed (as before), `'positive'` reads a one-sided grid 0 ... fs/2 and leaves it unchanged (the old rule made its upper half negative: a 20 Hz bin read -0.83 Hz). The five kernels take the half as a last argument.
  - **`wave`, `dampedwave`, `kleingordon` are the continuum forms.** `wave = cos(t*c*sqrt(lambda))`, `dampedwave = exp(-beta*|t|)*cos(t*c*sqrt(lambda))`, `kleingordon = cos(t*sqrt(c^2*lambda + m^2))` (`m` only through omega^2, a rest frequency in rad/s), in physical units, not normalised by lambda_max; on the unit sphere omega_l = c*sqrt(l(l+1))/R. The Chebyshev-angle form (`theta = acos(1 - 2*lambda/lmax)`, a discrete-time propagator whose frequency is not `c*sqrt(lambda)`) is the named variant `'Variant','discrete','LambdaMax',lmax`, never the default, numbers unchanged. The old positional `lmax` argument is refused with a message saying so. `diffusion(lambda, t, D) = exp(-D*lambda*|t|)` and `schrodinger(lambda, t, hbar) = exp(-1i*lambda*t/hbar)` lose their lambda_max normalisation (their time unit was relative to the largest eigenvalue).
  - **`travwave` has a sense.** New argument `sense` (+1 | -1 | 0): the ridge at `+/- c*sqrt(lambda)/(2*pi)` on the signed axis separates a wave travelling one way from the other on a complex (two-sided) coefficient; 0 (default) is `|f|`, the previous behaviour.
  - **`scatter`: scale is the centre wavenumber, the low-pass is not paired.** A member's scale is `1/frame.Centers`; the low-pass member (the scaling function of a mexhat bank, the first itersine window; not a heat bank, whose members are all low-pass) takes no part in the coarser-second pairing, so a mexhat bank of M members computes (M-1)(M-2)/2 pairs instead of (M-1)(M-2)/2 + (M-1), and an itersine bank, which had no second order, gets its band pairs. `S.Sigma` is NaN at the low-pass (was NaN as well) and 1/centre wavenumber elsewhere (was sqrt(2t)).
  - **`timefilterbank`: gain has one meaning, the peak.** New `'Gain'` (any positive number, default 2, recorded): every member's peak in-band gain is `Gain`, as in the Morse, Morlet and Bump families, so a real tone of amplitude A inside a passband gives `|coef| = (Gain/2)*A` (the default now gives `|coef| = A`; it was `0.71*A`). The frame is still tight at any gain: `sum |g|^2 = Gain^2`, `framebounds` reports `Gain^2/2` (2 by default, was 1), and for a real signal `sum ||w||^2 = (Gain^2/2) sum x^2`. `ingest.bank` builds its bank with `Gain = sqrt(2)` so that the stored tile energies stay Parseval-direct (unchanged).
  - **(lambda, t) kernels are lag kernels.** `diffusion`, `wave`, `dampedwave`, `kleingordon`, `schrodinger` act on a time lag and cannot multiply a joint (lambda, f) spectrum; the help of each says so and `jtv.bank` (type `'ts'`) refuses them with that reason.

- **One mesh and mass: the intrinsic, mollified Delaunay triangulation with the mixed-Voronoi mass is the default of every operator (Diellor, 2026-11-10; the Editor's `intrinsicDelaunay.ts` is the reference).** Galerkin / barycentric masses and the embedded mesh stay as named options, never silent defaults.
  - **One Delaunay.** `geom.measure.heatDistance` no longer has its own flip code: it builds on `operators.intrinsic_delaunay` (same mollify, FIFO flips at 1e-6 instead of LIFO at 1e-12, the mixed-Voronoi mass of the mollified input instead of the barycentric A/3 of the flipped mesh). Its solver records `.mesh`, `.mollify`, `.flipEps`; `mollify=` and `flipEps=` are options. The distances change by the mesh and mass difference (see the PR).
  - **`operators.laplace_beltrami(V, F, mesh=, mass=, mollify=, flipEps=)`**: default `mesh='intrinsic'` (`intrinsic_delaunay`'s `S.L`) and `mass='mixedVoronoi'` (`S.mass`); `mesh='embedded', mass='galerkin'` is the previous operator exactly. A third output records `mesh`, `mass`, `mollify`, `flipEps`, `mollifyEps`, `flips`. The positional mass type is gone (name-value only).
  - **`operators.mass`**: default `'mixedVoronoi'` (was `'galerkin'`); `'lumped'` is renamed `'barycentric'` (A/3), the name "lumped" meaning the mixed-Voronoi mass in nxr and the Editor. No alias.
  - **`operators.connection_laplacian(V, F, mesh=, mass=, mollify=, flipEps=)`**: one operator with a mesh option. Default `'intrinsic'` is the construction of `connection_laplacian_idt` (now with `mollify`/`flipEps` and a fourth output, the edge list) returned as the complex `A` and `Rt` on the intrinsic edges; `'embedded'` is the previous operator. The face-walking detectors and flow code (`detect.*`, `flow.*`, `operators.gauge`, `trivial_connection`, the imports' stored bases) pass `mesh='embedded'` explicitly. The result records the mesh, mass type, constants and flips.
  - **`operators.dirac_frame` / `eigen.dirac_frame`**: only the mass default changed (`mass='mixedVoronoi'`, was Galerkin; `mass=` and `mollify=` options, recorded as `MassType`). The relative-Dirac mesh, scaling, normals and domain are an open decision; D_int^2 and E stay on the embedded edge vectors and Gauss map.
  - **Named constants.** `intrinsic_delaunay`'s `eps` is renamed `flipEps` (no alias); the struct records `S.mollify` and `S.flipEps` beside `S.mollifyEps` and `S.flips`. Every operator that builds the intrinsic mesh takes `mollify=` and `flipEps=`.
  - Callers that need the previous operator pass `mesh='embedded', mass='galerkin'` (or `mass='barycentric'`) explicitly: `import.bases/cortex/sphere/surface`, `source.dirac`, `dynamics.opticalflow_*`, `differential.poisson`, `operators.lb_connectome`, `scale.measure_ambient`, the demos, and the fixtures and tests that compare against Brainstorm. New test `tOneDelaunay`.

- **`operators.hodge_dirac`: no star1 floor by default, as nxr (Diellor, 2026-11-10: "nxr has the right floor. Rheome should agree with nxr as a default but with option to change").** `star1Floor` now defaults to 0 (off): star1 is the raw intrinsic cotan edge weight, the numbers nxr and the Cortical Flow Editor build. A positive value (`star1Floor=1e-6` was the old default) floors star1 = max(star1, floor) and counts the edges floored in `H.star1Clamped`; it guards near-cocircular edges whose cotan weight is <= 0, where the mass is not positive and star1^-1 blows up. Tests: `tHodgeDirac` (`cocircularEdgeKeepsTheMassPositive` passes the floor explicitly, `defaultStar1IsTheRawCotanWeight`).

- **Eigen solves and stale help corrected (audit of operators and eigenmodes, 2026-11-10).**
  - **No direct `smallestabs` solve on a singular pencil.** `operators.gauge` (`Method="smoothest"`) and `geom.tree` (the Fiedler vector of each patch) now solve through `rheome.eigen.smallest` (Brainstorm's `bst_eigs_smallest` with the negative-shift ladder), as `eigen.modes` documents. On a pencil with a doubled null mode (two disconnected spheres, 4 modes) the direct `smallestabs` solve had a worst eigenpair residual of 0.21, `eigen.smallest` 8.9e-12; on a well-conditioned pencil both agree (gauge: same field and eigenvalue; tree: the root split is the sign partition of the dense Fiedler vector, the constant mode skipped as column 1 of the ascending result). The confusing comment on `gauge` ("eigs with ONE output returns the EIGENVALUE") is replaced. Test: `tEigenRouting`.
  - **Help that contradicted the code, now true.** `operators.hodge_dirac`: star0 is the mixed-Voronoi mass of the mollified input mesh, not a "lumped mass". `operators.intrinsic_delaunay`: mollification and signposts ARE implemented (the header said "not ported"; the entry below, "Not ported: mollification ...", predates them and stands as history). `operators.spin_dirac` cites `vertex_to_face`. `operators.registry` and `fieldtype.registry` no longer say there is no d0, d1 or Hodge star: `operators.hodge_dirac` builds them (they are not yet rows of the registry; the naming pass will add them). `operators/Contents.m` lists the ten operators it omitted.
  - **`eigen.dirac_frame` states its eigenvalue convention.** The default `Lambda` is a normalised, unitless eigenvalue (spectrum in [0,1]); `normalize=false` gives the physical lambda in 1/m^2 (second order; lambda = mu^2 for the first-order mu in 1/m). Nothing renamed.

- **Rheome matches nxr 1.5 (Diellor, 2026-11-10: "bring Rheome up to date so it matches nxr"; one pipeline, one mathematics, one set of names).**
  - **`rheome.register`, nxr's `nxr.register`, and `rheome.eigen` is the modes only.** Moved, no forwarding names (they were merged the same day): `rheome.eigen.transfer` -> `rheome.register.mode_registration` (nxr's: T = A_template R P'_subject, now taking `rheome.register.solve` structs and a `rheome.register.surfaces` struct as nxr's `Solve` and `Surfaces`), `rheome.eigen.transfercheck` -> `check_mode_registration`, `rheome.geom.sphereregister` -> `registration` (options `SrcCortex`, `SrcColumns`, `Sampling`, `Neighbours`, `Power`), `rheome.geom.hemispheres` -> `hemispheres`. New, all as nxr's: `shepard_weights`, `glued_registration`, `register_in_modes`, `to_glued`, `to_input`, `register_values`, `nearest_vertices`. `nxr.register.tangent` and `interpolation` are not ported. The word is "register" (the Cortical Flow Editor's "Register to template"), not carry or transfer, in names and help; `rheome.inverse.hodge_register` keeps its name. Help text now cites nxr, not the retired nxr-convert.
  - **`Sampling="shepard"`**: Brainstorm's Shepard weights (`bst_shepards`: the k nearest on the sphere, d2 with 0 -> eps, R2 the k-th, w = ((R2 - d2)/(R2 d2))^p over the first k - 1, normalised, a zero row -> nearest, k = 1 -> nearest) with `Neighbours` (8) and `Power` (2), equal to nxr's to 1e-12 on a synthetic sphere (`tRegister/samplingIsNxrs`, with barycentric and nearest).
  - **nxr datasets are read as nxr reads them.** The catalogue must be schema 50 (`PRAGMA user_version`); a 49 is refused with "run `nxr migrate --dataset <dataset folder>`", any other with "re-convert it from its source". `session` is gone (schema 50): `nxrrecording` and `nxrkernel` return `.Acquisition` (the acquisition row: id, label, `ses`, `task`, `run`, `acq`, `proc`, `split`), and a kernel's recording is the one of its acquisition. New `rheome.io.read.nxrregistration`: the stored `hodgeDirac_to_<template>` (kind 'registration'): T, the sampling record and both packages' record (K, mu, kinds, mesh, sha), followed by its rows as the Editor's `hodgeStore` reads it.
  - **`rheome.io.read.zarr` reads compressed arrays**: blosc (nxr compresses arrays over 4 MiB: blosc/zstd, level 3, byte shuffle), also zstd, gzip and a trailing crc32c, whole or as the inner chunks of shards. Only a chunk's decode is numcodecs' (MATLAB's in-process Python); the sharding, regions and placing stay MATLAB. Without numcodecs the error says what to install (MATLAB R2023b takes Python 3.9-3.11: a venv with numcodecs, then `pyenv(Version=...)` before the first Python call). A transposed array is still refused.
  - **Parity tests call nxr** (`$NXR`, the nxr checkout, with `uv`; skipped when absent): `tGlueHodgeParity` (the Editor's goldens from nxr's tests), `tResample`, `tNxrStoreParity` (still needs `$RHEOME_NXR_STORE`), `tests/nxr_convert_chain.py` -> `tests/nxr_chain.py`. New tests: `tRegister` (renamed from `tTransfer`, plus Shepard and nxr parity), `tNxrReaders` (a synthetic schema-50 dataset built by nxr), `tZarrRead` (blosc, zstd, gzip fixtures from zarr-python).

- **Subjects carried onto a template (fsaverage), as the editor's group view does.** `rheome.eigen.transfer` (the transfer matrix T: a template vertex pulled from the subject's cortex through the registered spheres, `rheome.geom.sphereregister`, a hemisphere at a time by `rheome.geom.hemispheres`; the template's bridge from its cut walls; analysed on the template's modes, block-diagonal by kind), `rheome.eigen.transfercheck` (T on spherical harmonics l = 1..3) and `rheome.inverse.hodge_register` (a subject's fused weights onto the template, flow modes as flow fields). `nxrsurface` returns the registered sphere. Checked against the editor's stored registration (nxr-convert registerkernel): the template's basis rebuilt (eigenvalues 1e-10) and the registered weights rebuilt from the subject's kernel to float32 rounding (4e-8). Tests: `tTransfer` (synthetic), `tNxrStoreParity`.
- **`rheome.filters.resample`: a general resampler** (any rate pair whose ratio is a small fraction, down or up), a port of scipy's `resample_poly` (a Kaiser-windowed sinc, beta 5, through `upfirdn`), the method nxr-convert makes a recording's resampled copy with. `tResample` holds it to scipy on six ratios (1e-12); on a real store it reproduces the stored 300 Hz copy from the 1200 Hz recording to float32 rounding (`tNxrStoreParity`). MATLAB's own `resample` designs its filter differently and is not this.
- **The editor's spectral chain, sensors to flow, in MATLAB.** `rheome.eigen.hodge_package` (the modes as vertex potentials, merged by mu, single), `rheome.inverse.hodge_weights` (the kernel carried onto the bridge and fused into the modes: B [2K x C]), `rheome.inverse.hodge_synthesis` (the current from the weights), `rheome.spectral.root` (a chunk's one-sided spectra, the editor's load), `rheome.spectral.column` (the real and analytic column at a time, under a gain) and `rheome.spectral.meyerband` (the band gain). Checked: B rebuilt from a real subject's kernel against the editor's stored B (whole matrix 4e-6, a typical row 3e-8: float32 storage; the worst rows at nearly coincident eigenvalues), the root and columns against the editor's own functions run on the same chunk (`tests/editor_spectral.mts`; equal to float32 rounding), and the synthesis read back to its weights (1e-15). Tests: `tSpectralRoot`, `tGlueHodgeParity/fusedWeightsReadBackTheirSynthesis`, `tNxrStoreParity`. The catalogue is now read in-process through MATLAB's Python (`sqlite3`, read-only), not the command-line tool.
- **The editor's spectral-chain inputs are read from an nxr store.** `rheome.io.read.zarr` reads any rectangular part of an uncompressed zarr v3 array, single-chunk, a regular grid or `sharding_indexed` shards (only the inner chunks a range touches, by the shard index; an absent chunk is the fill value; a compressed array is refused). On top of it, with the editor's own lookups: `nxrrecording` (any span of samples, fs and t0 from the time Line, channel names, bad flags), `nxrkernel` (the kernel, its channels as recording rows, the recording it applies to, and its stored fused weights B with the eigenvalues and glued mesh they were made with) and `nxrhodge` (a stored basis). `tZarrRead` holds the reader to arrays written by zarr-python (synthetic fixtures, `tests/fixtures/zarr/make.py`); `tNxrStoreParity` reads a recording chunk, the kernel and B against zarr-python on a real store and rebuilds B's basis (glue, flips, solve), its eigenvalues against those B records. `nxrsurface` now reads through `rheome.io.read.zarr`.
- **`rheome.io.read.nxrsurface` reads a subject's cortex from an nxr store** (the editor's dataset folder: `dataset.sqlite` and `<subject>.nxr.zarr`): the primary surface (or `Surface=`), its faces and its medial wall (the Desikan-Killiany `unassigned` vertices), by nxr-convert's own lookups. The catalogue is read with the `sqlite3` command-line tool (MATLAB's `sqlite` needs the Database Toolbox) and each array from its chunk file; only uncompressed single-chunk zarr v3 arrays are read, any other layout is refused. `tNxrStoreParity` runs the whole chain (read, glue, flips, solve of 200 modes) against nxr-convert's on every subject of `$RHEOME_NXR_STORE`, read-only; nothing from the store is kept in the repository.
- **The editor's surface preparation and Hodge-Dirac solve ported** (nxr-cortical-flow, via nxr-convert `eigen/`): `rheome.geom.glue` (the hemispheres glued by anatomy: each medial wall cut as a disk, a bridge of rings grown and stitched), `rheome.geom.meshkey` (the editor's mesh fingerprint) and `rheome.eigen.hodge_dirac` (the exact and coexact ladders, one lambda cut). Parity in `tGlueHodgeParity`: on the editor's golden mesh the same faces, kept vertices, bridge points (to 1e-15), mesh key, flips and pencil (K to 1e-12, M to 1e-9), and nxr-convert's solve (lambda to 1e-9, modes to 1e-6); on the editor's `cortex_pial_low` (38 678 glued faces) the glue is bit-identical to nxr-convert's, the flips the same, no edge clamped, and lambda to 1e-12.

- **`rheome.geom.draw.line` no longer takes the heat solver: drawing and measuring are separate (Diellor, 2026-11-10).** New `rheome.geom.draw.flipMesh(S)` builds the mesh FlipOut reuses (the Editor's `flipMesh`; what `flipOut` built and returned as `info.cache`). `draw.line(S, from, to, mesh)` and `draw.flipOut(S, from, to, mesh)` take it; `info.mesh` carries it for the next path (`flipOut` keeps `info.cache`, the same struct, for one release; anything else as `flipOut`'s `mesh`, or a mesh of another surface, is `rheome:flipOut:mesh`). For one release a heat solver (anything not a flip mesh) as `line`'s 4th argument, or `geodesicpath`'s, is ignored, warns once (`rheome:deprecated`, naming `draw.flipMesh`) and comes back as `info.solver` unchanged; `info.solver` no longer holds the halfedge mesh. Lines are unchanged. `rheome.geom.measure.heatDistanceEigen` without `opts.t` no longer crashes on a `face_gradient` bundle (it has no `.Vertices`): the mean edge length is read from the face gradients (|p2 - p1| = 2A |grad phi_3|, the same h to rounding), from `.Vertices` when given; with neither, error `rheome:heatDistanceEigen:noT` asks for `opts.t` (the mass-derived spacing fallback is gone).

- **`rheome.geom`'s surface functions grouped by job into three sub-packages, the Cortical Flow Editor's and nxr's words (Diellor, 2026-11-10): `measure` returns numbers, `select` returns a set of vertices, `draw` returns a curve on the surface.** `rheome.geom.geodesic` -> `rheome.geom.measure.heatDistance`, `heatdistance` -> `measure.heatDistanceEigen`, `karcher_mean` -> `measure.karcherMean`, `tilemean` -> `measure.tileMean`, `cortexmask` -> `select.cortex`, `flipOut` -> `draw.flipOut` (with `LICENSE-geometry-central` and its private helpers, now in `+rheome/+geom/+draw`), `geodesicpath` -> `draw.line`; same inputs and outputs, error and warning ids unchanged. Mesh structure and bookkeeping (icosphere, edgegraph, tree, tiles, the ladders, sphereframe, spherepatches, spiralorder, cortexnodes) stay in `rheome.geom`. The old names forward and warn once per session (`rheome:deprecated`); every caller in Rheome uses the new names (`tests/tGeomByJob`). The generated reference (`doc/reference`) is not regenerated here.

- **`rheome.geom.flipOut`: geometry-central's FlipOut geodesic (Sharp & Crane 2020, "You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges"), ported to pure MATLAB; `rheome.geom.geodesicpath` now traces it (Diellor, 2026-11-10: every tool draws the same line).** A port of geometry-central at 8579865 (nxr-compute's dependency, MIT; notice in `+rheome/+geom/LICENSE-geometry-central`): the mesh built from the face list in geometry-central's order, Dijkstra's edge path with its tie-breaking, `FlipEdgeNetwork::iterativeShorten` on a signpost intrinsic triangulation, and each intrinsic edge traced back onto the input faces (`traceGeodesic`), as nxr-compute's `tracePathSurface` calls them. Against the released viewer's WASM (`tests/goldens/flipout`, `tFlipOut`; 40 + 40 pairs on two folded spheres, 60 on fsaverage5 lh pial) every pair has the same points, within 1e-14 m (tolerance 1e-9), the same vertex/edge placement and the same length (3e-16 relative). fsaverage5: ~0.03 s a path once the mesh is built (`info.cache`). `geodesicpath` keeps its signature; `info` has `.length`, `.faces`, new `.vertices`/`.weights` (where each point is), `.solver` (now the halfedge mesh, pass it back); `.distance` (the heat distance) and `.fallbacks` are gone, and the heat-descent tracer is retired (`rheome.geom.geodesic`, the heat distance, stays). Paths move: FlipOut's path is locally shortest, the old descent followed the heat distance's gradient.

- **`rheome.snapshot` takes the viewer's phase-3 shot fields (nxr-cortical-flow #139; nxr's viewer-twin phase 3, T2; 2026-11-09), as `nxr.snapshot` does, job byte for byte.** New shot fields: `Cortex` (true/false/opacity), `Labelling` (a labelling of `Labellings=` painted), `Outlines` (`Labelling`, `Coarser`, `Width`, `CoarserWidth`, `Select`, `SelectStyle`), and RAW PICTURES `StagePx` [w h], `Fit` (`Margin`, `Shared`, or true), `Dpi`. A shot may paint no field (`Field=""`, or none given with a `Labelling` or no real field): no colourbar, and `Unit`, `Label`, `Reference`, `Quantity`, `CLim`, `Colormap` refused. `View` may be a list; a strip may have several views (a row a view; it was refused). New options `ScaleBar="auto"`, `ScaleBarKind` (`geodesic` | `plane`), `ScaleBarTile` (1-based). New `out` fields, present when a shot of the call asks: `labelling`, `outlines`, `stops`, `fit`, `raw` (each raw render: view, camera, 1-based frame, size, mm per pixel, fit, colourbar, and its `png` and `marks` files in `<Out>/<Name>/`). Raw pictures without `Width` have no panel (`out.svg` ''). New errors `rheome:snapshot:outlines`, `rheome:snapshot:raw`. A bad `ScaleBar` is `rheome:snapshot:input` (was the input parser's). Calls without the new options send the same job as before.
- **`rheome.snapshot` is a thin caller of the Cortical Flow viewer's snapshot panel (phase 2 of nxr's viewer-twin spec, decisions D1-D5; 2026-11-09).** It keeps its options, shot fields, checks and `rheome:snapshot:*` errors, writes the bundle and a job (`nxr-snapshot-job@1`, byte-identical to `nxr.snapshot`'s) and runs the viewer's own driver, `<viewer>/driver/snapshot-driver.mjs` (a viewer build of nxr-cortical-flow ad9d824 or later; `rheome:snapshot:driver` otherwise); limits, visibility, the scale bar, the paths, the selections, the inset and the panel's composition are the viewer's (by Rheome's rules). Changed: the panel is the viewer's print figure (the key band on top: colourbar and a key of every mark, then the views); a path is the viewer's exact geodesic (geometry-central), no longer the heat-method descent, and the bar's curve is that geodesic cut at L (its placement Rheome's rule); `Paths.Points` (a given path) is refused (`rheome:snapshot:path`) until the viewer draws given paths; an inset is the viewer's camera on the region (far above the old 4x cap); `out.cameras` and `out.settings` are gone, `out.inset.px` is gone, `.scalebar.px`/`.paths.px` are the raw render's pixels with the origin at its top-left corner; the marks are one `<Out>/<Name>.marks.json` per shot (was one per render), with `<Name>.report.json` and `<Out>/job.json` beside; a strip's renders are named by 0-based frames (`<view>_f0.png`, was `_f1`); `Scale` is ignored. `+rheome/private/snapshot_driver.mjs` is deleted. `tSnapshot` checks the reported pixels against where the viewer paints the vertices.
- **`rheome.cores` lagged pair, share, S3 and the Wishart null, ported from nxr-analysis (`src/nxr/vortex.py` pair/share/s3/wishart/cov/dof_fraction/sqrt_psd/analytic, commit 3da1320; Diellor, 2026-11-02: one pipeline, one mathematics; step 10a).** New: `pair` (reduced SVD, smallest-eigenvalue eigenvector, fb sign so Im(fa Sx fb) >= 0, u = p + i q), `pairgauge` (the canonical gauge of nxr's goldens, since the eigenvector phase is free and MATLAB's eig differs from numpy's), `pairshare`, `pairs3`, `pairnull` (k = max(half-even round(dof), 2); draws are an input, [C x k x 2], or a RandStream), `paircov`, `pairdof`, `pairsqrt`, `analytic`. Parity with nxr's synthetic goldens (`tests/goldens/vortex-pair`, `tVortexPair`): see the test for tolerances; raw pair values are compared in the canonical gauge, the plane, form, share, S3 and the Wishart p values with the stored draws replayed.
- **Every divergence / curl option nxr-analysis has, matched to nxr's goldens (Diellor, 2026-11-02: "keep all three. This is a toolbox. We need all the options"; one pipeline, one mathematics).** `rheome.differential.curl` and `divergence` take the 4th argument as before (curl `true` = tangent, `false` = ambient; divergence `false` = ambient, `true` = tangent; defaults unchanged, no number moves) and now also the names `"tan"`, `"3d"` and `"curv"`; `"curv"` is the operator applied to the normal part J_n n alone (J_n = J.n, area-weighted vertex normal). New: `rheome.differential.meancurvature` (H = -1/2 (DX.n), barycentric-area cotangent Laplacian, +1/R on a sphere with outward normals; nxr `FaceOperators.mean_curvature`; `rheome.scale.measure_ambient`'s inline H differs by the lumped mass it divides by, so it is not reused), `curvaturesource` (2 H J_n; nxr `div_curv_2h`) and `divcurl` (the six maps `div_tan, div_curv, div3d, curl_tan, curl_curv, curl3d` in one call; nxr `div_curl_all`). Exactly, to rounding (7e-16 / 1.4e-15): div3d = div_tan + div_curv, curl3d = curl_tan + curl_curv. Parity with nxr's synthetic goldens (`tests/goldens/div-options`, ellipsoid 70/55/45 mm, `tDivOptions`; positions V_mm/1000 cast to single as nxr does): all six maps, H and 2 H J_n agree to <= 7.5e-14 of the map's maximum (curl_curv); without the single cast <= 1.1e-5 (curl_curv; the rest <= 1e-5, where float32 positions of nxr differ by that much). The curvature term through the Dirac: `diracparts` of J_n n returns div_curv and curl_curv to 1e-15 (the intrinsic Dirac reads the ambient operators, so its divergence is div_tan + div_curv); and the extrinsic-minus-intrinsic spin Dirac (D_HY - D_ext = H_f, the dihedral-angle integrated mean curvature per face) gives face-mean div_curv = 2 (H_f/A_f) q.n_f with q = J_n n to 3.0 % (the cotan 2 H J_n differs from it by 1.5 %); the factor 2 is the divergence of a normal current being the trace kappa_1 + kappa_2 of the shape operator, where the Dirac's H is the mean (kappa_1 + kappa_2)/2.
- **`rheome.cores` surface ladder, tile crossings, path consistency R, repeat, tile dot, tile p/q, BH and place-based direction consistency, ported from nxr-analysis (`src/nxr/tile_connectome.py`, `src/nxr/vortex_cores.py`, `src/nxr/vortex.py`, commit 0db82e4; Diellor, 2026-11-02: one pipeline, one mathematics).** New: `ladder` (two-piece mesh + hemisphere, codes 2^level*h + t, Dijkstra on Euclidean edge lengths, stable split, one connected piece per tile), `crossings` (float32 steps, nearest support vertex, unrenormalised tangent projection, <= 1e-9 dropped, row 2*tile + (sense > 0) + 1, ordered by (row, track)), `pathR` (n >= 5), `pathchance` (shuffles among crossings; permutations are an input or MATLAB draws), `senses`, `repeat` and `tiledot` (n >= 10, senses pooled), `tilepq`, `bh`, `placeconsistency` (n >= 10, no per-track averaging). `rheome.geom.surfaceladder` (nsp's ladder), `ladder` and `tiles` are unchanged: nxr's ladder differs from surfaceladder (takes both hemispheres and a hemisphere code; has the one-piece step; a tile whose first vertex is out of the far end's reach sets b differently; edge lengths and areas by plain sqrt-of-squares as numpy, since ties in the bisection key are decided by the last bit). On the goldens' mesh surfaceladder agrees with nxr up to level 5-6 and then differs on 2-29 vertices per level (levels 7-8). Parity with nxr's synthetic goldens (`tests/goldens/vortex-tiles`, `tVortexTiles`): ladder codes exact (all 9 levels), crossings (counts, order, vectors <= 1e-12), R and shuffled chance <= 1e-12, place consistency <= 1e-12; repeat, tile dot, R from the stored float32 sums agree to float32 rounding (<= 1e-6; nxr computes them on float32 sums), tile p/q and BH exact.
- **`rheome.cores` burst stages, constant-Q band frames and the envelope-matched lag-free null, ported from nxr-analysis (`src/nxr/vortex.py`, commit 30b9622; Diellor, 2026-11-02: one pipeline, one mathematics).** New: `morletgain` (peak 2, 0 at f <= 0), `bandframes` (complex single frames at fps, bins 1..k-1 copied, m/n scaling), `boxsmooth` (zero-padded `np.convolve` 'same', even-k window [i-k/2, i+k/2-1]; not `movmean`), `linpercentile` (numpy's linear rule; `prctile` differs), `burstenvelope` (nxr `spindle_envelope`: hi 1.5 / lo 1.0 robust SD, gates in cycles, no zero-MAD guard, no maximum length), `burstfit` (nxr `risefall`, quadrature model, `lsqcurvefit` with 1e-10 tolerances), `burststages` (onset/within/offset masks and fits), `burstmask` (older percentile mask), `lagfreenull` and `envelopenull` (the white-noise draw is an input, or a `RandStream`, since numpy's generator cannot be reproduced). `rheome.detect.spindle` and `rheome.detect.risefall` (sleep-spindle detectors) are unchanged: nxr's versions depart from them (thresholds, cycles instead of seconds, banker's rounding, numpy percentile, ub - 1e-12 start clip) and live under their own names. Parity with nxr's synthetic goldens (`tests/goldens/vortex-signal`, `tVortexSignal`): band frames bit-identical, amplitude stages <= 1.2e-7 of peak, thresholds 9e-10, candidates, masks and the old burst mask frame for frame, the four fits within 1.1e-5 s on tOn/tOff (tau 6e-5 relative, r2 6e-10), null stages <= 1.4e-7 of peak.
- **`rheome.cores` parity extended to nxr's hard-case goldens (Diellor, 2026-11-02).** `tests/goldens/vortex-hard` (nine synthetic scenarios from nxr-analysis PR #7: greedy competition, distance tie, frame gaps and missing frames, recurrence, reversal, the float32 threshold bin edge) and `tVortexCoresHard` check patches, steps, scalars (null = NaN), histograms and occupancy; no change to `rheome.cores` was needed. The first set's manifest and README are refreshed from nxr (curl labelled per metre; no numbers changed).
- **Welch spectra use the periodic Hann window (Diellor, 2026-11-02).** `rheome.scale.cleanwelch` and the sensor spectrum in `rheome.modes.fields` now use `hann(nf,'periodic')`, w(i) = 1/2 - 1/2 cos(2 pi i/nf), as nxr-analysis and the Cortical Flow Editor do (was MATLAB's symmetric `hann(nf)`). Welch PSDs shift slightly (window length nf-1 vs nf; a few 1e-4 relative); `cleanwelch` still equals the matching `pwelch` call bit for bit. The remaining Hann-Welch callers now use it too (2026-11-02, same decision): `scale.measure_patterns`, `scale.bandsnr`, `scale.measure_aperiodic`, `scale.measure_bandperiodic`, `flow.vortexspectrum` (`hann(n,'periodic')`), and the reference spectra in `tSpectralSynth` and `tSpectralSplit`. Their PSDs shift by the same few 1e-4 relative. Not Welch/Hann-Welch, unchanged: `forward.simulate`, `scale.measure_noisefloor` envelopes, `filters.firbandpass`, `flow.psd` (Hamming), `flow.joint` (optional taper), `detect.spindle`, `ingest.peaks`.
- **`rheome.differential.curl` now projects the field onto the tangent plane first, as nxr-analysis and the Cortical Flow Editor do.** Same pipeline, same mathematics: each vertex's ambient J loses its normal part (J - (J.n)n, n = the area-weighted vertex normal) before the per-face curl and the area-weighted face-to-vertex mean (nxr `FaceOperators.div_curl`, Editor `flowFieldCpu` with `tangential`). Outputs of `curl` change slightly: on the vortex goldens the old raw curl differed from nxr's `curl_vertex` by 1.5e-4 relative, the new one by 4.4e-8; on a folded cortex expect changes of that order (more where vertex normals differ strongly between neighbours). Callers that move: `differential.helmholtz` (`.Curl`), `detect.criticalPoints` (vorticity, so `detect.vortex` / `detect.track` classification at the margin), `demos.sensor_topology`. `curl(J,S,fg,false)` gives the raw ambient curl (used by `differential.diracparts`, whose curl is the Dirac's and keeps its exact energy identity). `rheome.differential.divergence` is unchanged by default (the Editor's `tangential` is off by default and the ambient divergence deliberately carries the mean-curvature term); new optional `divergence(J,S,fg,true)` projects as nxr's `div_curl` does. Positions are not cast to float32 (nxr and the Editor do): the effect is ~1e-8 relative, inside the 1e-6 parity test.
- **`rheome.cores`: the threshold-on-curl vortex core detector and tracker, ported from nxr-analysis (`src/nxr/vortex.py`, commit 903fdc2).** Beside the critical-point route (`rheome.detect.vortex/criticalPoints/track`, unchanged), one pipeline and one mathematics with nxr and the Cortical Flow Editor. `heatgain`, `hatgain` (level hat, peak 1), `modecoords`, `operator` (P'((g/s)Z) in single), `threshold` (4096-bin log histogram, 98th percentile), `patches` (same-sense connected components, nxr order), `adjacency`, `stats` (greedy 15 mm linking, life/path/displacement/straightness/speed, half-cycle reversal and persistence, cycle recurrence, Jaccard, histograms, occupancy). Parity with nxr's synthetic goldens (`tests/goldens/vortex`, copied from nxr 9d97761): `tVortexCores`; every stage matches to float precision, patches and tracks exactly. `rheome.differential.curl` equals nxr's `curl_vertex` (4e-8) when the mesh is in metres and the flow is first projected onto the tangent plane.
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
  the heat distance (replaced by FlipOut, `rheome.geom.flipOut`, 2026-11-10). On fsaverage6 its median length is within 0.5 % of the exact geodesic, and the board's
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
