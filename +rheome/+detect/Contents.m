% DETECT  Feature detection on source vector fields.
%
% Algorithms that chain the module's operators into a detector. The first is the
% Dirac-connection route to critical points:
%
%   rheome.detect.operator(S)                    - PRECOMPUTE the mesh-only operators the detectors
%       reuse (whole-surface face gradient + per-hemisphere connection Laplacian). Build once,
%       pass as the trailing 'op' argument to rheome.detect.criticalPoints / rheome.detect.vortex to skip
%       rebuilding them every frame -- essential for time series. Result is unchanged; only the
%       repeated full-resolution operator construction is removed.
%
%   rheome.detect.criticalPoints(J, S [,types][,op]) - classified critical points of an ambient current.
%       Chains rheome.operators.connection_laplacian (Levi-Civita) -> complex tangent field ->
%       per-triangle winding number (integer charge, Poincare-Hopf exact) -> classify each
%       via local div/curl into vortex / source / sink / saddle, ranked by strength. The
%       'types' argument returns only the requested kinds. See the scripts vortex_analysis.m
%       (vortices on real MEG) and detect_validate.m (ground truth on analytic sphere fields).
%
%   rheome.detect.vortex(J, S [,props][,op])     - PHYSICAL PROPERTIES of the vortices (measurement).
%       Detects vortices (via rheome.detect.criticalPoints) then measures each in real units: size
%       (mm, radius of peak tangential speed = sigma), chirality (+1 CCW / -1 CW from outside,
%       via (r x v).n_outward), strength (circulation), position. The 'props' argument (default
%       'all') selects which to compute. Calibrated on the sphere: known sigma/handedness read
%       back correctly (10/15/20/30 mm -> 11/15/19/29 mm; CCW=+1, CW=-1).
%
%   rheome.detect.track(J, S [,opts])            - TRACK the critical points across time (Lagrangian).
%       Per frame rheome.detect.criticalPoints, then a ByteTrack-style two-stage association into
%       trajectories: match STRONG detections, RECOVER open tracks from WEAK ones, and hold a
%       missed track 'lost' for a few frames (lostTrackBuffer) before ending it -- so flickering
%       motifs stay one track instead of fragmenting; merge/split events are flagged (flow-
%       singularity pair bifurcation/annihilation). Each track carries velocity (m/s), lifetime,
%       straightness (1 traveling / 0 standing), state, birth/death, gaps and an intrinsic home
%       location (rheome.geom.karcher_mean). Association ported from the sibling repo nxr-tracking
%       (tracking/meshByteTrack.ts), adapting ByteTrack [Zhang et al., ECCV 2022] to the cortex.
%       Validated on a translating vortex + a gap-recovery case (feature_tracking_validate.m);
%       demo in feature_tracking.m (motifs on real MEG alpha).
%
%   rheome.detect.peaks(X, S [,opts])            - one-ring LOCAL MAXIMA of a scalar map per frame, with
%       sub-edge position, value, local prominence, half-height width, optional tile address, and
%       non-maximum suppression. Emits rheome.detect.track's frameFeatures. ⭐ Peaks are found at vertices
%       and THEN addressed to tiles: a tile's maximum is not a peak. plant_track_omega.m recovers a
%       planted speed at slope 0.99-1.01 from 0.05 to 1 m/s through 10 dB.
%   rheome.detect.tilepeaks / rheome.detect.tiletrack - the peak as a TILE and the trajectory as a walk on the
%       tile graph (adjacent steps only, global assignment, margin hysteresis). Tile size is chosen
%       from the feature (rheome.geom.jointcell); plant_tiletrack_omega.m is the validation.
%   rheome.detect.tilepath - ⭐ TRACK-BEFORE-DETECT: Viterbi on the tile graph within one observation window,
%       scoring Y - baseline along every admissible walk before deciding. At 10 dB it keeps the movers
%       detect-then-link loses (plant_tilepath_omega.m). Use it, not tiletrack, below ~20 dB.
%   rheome.detect.trackstats / tracknull / trackselect - per-track duration, net displacement, straightness,
%       net speed; thresholds from noise-only windows at a false-trajectory RATE alpha per second.
%   rheome.detect.blobscale  - the characteristic SIZE of structure, by scale selection over a
%                       band-pass ladder on the surface's own operator. ⚠ the response must
%                       be scale-normalised or every blob measures as small; rheome.filters.mexhat
%                       already carries the factor. ⚠ a negative eigenvalue detonates it, so
%                       they are clamped with a warning. MinWavelength floors the ladder at
%                       the instrument's resolution, so a size below it is never reported.
%
% Author: Diellor Basha, 2026
