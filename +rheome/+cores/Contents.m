% CORES  The threshold-on-curl vortex core detector and its tracker (a port of nxr src/nxr/vortex.py).
%
% Cores here are connected patches of the smoothed curl whose |curl| passes a high percentile, one rotation sense per
% patch, tracked frame to frame. This is NOT the critical-point route of rheome.detect (criticalPoints, vortex, track).
% Ported function for function from nxr-analysis (commit 903fdc2 and the goldens of 9d97761): rheome, nxr and the
% Cortical Flow Editor are one pipeline and one mathematics. Operator stages are float32, as in nxr.
% Burst stages, band frames and the lag-free envelope-matched null (commit 30b9622, goldens tests/goldens/vortex-signal) are
% the vortex pipeline's own; rheome.detect.spindle and rheome.detect.risefall (sleep spindles) are unchanged.
%
% Vectors from different places are CARRIED to the tile's centre / the place (the vector heat method's transport, rheome.geom.GeometryKit) before they are
% added (nxr 1.15.0; goldens regenerated). A raw 3-D sum on a folded tile cancels what points the same way along the surface.
%
% Tiles (commit 0db82e4, goldens tests/goldens/vortex-tiles): the surface ladder, tile crossings, path consistency R, repeat, tile dot, tile p/q, BH and
% the place-based direction consistency. Permutations are inputs (numpy's generator cannot be reproduced) or MATLAB draws.
%
% Functions:
%   rheome.cores.adjacency   - Symmetric vertex adjacency of a mesh (nxr support_adjacency).
%   rheome.cores.bandframes  - A band's analytic signal at fps frames a second, complex single (nxr band_at_frames).
%   rheome.cores.boxsmooth   - Zero-padded moving mean, numpy np.convolve(x, ones(k)/k, 'same') (not movmean).
%   rheome.cores.burstenvelope - Candidate bursts of an envelope, durations in cycles (nxr spindle_envelope).
%   rheome.cores.burstfit    - Quadrature double-sigmoid rise/fall fit of an episode (nxr risefall).
%   rheome.cores.burstmask   - Frames inside bursts by a percentile threshold (nxr burst_mask).
%   rheome.cores.burststages - Onset, within and offset masks and fits of the bursts (nxr burst_stages).
%   rheome.cores.envelopenull - Lag-free surrogate rescaled to the data's envelope (nxr envelope_matched_null).
%   rheome.cores.heatgain    - Heat gain exp(-lambda*sigma^2/2) (nxr core_operator).
%   rheome.cores.lagfreenull - Lag-free Gaussian surrogate of a band's real signal; the draw is an input (nxr lagfree_real).
%   rheome.cores.linpercentile - numpy's linear np.percentile (prctile differs).
%   rheome.cores.modecoords  - A vertex field in mass-orthonormal mode coordinates (nxr Cortex.coords).
%   rheome.cores.morletgain  - Morlet gain per FFT bin, peak 2, zero at f <= 0 (nxr morlet_gain).
%   rheome.cores.operator    - The core operator P'*((g/s).*Z) in single (nxr core_operator, hat_operator).
%   rheome.cores.patches     - Cores of a curl movie: frame, sense, centre, area (nxr _core_patches).
%   rheome.cores.stats       - Core tracks and track measures of a curl movie (nxr core_stats).
%   rheome.cores.threshold   - The 98th-percentile |curl| threshold from a log histogram (nxr _core_threshold).
%   rheome.cores.bh          - Benjamini-Hochberg q values (nxr bh).
%   rheome.cores.crossings   - Per tile-sense summed crossing directions of core steps at one ladder level (nxr Tiles.crossings).
%   rheome.cores.ladder      - Two-piece surface ladder, codes 2^level*h + t, one piece per tile (nxr surface_ladder).
%   rheome.cores.Pooling     - Vectors pooled at anchors: carried to the anchor's tangent plane, then added (nxr measure.vectors.Pooling).
%   rheome.cores.pathR       - Path consistency R, count-weighted mean resultant length (nxr mean_R).
%   rheome.cores.pathchance  - Chance R with the crossing directions shuffled among crossings (nxr stage_measures).
%   rheome.cores.placeconsistency - Direction consistency at places, no per-track averaging (nxr direction_consistency).
%   rheome.cores.repeat      - Count-weighted mean cosine of tile directions between two runs (nxr repeat).
%   rheome.cores.senses      - Pool the two senses of each tile (nxr senses).
%   rheome.cores.tilecentres - The member vertex nearest each tile's centroid (nxr tile_centres).
%   rheome.cores.tilepool    - The pooling of one ladder level's tiles: vectors carried to the tile centres.
%   rheome.cores.placepool   - The pooling of the direction consistency: vectors carried to the nearest place (nxr place_pooling).
%   rheome.cores.tiledot     - Per-tile dot of two runs' mean resultant vectors (nxr tile_dot).
%   rheome.cores.tilepq      - Tile permutation p and BH q against a null of runs (nxr vortex.py:203-204).
%
% Lagged pair (commit 3da1320, goldens tests/goldens/vortex-pair): the pair of a reading under a cross-spectrum, its share, S3 and the Wishart null.
%   rheome.cores.analytic    - Decimated Morlet analytic signal, complex double (nxr analytic).
%   rheome.cores.pair        - Lagged pair fa, fb, u of an operator under Sx (nxr pair); phase free, see pairgauge.
%   rheome.cores.paircov     - Masked cross-spectrum, no mean removal (nxr cov).
%   rheome.cores.pairdof     - Effective independent samples per sample (nxr dof_fraction).
%   rheome.cores.pairgauge   - Canonical gauge of a pair, for comparing ports (not in nxr).
%   rheome.cores.pairnull    - Lag-free Wishart null covariance; the draws are an input (nxr wishart).
%   rheome.cores.pairs3      - Poincare S3 of a pair (nxr s3).
%   rheome.cores.pairshare   - Lagged share ||L Im(Sx) L'||/||L Sx L'|| (nxr share).
%   rheome.cores.pairsqrt    - Symmetric PSD square root (nxr sqrt_psd).
%
% Author: Diellor Basha, 2026
