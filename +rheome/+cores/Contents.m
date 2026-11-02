% CORES  The threshold-on-curl vortex core detector and its tracker (a port of nxr src/nxr/vortex.py).
%
% Cores here are connected patches of the smoothed curl whose |curl| passes a high percentile, one rotation sense per
% patch, tracked frame to frame. This is NOT the critical-point route of rheome.detect (criticalPoints, vortex, track).
% Ported function for function from nxr-analysis (commit 903fdc2 and the goldens of 9d97761): rheome, nxr and the
% Cortical Flow Editor are one pipeline and one mathematics. Operator stages are float32, as in nxr.
% Burst stages, band frames and the lag-free envelope-matched null (commit 30b9622, goldens tests/goldens/vortex-signal) are
% the vortex pipeline's own; rheome.detect.spindle and rheome.detect.risefall (sleep spindles) are unchanged.
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
%
% Author: Diellor Basha, 2026
