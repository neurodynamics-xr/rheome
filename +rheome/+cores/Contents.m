% CORES  The threshold-on-curl vortex core detector and its tracker (a port of nxr src/nxr/vortex.py).
%
% Cores here are connected patches of the smoothed curl whose |curl| passes a high percentile, one rotation sense per
% patch, tracked frame to frame. This is NOT the critical-point route of rheome.detect (criticalPoints, vortex, track).
% Ported function for function from nxr-analysis (commit 903fdc2 and the goldens of 9d97761): rheome, nxr and the
% Cortical Flow Editor are one pipeline and one mathematics. Operator stages are float32, as in nxr.
%
% Functions:
%   rheome.cores.adjacency   - Symmetric vertex adjacency of a mesh (nxr support_adjacency).
%   rheome.cores.hatgain     - Level hat gain on the Laplace-Beltrami eigenvalues, peak 1 (nxr hat_operator).
%   rheome.cores.heatgain    - Heat gain exp(-lambda*sigma^2/2) (nxr core_operator).
%   rheome.cores.modecoords  - A vertex field in mass-orthonormal mode coordinates (nxr Cortex.coords).
%   rheome.cores.operator    - The core operator P'*((g/s).*Z) in single (nxr core_operator, hat_operator).
%   rheome.cores.patches     - Cores of a curl movie: frame, sense, centre, area (nxr _core_patches).
%   rheome.cores.stats       - Core tracks and track measures of a curl movie (nxr core_stats).
%   rheome.cores.threshold   - The 98th-percentile |curl| threshold from a log histogram (nxr _core_threshold).
%
% Author: Diellor Basha, 2026
