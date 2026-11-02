function p = linpercentile(x, q)
% CORES.LINPERCENTILE  Percentile with numpy's default linear rule, np.percentile(x, q).
%
%   p = rheome.cores.linpercentile(x, q)       % q in [0, 100]
%
% The q-th value sits at sorted position (n-1)*q/100 (0-based) and is interpolated linearly between its neighbours.
% MATLAB's prctile places the i-th sorted value at (i-0.5)/n instead and differs, most on short vectors. The
% interpolation is numpy's lerp (from the upper neighbour when the fraction is >= 0.5).
%
% Ports the percentile of nxr-analysis src/nxr/vortex.py (burst_mask, risefall; commit 30b9622).
%
% See also: rheome.cores.burstfit, rheome.cores.burstmask
%
% Author: Diellor Basha, 2026
    arguments
        x {mustBeNumeric, mustBeNonempty}
        q (1,1) double {mustBeInRange(q, 0, 100)}
    end
    s = sort(double(x(:)));
    n = numel(s);
    pos = (n - 1) * q / 100;
    lo = floor(pos);  f = pos - lo;
    a = s(lo + 1);  b = s(min(lo + 2, n));
    d = b - a;
    if f >= 0.5, p = b - d * (1 - f); else, p = a + d * f; end
end
