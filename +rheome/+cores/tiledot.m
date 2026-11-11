function d = tiledot(sumsA, nA, sumsB, nB, nMin)
% CORES.TILEDOT  Per tile (senses pooled), the dot product of the two runs' mean resultant vectors (nxr tile_dot).
%
%   d = rheome.cores.tiledot(sumsA, nA, sumsB, nB)          % nMin = 10 (nxr NMIN_TILE)
%
% sumsA, sumsB: the runs' tile-sense sums as rheome.cores.crossings gives them (tangent vectors AT THE TILE CENTRES since nxr 1.15.0: the two runs' vectors lie
% in the same plane). d [T x 1] = (sA/nA) . (sB/nB), each run's mean resultant vector (not unit); NaN where either run has n < nMin. One null draw at a time. Ports nxr-analysis src/nxr/vortex_cores.py tile_dot (line 90, commit 0db82e4).
%
% See also: rheome.cores.repeat, rheome.cores.tilepq
%
% Author: Diellor Basha, 2026
    if nargin < 5, nMin = 10; end
    [sa, na] = rheome.cores.senses(sumsA, nA);  [sb, nb] = rheome.cores.senses(sumsB, nB);
    ok = na >= nMin & nb >= nMin;
    d = sum(sa ./ max(na, 1) .* sb ./ max(nb, 1), 2);
    d(~ok) = NaN;
end
