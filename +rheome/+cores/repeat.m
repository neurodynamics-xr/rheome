function r = repeat(sumsA, nA, sumsB, nB, nMin)
% CORES.REPEAT  Repeatability of the tile directions between two runs: the count-weighted mean cosine, senses pooled (nxr repeat).
%
%   r = rheome.cores.repeat(sumsA, nA, sumsB, nB)          % nMin = 10 (nxr NMIN_TILE)
%
% Tiles with n >= nMin in BOTH runs, the two senses pooled (rheome.cores.senses); per tile the cosine of the two pooled sum vectors,
% weight min(nA, nB); r = sum(c*w)/sum(w), NaN when no tile qualifies. Ports nxr-analysis src/nxr/vortex_cores.py repeat (line 83,
% commit 0db82e4).
%
% See also: rheome.cores.tiledot, rheome.cores.tilepq
%
% Author: Diellor Basha, 2026
    if nargin < 5, nMin = 10; end
    [sa, na] = rheome.cores.senses(sumsA, nA);  [sb, nb] = rheome.cores.senses(sumsB, nB);
    m = na >= nMin & nb >= nMin;
    if ~any(m), r = NaN; return; end
    c = sum(sa(m, :) .* sb(m, :), 2) ./ (vecnorm(sa(m, :), 2, 2) .* vecnorm(sb(m, :), 2, 2));
    w = min(na(m), nb(m));
    r = sum(c .* w) / sum(w);
end
