function [p, q] = tilepq(dot, nullDots)
% CORES.TILEPQ  Permutation p and BH q of the tile dots against a null of independent runs (nxr vortex.py:203-204, inline there).
%
%   [p, q] = rheome.cores.tilepq(dot, nullDots)
%
% dot [T x 1] (rheome.cores.tiledot of the two runs), nullDots [D x T] the tile dots of D null pairs. ok = finite dot;
% p = (#null dots >= dot + 1) / (#finite null dots + 1), NaN where not ok; q = rheome.cores.bh(p(ok)), NaN elsewhere. Significant at q < 0.05.
% Ports nxr-analysis src/nxr/vortex_cores.py (lines 203-204, commit 0db82e4; the golden replicates it).
%
% See also: rheome.cores.tiledot, rheome.cores.bh
%
% Author: Diellor Basha, 2026
    dot = dot(:);  ok = isfinite(dot);
    ge = nullDots >= dot.';                                 % NaN compares false
    p = (sum(ge, 1).' + 1) ./ (sum(isfinite(nullDots), 1).' + 1);
    p(~ok) = NaN;
    q = NaN(size(p));  q(ok) = rheome.cores.bh(p(ok));
end
