function R = pathR(sums, n, nMin)
% CORES.PATHR  Path consistency R: the count-weighted mean resultant length over tile-senses with n >= nMin (nxr mean_R).
%
%   R = rheome.cores.pathR(sums, n)          % nMin = 5 (nxr NMIN)
%   R = rheome.cores.pathR(sums, n, nMin)
%
% sums: the tile-senses' summed crossings, tangent vectors at the tile centres [2T x 3] (rheome.cores.crossings) or the same as complex
% numbers in the centres' gauges [2T x 1]. R = sum |sums| / sum n over the rows (tile-senses) with n >= nMin; NaN when no row qualifies. The senses
% are not pooled. Ports nxr mean_R (vortex_cores.py).
%
% See also: rheome.cores.crossings, rheome.cores.pathchance
%
% Author: Diellor Basha, 2026
    if nargin < 3, nMin = 5; end
    n = n(:);  m = n >= nMin;
    if ~any(m), R = NaN; else, R = sum(vecnorm(sums(m, :), 2, 2)) / sum(n(m)); end   % vecnorm of a complex column is its abs
end
