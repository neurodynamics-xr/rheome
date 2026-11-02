function R = pathR(sums, n, nMin)
% CORES.PATHR  Path consistency R: the count-weighted mean resultant length over tile-senses with n >= nMin (nxr mean_R).
%
%   R = rheome.cores.pathR(sums, n)          % nMin = 5 (nxr NMIN)
%   R = rheome.cores.pathR(sums, n, nMin)
%
% R = sum |sums| / sum n over the rows (tile-senses) with n >= nMin; NaN when no row qualifies. The senses are not pooled. Ports nxr-analysis
% src/nxr/vortex_cores.py mean_R (line 75, commit 0db82e4).
%
% See also: rheome.cores.crossings, rheome.cores.pathchance
%
% Author: Diellor Basha, 2026
    if nargin < 3, nMin = 5; end
    n = n(:);  m = n >= nMin;
    if ~any(m), R = NaN; else, R = sum(vecnorm(sums(m, :), 2, 2)) / sum(n(m)); end
end
