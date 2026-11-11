function [s, n] = senses(sums, n)
% CORES.SENSES  Pool the two senses of each tile: rows 2k-1 and 2k of the tile-sense table are added (nxr senses).
%
%   [s, n] = rheome.cores.senses(sums, n)
%
% sums [2T x 3] (tangent vectors at the tile centres: both senses of a tile share one centre and one plane, so they add), n [2T x 1] from rheome.cores.crossings; output [T x 3], [T x 1], tile k = code k-1. Ports nxr-analysis
% src/nxr/vortex_cores.py senses (line 79, commit 0db82e4).
%
% Author: Diellor Basha, 2026
    s = sums(1:2:end, :) + sums(2:2:end, :);
    n = n(1:2:end) + n(2:2:end);
end
