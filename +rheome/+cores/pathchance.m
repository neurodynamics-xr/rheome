function [chance, R] = pathchance(key, vec, nt, perms, nMin)
% CORES.PATHCHANCE  The chance level of path consistency: R after shuffling the crossing directions among the crossings (nxr stage_measures).
%
%   [chance, R] = rheome.cores.pathchance(key, vec, nt, perms)               % perms: cell of permutations (1-based), replayed
%   [chance, R] = rheome.cores.pathchance(key, vec, nt, nshuf)               % or nshuf draws from the default stream
%   [chance, R] = rheome.cores.pathchance(key, vec, nt, nshuf, nMin)         % nMin default 5
%
% key, vec [c x 1]: the crossings in rheome.cores.crossings order, vec complex (each in its own tile's gauge). For each permutation p the shuffled
% numbers vec(p) are re-added to the SAME rows (the counts are kept; each crossing takes another's angle in its own gauge) and rheome.cores.pathR is taken; chance is the mean of those, R [nshuf x 1] each one. nxr draws
% rng.permutation(len(u)) NSHUF = 5 times per level and stage; numpy's generator cannot be reproduced, so the permutations are an input
% (the nxr goldens store them) or come from MATLAB's generator. Ports nxr-analysis src/nxr/vortex_cores.py stage_measures / _sum
% (nxr 1.15.0; commit 0db82e4 added 3-D vectors).
%
% See also: rheome.cores.pathR, rheome.cores.crossings
%
% Author: Diellor Basha, 2026
    if nargin < 5, nMin = 5; end
    if isnumeric(perms) && isscalar(perms)
        c = cell(perms, 1);
        for i = 1:perms, c{i} = randperm(numel(vec)); end
        perms = c;
    end
    R = zeros(numel(perms), 1);  n = accumarray(key, 1, [nt 1]);
    for i = 1:numel(perms)
        u = vec(perms{i}(:));
        s = accumarray(key, real(u), [nt 1]) + 1i * accumarray(key, imag(u), [nt 1]);
        R(i) = rheome.cores.pathR(s, n, nMin);
    end
    chance = mean(R);
end
