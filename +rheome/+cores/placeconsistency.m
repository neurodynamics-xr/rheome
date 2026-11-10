function [R, chance, nUsed] = placeconsistency(S, V, normals, places, perms, nMin)
% CORES.PLACECONSISTENCY  Direction consistency at places: every step as a unit tangent vector, no per-track averaging (nxr direction_consistency).
%
%   [R, chance, nUsed] = rheome.cores.placeconsistency(S, V, normals, places, perms)           % nMin = 10
%   [R, chance, nUsed] = rheome.cores.placeconsistency(S, V, normals, places, nshuf, nMin)
%
% S [s x 9] core steps (sense, from xyz, displacement xyz, track, frame; double, not cast). Each step goes to the place nearest its origin
% among `places` (1-based vertices of V / normals, not of a support set), is projected on the tangent plane there (not renormalised first),
% dropped if the length is <= 1e-9, else made a unit vector. key = place position*2 + (sense > 0). R = sum |sums| / sum n over keys with
% n >= nMin; chance = the mean over the shuffles of the same with the unit vectors of ALL usable steps permuted among the keys (one
% permutation of the usable steps per shuffle, in the original step order). perms: a cell of 1-based permutations (replayed) or nshuf
% (default 20, MATLAB's generator). nUsed = usable steps. NaN, NaN, 0 when S is empty. Ports nxr-analysis src/nxr/vortex.py
% direction_consistency (line 760, commit 0db82e4).
%
% See also: rheome.cores.crossings, rheome.cores.pathR
%
% Author: Diellor Basha, 2026
    if nargin < 6, nMin = 10; end
    if nargin < 5 || isempty(perms), perms = 20; end
    R = NaN;  chance = NaN;  nUsed = 0;
    if isempty(S), return; end
    S = double(S);  P = double(places(:));
    near = knnsearch(double(V(P, :)), S(:, 2:4));
    nrm = double(normals(P(near), :));
    d = S(:, 5:7) - sum(S(:, 5:7) .* nrm, 2) .* nrm;
    ln = vecnorm(d, 2, 2);  ok = ln > 1e-9;
    u = d(ok, :) ./ ln(ok);
    key = near(ok) * 2 + (S(ok, 1) > 0) - 1;                 % (0-based position)*2 + sense, plus 1
    nk = max(key);
    nUsed = nnz(ok);
    R = i_meanR(key, u, nk, nMin);
    if isnumeric(perms) && isscalar(perms)
        c = cell(perms, 1);
        for i = 1:perms, c{i} = randperm(nUsed); end
        perms = c;
    end
    chance = mean(cellfun(@(p) i_meanR(key, u(p, :), nk, nMin), perms));
end

function R = i_meanR(key, U, nk, nMin)
    s = [accumarray(key, U(:,1), [nk 1]), accumarray(key, U(:,2), [nk 1]), accumarray(key, U(:,3), [nk 1])];
    R = rheome.cores.pathR(s, accumarray(key, 1, [nk 1]), nMin);
end
