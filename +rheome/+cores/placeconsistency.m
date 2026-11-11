function [R, chance, nUsed] = placeconsistency(S, pool, support, V, normals, places, perms, nMin)
% CORES.PLACECONSISTENCY  Direction consistency at places: every step a unit tangent vector CARRIED to its place, no per-track averaging (nxr direction_consistency).
%
%   [R, chance, nUsed] = rheome.cores.placeconsistency(S, pool, support, V, normals, places, perms)           % nMin = 10
%   [R, chance, nUsed] = rheome.cores.placeconsistency(S, pool, support, V, normals, places, nshuf, nMin)
%
% S [s x 9] core steps (sense, from xyz, displacement xyz, track, frame; double, not cast). pool = rheome.cores.placepool(kit, Vfull, support, places);
% support [nk x 1] the 1-based vertices of the full surface; V [nk x 3] and normals [nk x 3] those of the SUPPORT; places: indices INTO the support.
% Each step goes to the support vertex nearest its origin, is projected on the tangent plane there (not renormalised first), dropped if the length
% is <= 1e-9, else made a unit vector, CARRIED to the place its vertex is pooled at (the vector heat method's transport; never added in 3-D across
% a fold) and added there. key = place position*2 + (sense > 0). R = sum |sums| / sum n over keys with n >= nMin; chance = the mean over the
% shuffles of the same with the carried vectors of ALL usable steps permuted among the keys (one permutation of the usable steps per shuffle, in
% the original step order). perms: a cell of 1-based permutations (replayed) or nshuf (default 20, MATLAB's generator). nUsed = usable steps.
% NaN, NaN, 0 when S is empty. Ports nxr direction_consistency (nxr 1.15.0; before, steps went to the place nearest their origin and were added in 3-D).
%
% See also: rheome.cores.placepool, rheome.cores.crossings, rheome.cores.pathR
%
% Author: Diellor Basha, 2026
    if nargin < 8, nMin = 10; end
    if nargin < 7 || isempty(perms), perms = 20; end
    R = NaN;  chance = NaN;  nUsed = 0;
    if isempty(S), return; end
    S = double(S);
    near = knnsearch(double(V), S(:, 2:4));
    nrm = double(normals(near, :));
    d = S(:, 5:7) - sum(S(:, 5:7) .* nrm, 2) .* nrm;
    ln = vecnorm(d, 2, 2);  ok = ln > 1e-9;
    u = pool.carry(d(ok, :) ./ ln(ok), support(near(ok)));
    key = (pool.PlaceOf(support(near(ok))) - 1) * 2 + (S(ok, 1) > 0) + 1;
    nk = max(key);
    nUsed = nnz(ok);
    R = i_meanR(key, u, nk, nMin);
    if isnumeric(perms) && isscalar(perms)
        c = cell(perms, 1);
        for i = 1:perms, c{i} = randperm(nUsed); end
        perms = c;
    end
    chance = mean(cellfun(@(p) i_meanR(key, u(p(:)), nk, nMin), perms));
end

function R = i_meanR(key, U, nk, nMin)
    s = accumarray(key, real(U), [nk 1]) + 1i * accumarray(key, imag(U), [nk 1]);
    R = rheome.cores.pathR(s, accumarray(key, 1, [nk 1]), nMin);
end
