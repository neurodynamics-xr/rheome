function [sums, n, key, vec] = crossings(S, frames, lab, V, normals, support)
% CORES.CROSSINGS  Per tile-sense, the summed crossing directions of the core steps at one ladder level (nxr Tiles.crossings).
%
%   [sums, n, key, vec] = rheome.cores.crossings(S, frames, lab, V, normals, support)
%
% S [s x 9] core steps (rheome.cores.stats rows: sense, from xyz, displacement xyz, track, frame); cast to single as nxr does. frames: the
% frames kept (ismember on the frame column). lab [nV x 1] tile codes of ONE level (a column of rheome.cores.ladder), V [nV x 3] mm,
% normals [nV x 3] outward unit, support: 1-based vertices cores may snap to. Each step goes to the Euclidean nearest SUPPORT vertex of its
% origin (3-D, not geodesic); its tile is lab there and its displacement is projected on the tangent plane at that vertex's normal (not
% renormalised first); steps of length <= 1e-9 are dropped, the rest become unit vectors. Row = 2*tile + (sense > 0) + 1, nt = 2*(max(lab)+1).
% A crossing = one (row, track) pair: that track's unit vectors summed, then normalised (norm floored at 1e-12). Crossings are ordered by
% (row, track) lexicographically; key [c x 1] is the row of each, vec [c x 3] its unit vector, sums [nt x 3] their sum per row, n [nt x 1]
% their count per row. Ties in nearest distance are broken by the search tree (avoid). Ports nxr-analysis src/nxr/vortex_cores.py
% Tiles.crossings (line 61, commit 0db82e4).
%
% See also: rheome.cores.ladder, rheome.cores.pathR, rheome.cores.pathchance
%
% Author: Diellor Basha, 2026
    S = double(single(S));
    S = S(ismember(round(S(:, 9)), frames), :);
    nt = 2 * (double(max(lab)) + 1);
    Vs = double(V(support, :));  Ns = double(normals(support, :));
    if isempty(S), sums = zeros(nt, 3); n = zeros(nt, 1); key = zeros(0, 1); vec = zeros(0, 3); return; end
    v = knnsearch(Vs, S(:, 2:4));
    nrm = Ns(v, :);
    d = S(:, 5:7) - sum(S(:, 5:7) .* nrm, 2) .* nrm;
    ln = vecnorm(d, 2, 2);  ok = ln > 1e-9;
    tile = double(lab(support));
    row = tile(v) * 2 + (S(:, 1) > 0) + 1;
    [g, ~, inv] = unique([row(ok) round(S(ok, 8))], 'rows');
    u = d(ok, :) ./ ln(ok);
    vec = [accumarray(inv, u(:,1)), accumarray(inv, u(:,2)), accumarray(inv, u(:,3))];
    vec = vec ./ max(vecnorm(vec, 2, 2), 1e-12);
    key = g(:, 1);
    sums = [accumarray(key, vec(:,1), [nt 1]), accumarray(key, vec(:,2), [nt 1]), accumarray(key, vec(:,3), [nt 1])];
    n = accumarray(key, 1, [nt 1]);
end
