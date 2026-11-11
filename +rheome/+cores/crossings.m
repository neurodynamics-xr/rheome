function [sums, n, key, vec] = crossings(S, frames, lab, V, normals, support, pool)
% CORES.CROSSINGS  Per tile-sense, the summed crossing directions of the core steps at one ladder level, CARRIED to the tile's centre (nxr Tiles.crossings).
%
%   [sums, n, key, vec] = rheome.cores.crossings(S, frames, lab, V, normals, support, pool)
%
% S [s x 9] core steps (rheome.cores.stats rows: sense, from xyz, displacement xyz, track, frame); cast to single as nxr does. frames: the
% frames kept (ismember on the frame column). lab [nV x 1] tile codes of ONE level (a column of rheome.cores.ladder), V [nV x 3] mm,
% normals [nV x 3] outward unit, support: 1-based vertices cores may snap to, pool = rheome.cores.tilepool(kit, V, lab, support).
% Each step goes to the Euclidean nearest SUPPORT vertex of its origin (3-D, not geodesic); its tile is lab there and its displacement is
% projected on the tangent plane at that vertex's normal (not renormalised first); steps of length <= 1e-9 (or off the tiles) are dropped,
% the rest become unit vectors, CARRIED to their tile's centre (the vector heat method's transport: vectors at different places of a folded
% tile are never added in 3-D). Row = 2*tile + (sense > 0) + 1, nt = 2*(max(lab)+1).
% A crossing = one (row, track) pair: that track's carried vectors summed, then made a unit complex number (norm floored at 1e-12, in the
% centre's gauge, rheome.geom.GeometryKit). Crossings are ordered by (row, track) lexicographically; key [c x 1] is the row of each,
% vec [c x 1] complex its unit direction, n [nt x 1] their count per row, and sums [nt x 3] their sum per row as a tangent vector AT THE
% CENTRE (zeros for an empty tile). Ties in nearest distance are broken by the search tree (avoid). Ports nxr Tiles.crossings (nxr 1.15.0;
% before, the vectors were added in 3-D: nxr-analysis src/nxr/vortex_cores.py line 61, commit 0db82e4).
%
% See also: rheome.cores.tilepool, rheome.cores.ladder, rheome.cores.pathR, rheome.cores.pathchance
%
% Author: Diellor Basha, 2026
    S = double(single(S));
    S = S(ismember(round(S(:, 9)), frames), :);
    nt = 2 * (double(max(lab)) + 1);
    Vs = double(V(support, :));  Ns = double(normals(support, :));
    centre = zeros(nt / 2, 1);  centre(1:numel(pool.Centre)) = pool.Centre;
    if isempty(S), sums = zeros(nt, 3); n = zeros(nt, 1); key = zeros(0, 1); vec = zeros(0, 1); return; end
    v = knnsearch(Vs, S(:, 2:4));
    nrm = Ns(v, :);
    d = S(:, 5:7) - sum(S(:, 5:7) .* nrm, 2) .* nrm;
    ln = vecnorm(d, 2, 2);
    tile = double(lab(support));  tv = tile(v);
    ok = ln > 1e-9 & tv >= 0;
    row = tv * 2 + (S(:, 1) > 0) + 1;
    w = pool.carry(d(ok, :) ./ ln(ok), support(v(ok)));
    [g, ~, inv] = unique([row(ok) round(S(ok, 8))], 'rows');
    s = accumarray(inv, real(w)) + 1i * accumarray(inv, imag(w));
    vec = s ./ max(abs(s), 1e-12);
    key = g(:, 1);
    sc = accumarray(key, real(vec), [nt 1]) + 1i * accumarray(key, imag(vec), [nt 1]);
    n = accumarray(key, 1, [nt 1]);
    sums = pool.embed(sc, reshape(repmat(centre(:)', 2, 1), [], 1));
end
