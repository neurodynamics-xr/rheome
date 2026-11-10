function labels = ladder(V, F, hemi)
% CORES.LADDER  The surface ladder of the whole two-piece cortex, tile codes 2^level*h + t, one connected piece per tile (nxr surface_ladder).
%
%   labels = rheome.cores.ladder(V, F, hemi)
%
% V [n x 3] positions in mm, F [m x 3] faces (1-based), hemi [n x 1]: 0 / 'L' = left, 1 / 'R' = right (nxr hemispheres(V, F, '+y'):
% the two largest components, L the one with the larger mean y). Each hemisphere is laddered alone on the sub-mesh of the faces whose three
% vertices it holds. A tile is halved across its geodesic diameter (Dijkstra on the edge graph, Euclidean edge lengths; the distances
% from the tile's first vertex give the far end a, from a the far end b; vertices STABLY sorted by d_a - d_b, 0 where either is
% unreachable), cut where the running vertex mass (a third of each face) first reaches half, searchsorted 'left' + 1, clipped to 1..n-1;
% a tile of < 2 vertices goes whole to child 2t. Then each child keeps its heaviest connected piece (parity 0, then parity 1; the other
% pieces go to the sibling, code xor 1; ties by first vertex). The ladder stops at the first level that leaves a tile empty, or past
% 2^floor(log2 n) tiles. OUTPUT: labels [n x L] int32, column l+1 = level l (the levels BOTH hemispheres reach), code = 2^l*h + t, 0-based,
% -1 off both hemispheres. A parent is bitshift(code, -1) within the hemisphere.
%
% Differs from rheome.geom.surfaceladder (unchanged), which ports nsp's ladder: that one takes one connected surface, has no one-piece
% step, no hemisphere code, and treats a tile whose first vertex is out of a's reach differently (see CHANGELOG). Ports nxr-analysis
% src/nxr/tile_connectome.py surface_ladder / ladder / one_piece (commit 0db82e4).
%
% See also: rheome.cores.crossings, rheome.geom.surfaceladder
%
% Author: Diellor Basha, 2026
    arguments
        V (:,3) double
        F (:,3) double
        hemi
    end
    n = size(V, 1);
    if isstring(hemi) || ischar(hemi) || iscellstr(hemi), hemi = double(string(hemi(:)) == "R"); end
    hemi = double(hemi(:));
    per = cell(1, 2);
    for h = 0:1
        idx = find(hemi == h);
        remap = zeros(n, 1);  remap(idx) = 1:numel(idx);
        keep = all(remap(F) > 0, 2);
        per{h + 1} = {idx, i_ladder(V(idx, :), remap(F(keep, :)))};
    end
    L = min(size(per{1}{2}, 2), size(per{2}{2}, 2));
    labels = -ones(n, L, 'int32');
    for h = 0:1
        idx = per{h + 1}{1};  lab = per{h + 1}{2};
        for l = 0:L-1, labels(idx, l + 1) = int32(2^l * h + lab(:, l + 1)); end
    end
end

function labels = i_ladder(V, F)
    n = size(V, 1);
    a = i_norm(cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2)) / 2;
    mass = accumarray(reshape(F.', [], 1), repelem(a / 3, 3), [n 1]);
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    w = i_norm(V(E(:,1),:) - V(E(:,2),:));
    labels = zeros(n, 1);  tiles = 1;
    while 2 ^ size(labels, 2) <= 2 ^ floor(log2(n))
        prev = labels(:, end);
        counts = accumarray(prev + 1, 1, [tiles 1]);
        first = zeros(tiles, 1);  [pr, ix] = unique(prev);  first(pr + 1) = ix;
        inner = prev(E(:,1)) == prev(E(:,2));
        Ei = E(inner, :);  wi = w(inner);
        d0 = i_sweep(Ei, wi, n, first);
        a0 = i_farthest(prev, d0, tiles);
        da = i_sweep(Ei, wi, n, a0);
        b0 = i_farthest(prev, da, tiles);
        db = i_sweep(Ei, wi, n, b0);
        key = zeros(n, 1);  fin = isfinite(da) & isfinite(db);  key(fin) = da(fin) - db(fin);
        nxt = zeros(n, 1);
        small = counts(prev + 1) < 2;
        nxt(small) = 2 * prev(small);
        [~, o] = sortrows([prev key]);                      % stable: ties by vertex index
        ends = cumsum(counts);
        for t = find(counts >= 2).' - 1
            ord = o(ends(t + 1) - counts(t + 1) + 1 : ends(t + 1));
            acc = cumsum(mass(ord));  m = numel(ord);
            cut = min(max(find(acc >= acc(end) / 2, 1), 1), m - 1);
            nxt(ord(1:cut)) = 2 * t;  nxt(ord(cut+1:end)) = 2 * t + 1;
        end
        nxt = i_onepiece(i_onepiece(nxt, mass, E, 0), mass, E, 1);
        if any(accumarray(nxt + 1, 1, [2 * tiles 1]) == 0), break; end
        labels(:, end + 1) = nxt;  %#ok<AGROW>
        tiles = tiles * 2;
    end
end

% Children of `parity` keep their heaviest connected piece; the others go to the sibling. Edges count only inside a tile.
function out = i_onepiece(lab, mass, E, parity)
    keep = lab(E(:,1)) == lab(E(:,2));
    comp = conncomp(graph(E(keep,1), E(keep,2), [], numel(lab))).';
    [~, ~, comp] = unique(comp, 'stable');                  % components numbered by first vertex
    cm = accumarray(comp, mass);
    tileOf = accumarray(comp, lab, [], @max);
    [~, ord] = sort(-cm);                                   % stable: ties by component label
    best = -ones(max(lab) + 1, 1);
    for k = ord.'
        if best(tileOf(k) + 1) == -1, best(tileOf(k) + 1) = k; end
    end
    out = lab;
    move = (bitand(lab, 1) == parity) & (best(lab + 1) ~= comp);
    out(move) = bitxor(lab(move), 1);
end

% Per tile, the first vertex (vertex order) at the tile's largest finite value.
function out = i_farthest(lab, v, tiles)
    v(~isfinite(v)) = -1;
    best = accumarray(lab + 1, v, [tiles 1], @max, -Inf);
    hit = find(v == best(lab + 1));
    out = zeros(tiles, 1);
    [t, f] = unique(lab(hit));
    out(t + 1) = hit(f);
end

% Multi-source Dijkstra on within-tile edges: one source per tile (0 = none) via a super source at zero weight.
function d = i_sweep(E, w, n, src)
    src = src(src > 0);
    G = graph([E(:,1); repmat(n + 1, numel(src), 1)], [E(:,2); src], [w; zeros(numel(src), 1)], n + 1);
    d = distances(G, n + 1, 'Method', 'positive').';
    d = d(1:n);
end

% Row norm as numpy's np.linalg.norm(axis=1): sqrt of the sequential sum of squares (vecnorm may round differently; ties in the bisection key
% da - db are decided by the last bit of the edge lengths).
function h = i_norm(D)
    h = sqrt(D(:,1).^2 + D(:,2).^2 + D(:,3).^2);
end
