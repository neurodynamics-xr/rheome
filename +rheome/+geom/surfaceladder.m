function [labels, mass] = surfaceladder(V, F, maxTiles)
% GEOM.SURFACELADDER  The atlas's dyadic surface ladder: equal-mass geodesic bisection, bit-exact to nsp.
%
%   [labels, mass] = rheome.geom.surfaceladder(V, F)
%   [labels, mass] = rheome.geom.surfaceladder(V, F, maxTiles)
%
% A port, step for step, of nsp's atlas ladder (nxr_convert/atlas/ladder.py, itself a port of
% cortical-flow's surfaceLadder), so labels(:, l+1) is the stored atlas code of every vertex at
% level l. Each level halves every tile by EQUAL MASS across its geodesic diameter: two Dijkstra
% sweeps inside the tile find the far ends a and b, vertices are STABLY sorted by d_a - d_b and cut
% where the running mass first reaches half the tile's. Level l has 2^l tiles; tile t has children
% 2t and 2t+1, so a vertex's level-l tile is floor(code_L / 2^(L-l)). The ladder stops before a
% level that would leave a tile empty.
%
% ⚠ BIT-EXACT MEANS THE SAME FLOATING-POINT ORDER. Mass is accumulated face by face (a, b, c), edge
% lengths use V8's Math.hypot (scaled, Kahan-summed), tile masses are summed in vertex order: each
% is written as the reference does it. Dijkstra distances are the fixed point min_u fl(d_u + w_uv),
% so any exact Dijkstra (MATLAB's distances, scipy's) returns the same bits.
%
% INPUTS   V [n x 3] positions, F [m x 3] faces (1-based), maxTiles (default 2^floor(log2 n))
% OUTPUTS  labels [n x (finest+1)] int32 tile codes, 0-based, column l+1 = level l;
%          mass [n x 1] barycentric lumped vertex areas (a third of each incident triangle)
%
% See also: rheome.geom.tree (Method="Ladder"), tests/tSurfaceLadder.m
%
% Author: Diellor Basha, 2026

    V = double(V);  F = double(F);  n = size(V, 1);
    if n < 1 || isempty(F), error('rheome:geom:surfaceladder', 'a mesh of %d vertices and %d faces', n, size(F, 1)); end
    cap = 2 ^ floor(log2(n));
    if nargin < 3 || isempty(maxTiles), maxTiles = cap; end
    maxTiles = min(maxTiles, cap);

    % mass: face by face, a then b then c (the reference's sequential np.add.at)
    ar = 0.5 * i_hypot(cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2)) / 3;
    mass = zeros(n, 1);  fv = reshape(F', [], 1);  av = repelem(ar, 3);
    for k = 1:numel(fv), mass(fv(k)) = mass(fv(k)) + av(k); end

    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    w = i_hypot(V(E(:,1),:) - V(E(:,2),:));

    labels = zeros(n, 1, 'int32');  tiles = 1;
    while tiles * 2 <= maxTiles
        prev = double(labels(:, end));
        counts = accumarray(prev + 1, 1, [tiles 1]);
        first = zeros(tiles, 1);  [pr, ix] = unique(prev);  first(pr + 1) = ix;   % 0 = absent
        inner = prev(E(:,1)) == prev(E(:,2));
        Ei = E(inner, :);  wi = w(inner);
        d0 = i_sweep(Ei, wi, n, first);
        a = i_firstargmax(prev, d0, tiles);  a(a == 0) = first(a == 0);
        da = i_sweep(Ei, wi, n, a);
        b = i_firstargmax(prev, da, tiles);  b(b == 0) = first(b == 0);
        stuck = first > 0 & ~isfinite(da(max(first, 1)));   % first vertex out of a's reach stays the far end
        b(stuck) = first(stuck);
        db = i_sweep(Ei, wi, n, b);
        key = da - db;                                       % NaN where both ends are out of reach

        nxt = zeros(n, 1);
        small = counts(prev + 1) < 2;
        nxt(small) = 2 * prev(small);
        [~, byTile] = sort(prev);  ends = cumsum(counts);    % stable: each tile's vertices in vertex order
        for t = find(counts >= 2)' - 1
            mem = byTile(ends(t + 1) - counts(t + 1) + 1 : ends(t + 1));
            [order, total] = i_order(mem, key(mem), mass(mem));
            acc = cumsum(mass(order));
            m = numel(order);
            cut = find(acc >= total / 2, 1);
            if isempty(cut), cut = m - 1; end
            cut = min(max(cut, 1), m - 1);
            nxt(order(1:cut)) = 2 * t;  nxt(order(cut+1:end)) = 2 * t + 1;
        end
        if any(accumarray(nxt + 1, 1, [2 * tiles 1]) == 0), break; end
        labels(:, end + 1) = int32(nxt); %#ok<AGROW>
        tiles = tiles * 2;
    end
end

% A tile's members in bisection order, and its mass summed in vertex order.
function [order, total] = i_order(mem, k, ms)
    acc = cumsum(ms);  total = acc(end);
    if all(isfinite(k))
        [~, o] = sort(k);                                    % MATLAB's sort is stable
        order = mem(o);
    else
        % the comparator may return NaN (read as 0): the reference's own sort, step for step
        if numel(k) < 64, order = mem(i_v8sort(k)); else, order = mem(rheome.utils.pysort(k)); end
    end
end

% Array.prototype.sort as V8 runs it (below 64 elements) on 1:numel(k) with compare (x,y) -> k(x)-k(y), NaN read as 0.
function a = i_v8sort(k)
    n = numel(k);  a = (1:n)';
    cmp = @(x, y) i_nan0(k(x) - k(y));
    if n < 2, return; end
    run = 1;
    if n >= 8
        run = 2;  desc = cmp(a(2), a(1)) < 0;  p = a(2);
        for i = 3:n
            o = cmp(a(i), p);
            if (desc && o >= 0) || (~desc && o < 0), break; end
            p = a(i);  run = run + 1;
        end
        if desc, a(1:run) = a(run:-1:1); end
    end
    for s = run+1:n
        pv = a(s);  lo = 1;  hi = s;                         % 1-based: insert into a(1:s-1)
        while lo < hi
            mid = lo + floor((hi - lo) / 2);
            if cmp(pv, a(mid)) < 0, hi = mid; else, lo = mid + 1; end
        end
        a(lo+1:s) = a(lo:s-1);  a(lo) = pv;
    end
end

function c = i_nan0(c)
    if isnan(c), c = 0; end
end

% Per tile, the first vertex (vertex order) at the tile's largest finite value; 0 for none.
function out = i_firstargmax(lab, v, tiles)
    v(~isfinite(v)) = -Inf;
    best = accumarray(lab + 1, v, [tiles 1], @max, -Inf);
    hit = find(v == best(lab + 1) & isfinite(v));
    out = zeros(tiles, 1);
    [t, f] = unique(lab(hit));            % first occurrence
    out(t + 1) = hit(f);
end

% Multi-source Dijkstra on the within-tile edges: one source per tile (0 = none), via a super
% source at zero weight (fl(0 + w) = w, so the bits are those of the per-tile sweeps).
function d = i_sweep(E, w, n, src)
    src = src(src > 0);
    G = graph([E(:,1); repmat(n + 1, numel(src), 1)], [E(:,2); src], [w; zeros(numel(src), 1)], n + 1);
    d = distances(G, n + 1, 'Method', 'positive')';
    d = d(1:n);
end

% Row-wise Math.hypot as V8 computes it: scale by the largest |x|, Kahan-sum the squares, sqrt, rescale.
function h = i_hypot(D)
    x = abs(D);  mx = max(x, [], 2);
    sf = mx;  sf(mx <= 0) = 1;
    s = zeros(size(x, 1), 1);  c = s;
    for k = 1:size(x, 2)
        q = x(:, k) ./ sf;
        summand = q .* q - c;  pre = s + summand;
        c = (pre - s) - summand;  s = pre;
    end
    h = sqrt(s) .* mx;  h(mx <= 0) = 0;
end
