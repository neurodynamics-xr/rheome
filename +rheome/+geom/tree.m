function T = tree(S, opts)
% GEOM.TREE  Recursive spectral bisection of a surface: the cortical position pyramid.
%
%   T = rheome.geom.tree(S)
%   T = rheome.geom.tree(S, MinArea=2e-3, MaxDepth=6)
%   T = rheome.geom.tree(S, L=L, M=M)                 % reuse a built Laplace-Beltrami pair
%   B = rheome.load.bases(name);  T = rheome.geom.tree(B.L.S, L=B.L.lbo.L, M=B.L.lbo.M)   % from the cache
%
% The cortical twin of rheome.sensors.tree: split each patch by the sign of its own Fiedler
% function, down to patches of a chosen size. A NESTED, DISJOINT partition, which is the
% property that matters -- a parent's members are exactly its two children's, so sums over
% vertices roll up the tree the way sums over samples roll up the dyadic time grid.
%
% ⭐ WHY THIS IS NOT A DIFFUSION SCALE-SPACE, AND WHY YOU WANT BOTH. Diffused deltas form a
% partition of UNITY: smooth, overlapping, scale-covariant, and not additive -- a parent
% window is not the sum of child windows, so nothing merges. This tree is a partition: no
% smoothness, exact merging. Use the tree for bookkeeping and for apertures you want to
% compare across sizes; use a diffusion bank (rheome.graphfilterbank.fromOperator) for measurement.
%
% ⚠ THE FIEDLER VECTOR HERE IS THE SURFACE'S, NOT SPACE'S. Two banks of a sulcus are
% millimetres apart in the volume and far apart on the cortex, so a split on coordinates
% would cut across the sheet. The second eigenfunction of the patch's own Laplace-Beltrami
% operator respects the geodesic structure and cuts along it.
%
% ⭐ TAKE L AND M FROM THE CACHE, NOT FROM A REBUILD. rheome.load.bases holds the surface and the
% solved Laplace-Beltrami pair per hemisphere (B.L.S, B.L.lbo.L, B.L.lbo.M), so nothing here
% needs rheome.operators.laplace_beltrami again.
%
% ⚠ BUT THE CACHED EIGENVECTORS ONLY COVER THE FIRST SPLIT. The cached Phi(:,2) IS the first
% cut -- measured on the reference left hemisphere, the sign of the cached Fiedler function and
% this function's own root split agree on 1.000 of vertices -- and below that they diverge,
% because a patch's Fiedler function is an eigenfunction of the RESTRICTED operator with a
% free boundary on the cut, which no global basis contains. The per-patch solves are real
% work and are kept. They are also cheap: one hemisphere to depth 5 (63 nodes) takes 0.65 s,
% because every solve is sparse on a shrinking submatrix.
%
% ⚠ THE CORTEX IS TWO COMPONENTS. Its Fiedler value is zero with the two hemispheres on
% opposite signs, which is the right first cut but a degenerate eigenproblem; a node spanning
% more than one connected component is therefore split by component before any eigensolve.
%
% ⚠ AREA, NOT EXTENT. A patch's size is the area its vertices carry in the mass matrix, and
% its diameter is the equivalent disc's, 2*sqrt(area/pi). Euclidean extent across a folded
% sheet measures the folding, not the patch.
%
% INPUTS
%   S        surface struct (.Vertices, .Faces) or a struct with .L and .M supplied
%   L, M     precomputed rheome.operators.laplace_beltrami output (optional)
%   MinVertices, MinArea, MaxDepth   stopping rules (1, 0, Inf); any one of them stops a node
%
% OUTPUT (table T, root first, children after parents):
%   node_id, parent_id, depth, n_vertices, is_leaf, members {row indices},
%   area m^2, diameter m (equivalent disc), centroid [x y z] m, wavelength_mm (the diameter
%   as a wavelength, for comparison with the instrument's resolution floor)
%
% See also: rheome.sensors.tree, rheome.ingest.groups, rheome.detect.blobscale, rheome.operators.laplace_beltrami
%
% Author: Diellor Basha, 2026

    arguments
        S (1,1) struct
        opts.L = []
        opts.M = []
        opts.MinVertices (1,1) double {mustBeInteger, mustBePositive} = 1
        opts.MinArea     (1,1) double {mustBeNonnegative} = 0
        opts.MaxDepth    (1,1) double = Inf
        opts.Verbose     (1,1) logical = false
    end
    L = opts.L;  M = opts.M;
    if isempty(L) || isempty(M)
        [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);
    end
    P = S.Vertices;  n = size(P, 1);
    a = full(sum(M, 2));                                      % vertex areas (row sums: exact)
    A = abs(L) > 0;  A = A - diag(diag(A));                   % adjacency, for components
    g = graph(A ~= 0);

    rows = {};  queue = {struct('idx', (1:n)', 'parent', 0, 'depth', 0)};  nextId = 1;
    while ~isempty(queue)
        nd = queue{1};  queue(1) = [];
        idx = nd.idx;  id = nextId;  nextId = nextId + 1;
        ar = sum(a(idx));
        stop = numel(idx) <= opts.MinVertices || nd.depth >= opts.MaxDepth || ...
               (opts.MinArea > 0 && ar <= opts.MinArea);
        isLeaf = true;
        if ~stop
            [p, q] = i_split(L, M, g, idx, P);
            if ~isempty(p) && ~isempty(q)
                isLeaf = false;
                queue{end+1} = struct('idx', p, 'parent', id, 'depth', nd.depth + 1); %#ok<AGROW>
                queue{end+1} = struct('idx', q, 'parent', id, 'depth', nd.depth + 1); %#ok<AGROW>
            end
        end
        d = 2 * sqrt(ar / pi);
        rows(end+1, :) = {id, nd.parent, nd.depth, numel(idx), isLeaf, idx(:)', ar, d, ...
                          mean(P(idx, :), 1), 1e3 * d};                          %#ok<AGROW>
        if opts.Verbose && mod(id, 8) == 1
            fprintf('  node %d: depth %d, %d vertices, %.1f cm^2\n', id, nd.depth, numel(idx), 1e4*ar);
        end
    end
    T = cell2table(rows, 'VariableNames', {'node_id','parent_id','depth','n_vertices', ...
        'is_leaf','members','area','diameter','centroid','wavelength_mm'});
end

% Components first, then the Fiedler function of the patch's own operator. ⚠ eigs, never eig:
% a dense eigendecomposition of a 20484-vertex patch is not an option, which is the one thing
% that stops rheome.sensors.tree being reused here verbatim.
function [p, q] = i_split(L, M, g, idx, P)
    c = conncomp(subgraph(g, idx));
    if max(c) > 1
        [p, q] = i_bycomponent(L, M, g, idx, c, P);
        return
    end
    Ls = L(idx, idx);  Ms = M(idx, idx);
    Ls = (Ls + Ls') / 2;  Ms = (Ms + Ms') / 2;
    f = [];
    try
        wr = warning('off', 'all');
        [V, D] = eigs(Ls, Ms, 2, 'smallestabs', 'Tolerance', 1e-6, 'MaxIterations', 500);
        warning(wr);
        [~, o] = sort(diag(D));
        f = V(:, o(2));
    catch
        f = [];
    end
    if isempty(f) || all(f >= 0) || all(f < 0)
        [p, q] = i_bylongest(idx, P);                          % degenerate: fall back
        return
    end
    p = idx(f >= 0);  q = idx(f < 0);
end

% Split a multi-component patch into two groups of whole components, balanced by area.
% ⚠ A DOMINANT COMPONENT IS CUT, NOT KEPT WHOLE. A Fiedler cut can leave a few vertices joined to
% their patch only through vertices outside it; the next split then sees them as components, and
% whole-component balancing puts the main body on one side and a one-vertex fragment on the other --
% a 1-vertex node that cannot split again, so the tree stops being full (another subject's left hemisphere, node 23 at
% depth 4, 2.3 mm^2). When one component holds more than 2/3 of the area it is split by its own
% Fiedler vector and the fragments go, largest first, to the lighter side.
function [p, q] = i_bycomponent(L, M, g, idx, c, P)
    a = full(sum(M(idx, :), 2));
    k = max(c);  ar = arrayfun(@(i) sum(a(c == i)), 1:k);
    [~, o] = sort(ar, 'descend');
    if ar(o(1)) > 2/3 * sum(ar)
        [p, q] = i_split(L, M, g, idx(c == o(1)), P);
        tot = [sum(a(ismember(idx, p))) sum(a(ismember(idx, q)))];
        for i = o(2:end)
            fr = idx(c == i);
            if tot(1) <= tot(2), p = [p; fr(:)]; tot(1) = tot(1) + ar(i);   %#ok<AGROW>
            else,                q = [q; fr(:)]; tot(2) = tot(2) + ar(i); end %#ok<AGROW>
        end
        p = sort(p);  q = sort(q);
        return
    end
    side = zeros(1, k);  tot = [0 0];
    for i = o                                                   % greedy: the larger pan first
        [~, s] = min(tot);  side(i) = s;  tot(s) = tot(s) + ar(i);
    end
    p = idx(ismember(c, find(side == 1)));
    q = idx(ismember(c, find(side == 2)));
end

function [p, q] = i_bylongest(idx, P)
    ext = max(P(idx, :), [], 1) - min(P(idx, :), [], 1);
    [~, k] = max(ext);
    m = median(P(idx, k));
    p = idx(P(idx, k) >= m);  q = idx(P(idx, k) < m);
    if isempty(q) || isempty(p)
        h = floor(numel(idx)/2);  p = idx(1:h);  q = idx(h+1:end);
    end
end
% Author: Diellor Basha, 2026
