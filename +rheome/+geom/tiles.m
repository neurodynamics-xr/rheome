function G = tiles(T, S, depth, opts)
% GEOM.TILES  The tiles of one rheome.geom.tree depth: vertex lookup, membership and adjacency.
%
%   G = rheome.geom.tiles(T, S, depth)
%   T = rheome.geom.tree(S, L=L, M=M, MaxDepth=3);  G = rheome.geom.tiles(T, S, 3);
%   G.A(G.tileOf(v1), G.tileOf(v2))          % > 0 iff the two vertices sit in touching tiles
%   G = rheome.geom.tiles(T, S, 3, Ruler=rheome.geom.edgegraph(S))   % also tile centres and the distances between
%
% A tile is a node of the position pyramid. This turns one level of it into a graph: the tiles are
% its nodes, and two tiles are adjacent when a mesh edge crosses between them. The adjacency is one
% sparse product, Atile = P' * A * P, with P the vertex-to-tile membership and A the mesh adjacency.
%
% ⭐ WHAT IT IS FOR IN TRACKING: AN EXACT CANDIDATE PRUNING. A tracker gated at `maxStep` smaller than
% a tile can only move a feature into its own tile or an adjacent one, so the adjacency lists every
% legal destination without computing a distance. The tile sequence of a track is then a coarse
% trajectory that can be queried -- derived from the fine track, never a substitute for it: at depth 3
% a 133 mm tile takes a 0.1 m/s peak over a second to cross, and a vertex-level peak within ~50 mm of a
% boundary (the localisation error at the pipeline's SNR) will flicker across it.
%
% ⭐ `.A` COUNTS CUT EDGES, not 0/1. The count is a discrete length of shared boundary, so a pair of
% tiles touching along a whole sulcus is distinguishable from a pair meeting at a corner. `.A > 0` is
% the plain adjacency.
%
% ⚠ THE NODE ARITHMETIC IS CONDITIONAL. rheome.geom.tree numbers nodes breadth-first, so the children of node
% k are 2k and 2k+1 -- and a tile's ancestor at any depth is floor(k / 2^steps) -- ONLY if every node
% down to this depth split. One early MinArea or MinVertices stop shifts every later id. `.heap` says
% whether the identity holds on this tree; check it before indexing by arithmetic.
%
% ⚠ Tiles are taken at EXACTLY `depth`. A leaf that stopped shallower is carried as its own tile at
% its own depth (`.depth` says which), so the tiles still partition every vertex.
%
% ⚠ The whole cortex is two components: tiles in different hemispheres are never adjacent. A feature
% that changes hemisphere between frames has jumped, not moved.
%
% ⭐ A TRAJECTORY ON TILES IS MEASURED BETWEEN TILE CENTRES. `.centre` is the member vertex nearest
% the tile's centroid (a centroid of a folded patch can sit off the sheet), and `.D` the Ruler's
% distance between every pair of centres -- so a track's displacement and path length are read on
% the same ruler as the planted truth (rheome.geom.edgegraph, +1.4% against the analytic geodesic).
%
% INPUTS
%   T      rheome.geom.tree table            S  the surface T was built on (.Faces, .Vertices)
%   depth  tree depth to tile at
%   Ruler  (optional) a graph on the vertices, e.g. rheome.geom.edgegraph(S); enables .centre and .D
%
% OUTPUT (struct G)
%   .node_id [nTile x 1]   .depth [nTile x 1]   .diameterMM [nTile x 1]   .centroid [nTile x 3]
%   .tileOf  [nV x 1] tile index (1..nTile) of every vertex
%   .P       [nV x nTile] sparse logical membership
%   .A       [nTile x nTile] sparse, cut-edge count between tiles, zero diagonal
%   .heap    true iff parent_id == floor(node_id/2) for every node of T
%   .centre  [nTile x 1] centre vertex; .D [nTile x nTile] Ruler distance between centres (m)
%            (both [] without a Ruler)
%
% See also: rheome.geom.tree, rheome.detect.peaks, rheome.detect.track
%
% Author: Diellor Basha, 2026

    arguments
        T table
        S struct
        depth (1,1) double
        opts.Ruler = []
    end
    sel = T.depth == depth | (T.is_leaf & T.depth < depth);
    Tt = T(sel, :);
    nV = size(S.Vertices, 1);  nTile = height(Tt);

    tileOf = zeros(nV, 1);
    for k = 1:nTile, tileOf(Tt.members{k}) = k; end
    if any(tileOf == 0)
        error('geom:tiles:partition', '%d vertices belong to no tile at depth %d.', ...
            sum(tileOf == 0), depth);
    end

    F = double(S.Faces);
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    ti = tileOf(E(:,1));  tj = tileOf(E(:,2));  cut = ti ~= tj;
    A = sparse([ti(cut); tj(cut)], [tj(cut); ti(cut)], 1, nTile, nTile);

    G.node_id = Tt.node_id;
    G.depth = Tt.depth;
    G.diameterMM = 1e3 * Tt.diameter;
    G.centroid = Tt.centroid;
    G.tileOf = tileOf;
    G.P = sparse((1:nV)', tileOf, true, nV, nTile);
    G.A = A;
    nr = T.parent_id ~= 0;
    G.heap = all(T.parent_id(nr) == floor(T.node_id(nr) / 2));
    G.centre = [];  G.D = [];
    if ~isempty(opts.Ruler)
        G.centre = zeros(nTile, 1);
        for k = 1:nTile
            mk = Tt.members{k}(:);
            [~, j] = min(vecnorm(double(S.Vertices(mk,:)) - Tt.centroid(k,:), 2, 2));
            G.centre(k) = mk(j);
        end
        G.D = distances(opts.Ruler, G.centre, G.centre);
    end
end

% Author: Diellor Basha, 2026
