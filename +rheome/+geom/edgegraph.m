function g = edgegraph(S, opts)
% GEOM.EDGEGRAPH  The mesh as a weighted graph whose shortest paths are a geodesic ruler.
%
%   g = rheome.geom.edgegraph(S)                     % edges + intrinsic diagonals (default)
%   g = rheome.geom.edgegraph(S, Diagonals=false)    % mesh edges only
%   d = distances(g, a, b);                   % metres, symmetric, any range
%
% Dijkstra on the mesh edges is intrinsic, symmetric and never catastrophically wrong, but a path
% confined to edges zig-zags: on a head-sized icosphere it overestimates the analytic great-circle
% distance by 3-13%. Adding, across every pair of triangles sharing an edge, the straight segment
% between their two opposite vertices -- measured in the UNFOLDED pair, so intrinsic, and only where
% that segment stays inside the pair -- straightens the path. See the measured table below.
%
% ⭐⭐ WHY THIS EXISTS WHEN rheome.geom.geodesic DOES. The heat method is the better ruler at short range and
% the worse one at long range, and on this cortex long range starts sooner than you would think.
% Measured 2026-09-29 (scratch diag_geo.m, 40 random pairs each):
%   icosphere-5, R = 70 mm, edge 2.64 mm: heat error -0.5..-2.7% out to 130 mm, then it SATURATES --
%     every pair beyond reads 128.4 mm (truth 151-203 mm, error -15..-37%). That is ~49 edges, where
%     the backward-Euler heat kernel (t = h^2, decay ~exp(-d/h)) falls below 1e-16 of its peak and
%     the gradient direction is roundoff.
%   reference left hemisphere, edge 3.85 mm: d(a->b) and d(b->a) DISAGREE -- by 2-11% typically and by
%     67% at 98 mm (87.6 vs 43.7 mm) and 75% at 153 mm (121.1 vs 55.0 mm). An irregular mesh has
%     obtuse cotan weights, the heat solution goes negative, and the direction field flips.
% So a SPEED, which divides a distance of 50-150 mm by a time, must not be read off rheome.geom.geodesic on
% this cortex. rheome.inverse.resolution already uses Dijkstra for r50 and is unaffected.
%
% ⚠ Diagonals can still UNDER-estimate on a sharply folded pair: the unfolded segment is intrinsic,
% but the check that it crosses the shared edge uses the two triangles alone, so across a crease it
% is exact for the pair and nothing more. Paths are built from these, so the bias stays local.
%
% INPUTS
%   S          surface (.Vertices, .Faces)
%   Diagonals  add the intrinsic diagonals of adjacent triangle pairs (true)
%
% OUTPUT
%   g          MATLAB graph, Weight = length in metres (mesh units)
%
% See also: rheome.geom.geodesic, rheome.inverse.resolution, rheome.geom.tiles
%
% Author: Diellor Basha, 2026

    arguments
        S struct
        opts.Diagonals (1,1) logical = true
    end
    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V, 1);
    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
    w = vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2);
    I = E(:,1);  J = E(:,2);

    if opts.Diagonals
        % every face-edge with the vertex opposite it; pairs of faces share a sorted edge key
        he = [F(:,[1 2 3]); F(:,[2 3 1]); F(:,[3 1 2])];          % [a b opposite]
        key = sort(he(:,1:2), 2);
        [~, ~, id] = unique(key, 'rows');
        [id, o] = sort(id);  he = he(o,:);  key = key(o,:);
        pr = find(id(1:end-1) == id(2:end));                      % interior edges: two faces
        a = key(pr,1);  b = key(pr,2);  c = he(pr,3);  d = he(pr+1,3);
        lab = vecnorm(V(a,:) - V(b,:), 2, 2);
        [xc, yc] = i_place(vecnorm(V(a,:) - V(c,:), 2, 2), vecnorm(V(b,:) - V(c,:), 2, 2), lab);
        [xd, yd] = i_place(vecnorm(V(a,:) - V(d,:), 2, 2), vecnorm(V(b,:) - V(d,:), 2, 2), lab);
        yd = -yd;                                                 % unfold: opposite side of ab
        s = yc ./ (yc - yd);                                      % where cd crosses the ab line
        xs = xc + s .* (xd - xc);
        ok = xs > 0 & xs < lab & c ~= d;                          % convex pair: segment stays inside
        I = [I; c(ok)];  J = [J; d(ok)];
        w = [w; hypot(xd(ok) - xc(ok), yd(ok) - yc(ok))];
    end
    g = simplify(graph(I, J, w, nV), 'min');
end

% Place a vertex at distances (ra, rb) from a = (0,0) and b = (L,0), above the axis.
function [x, y] = i_place(ra, rb, L)
    x = (ra.^2 - rb.^2 + L.^2) ./ (2 * L);
    y = sqrt(max(ra.^2 - x.^2, 0));
end

% Author: Diellor Basha, 2026
