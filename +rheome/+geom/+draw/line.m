function [P, info] = line(S, from, to, solver)
% GEOM.DRAW.LINE  A geodesic on the surface between two vertices, traced across the faces (FlipOut).
%
%   [P, info] = rheome.geom.draw.line(S, from, to)
%   [P, info] = rheome.geom.draw.line(S, from, to, solver)   % reuse info.solver from a path on the same S
%
% A path ON the folded surface: rheome.geom.draw.flipOut, geometry-central's FlipOut (Sharp & Crane 2020) ported to
% MATLAB -- the shortest path along the edges, straightened by edge flips until it is locally shortest, traced back
% onto the faces. It is the line every lab tool draws for "the geodesic between two vertices" (the viewer through
% nxr-compute, nxr, the Editor), matched point for point (tests/tFlipOut.m). Its points lie on the mesh, so its
% length is a length along the surface. It is locally shortest, not always the shortest: the exact distance (MMP)
% can be a few percent less on a folded cortex (rheome.snapshot reports both).
%
% INPUTS:
%   S       surface (.Vertices, .Faces)        from, to   vertex indices
%   solver  (optional) info.solver from an earlier path on the same S (the halfedge mesh, built once)
%
% OUTPUTS:
%   P     [n x 3] the path's points, from `from` to `to`: vertices, edge crossings and (rarely) face points
%   info  .length (mesh units, along the surface)  .faces [n-1 x 1] the face each segment lies in
%         .vertices, .weights [n x 3] where each point is: P(i,:) = weights(i,:) * S.Vertices(vertices(i,:),:)
%         .solver (pass it back for the next path on S)
%
% See also: rheome.geom.draw.flipOut, rheome.geom.measure.heatDistance, rheome.geom.edgegraph
%
% Author: Diellor Basha, 2026

    if nargin < 4 || ~isstruct(solver) || ~isfield(solver, 'signpost'), solver = []; end
    [P, f] = rheome.geom.draw.flipOut(S, from, to, solver);
    info = struct('length', f.length, 'faces', f.faces, 'vertices', f.vertices, 'weights', f.weights, 'solver', f.cache);
end

% Author: Diellor Basha, 2026
