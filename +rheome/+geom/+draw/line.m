function [P, info] = line(S, from, to, mesh)
% GEOM.DRAW.LINE  A geodesic on the surface between two vertices, traced across the faces (FlipOut).
%
%   [P, info] = rheome.geom.draw.line(S, from, to)
%   [P, info] = rheome.geom.draw.line(S, from, to, mesh)   % mesh = rheome.geom.draw.flipMesh(S), or info.mesh
%
% A path ON the folded surface: rheome.geom.draw.flipOut, geometry-central's FlipOut (Sharp & Crane 2020) ported to
% MATLAB -- the shortest path along the edges, straightened by edge flips until it is locally shortest, traced back
% onto the faces. It is the line every lab tool draws for "the geodesic between two vertices" (the viewer through
% nxr-compute, nxr, the Editor), matched point for point (tests/tFlipOut.m). Its points lie on the mesh, so its
% length is a length along the surface. It is locally shortest, not always the shortest: the exact distance (MMP)
% can be a few percent less on a folded cortex (rheome.snapshot reports both). Drawing is not measuring: the line
% does not use the heat distance (rheome.geom.measure.heatDistance).
%
% INPUTS:
%   S       surface (.Vertices, .Faces)        from, to   vertex indices
%   mesh    (optional) rheome.geom.draw.flipMesh(S), or info.mesh from an earlier path on S (built once)
%
% OUTPUTS:
%   P     [n x 3] the path's points, from `from` to `to`: vertices, edge crossings and (rarely) face points
%   info  .length (mesh units, along the surface)  .faces [n-1 x 1] the face each segment lies in
%         .vertices, .weights [n x 3] where each point is: P(i,:) = weights(i,:) * S.Vertices(vertices(i,:),:)
%         .mesh (the flip mesh: pass it back for the next path on S)
%         .solver (deprecated, one release: the 4th argument as passed, unchanged; [] if none)
%
% Deprecated for one release: a heat solver (anything that is not a flip mesh) as the 4th argument. It was the
% heat-descent tracer's; the line no longer uses it. It is ignored, warns once per session (rheome:deprecated), and
% comes back as info.solver, so callers that chain heat solvers keep working.
%
% See also: rheome.geom.draw.flipMesh, rheome.geom.draw.flipOut, rheome.geom.measure.heatDistance
%
% Author: Diellor Basha, 2026

    persistent warned
    if nargin < 4, mesh = []; end
    solver = mesh;
    if ~isempty(mesh) && ~isFlipMesh(mesh)
        if isempty(warned)
            warned = true;
            warning('rheome:deprecated', ['rheome.geom.draw.line: the 4th argument is a flip mesh; a heat solver is ' ...
                'ignored (the line no longer uses the heat distance) and will be refused in the next release. ' ...
                'Pass rheome.geom.draw.flipMesh(S), or nothing.']);
        end
        mesh = [];
    end
    [P, f] = rheome.geom.draw.flipOut(S, from, to, mesh);
    info = struct('length', f.length, 'faces', f.faces, 'vertices', f.vertices, 'weights', f.weights, ...
        'mesh', f.mesh, 'solver', {solver});
end
