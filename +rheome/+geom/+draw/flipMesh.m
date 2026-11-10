function M = flipMesh(S)
% GEOM.DRAW.FLIPMESH  The mesh FlipOut works on, built once for a surface (the Editor's flipMesh).
%
%   M = rheome.geom.draw.flipMesh(S)
%   [P, info] = rheome.geom.draw.line(S, from, to, M)      % many paths on S share it
%   [P, info] = rheome.geom.draw.flipOut(S, from, to, M)
%
% The halfedge mesh, the input geometry and the start of the intrinsic triangulation that
% rheome.geom.draw.flipOut reads (geometry-central's ManifoldSurfaceMesh, VertexPositionGeometry and
% SignpostIntrinsicTriangulation, ported in private/flipoutMesh.m). Each path flips its own copy, so one M serves
% every path on S. It is the mesh only: drawing a line does not use the heat distance
% (that is rheome.geom.measure.heatDistance).
%
% INPUTS:
%   S   surface (.Vertices [nV x 3], .Faces [nF x 3], 1-based, oriented, manifold; a boundary is allowed)
%
% OUTPUTS:
%   M   the flip mesh (a struct; also info.mesh of any line / flipOut call on S)
%
% See also: rheome.geom.draw.line, rheome.geom.draw.flipOut
%
% Author: Diellor Basha, 2026 (a port of geometry-central, MIT; licence: +rheome/+geom/+draw/LICENSE-geometry-central)

    M = flipoutMesh(S.Vertices, S.Faces);
end
