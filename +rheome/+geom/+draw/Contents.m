% DRAW  Surface geometry that returns a curve on the surface.
%
% One of rheome.geom's three jobs (measure, select, draw: the Cortical Flow Editor's and nxr's words).
%
% Functions:
%   rheome.geom.draw.flipOut  - the FlipOut geodesic between two vertices (Sharp & Crane 2020): a port of
%                               geometry-central's, the line the viewer, nxr and the Editor draw (to 1e-9 m)
%   rheome.geom.draw.line     - a geodesic between two vertices on the faces (calls draw.flipOut)
%   rheome.geom.draw.flipMesh - the mesh FlipOut works on, built once for a surface (the Editor's flipMesh): pass it
%                               to line / flipOut as `mesh` for many paths on one surface
%
% Drawing is not measuring: the line does not use the heat distance (rheome.geom.measure.heatDistance). A heat
% solver passed as line's 4th argument (as the old geodesicpath took) is ignored and warns for one release.
%
% The FlipOut port is geometry-central's code (MIT, Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
% contributors): its notice is in +rheome/+geom/+draw/LICENSE-geometry-central and covers flipOut.m, flipMesh.m and
% the helpers private/flipoutMesh.m and private/flipoutTrace.m.
%
% Formerly rheome.geom.flipOut, .geodesicpath (those names forward and warn once).
%
% See also: rheome.geom, rheome.geom.measure, rheome.geom.select
%
% Author: Diellor Basha, 2026
