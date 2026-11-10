% DRAW  Surface geometry that returns a curve on the surface.
%
% One of rheome.geom's three jobs (measure, select, draw: the Cortical Flow Editor's and nxr's words).
%
% Functions:
%   rheome.geom.draw.flipOut  - the FlipOut geodesic between two vertices (Sharp & Crane 2020): a port of
%                               geometry-central's, the line the viewer, nxr and the Editor draw (to 1e-9 m)
%   rheome.geom.draw.line     - a geodesic between two vertices on the faces (calls draw.flipOut)
%
% The FlipOut port is geometry-central's code (MIT, Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
% contributors): its notice is in +rheome/+geom/+draw/LICENSE-geometry-central and covers flipOut.m and its helpers
% private/flipoutMesh.m and private/flipoutTrace.m.
%
% Formerly rheome.geom.flipOut, .geodesicpath (those names forward and warn once).
%
% See also: rheome.geom, rheome.geom.measure, rheome.geom.select
%
% Author: Diellor Basha, 2026
