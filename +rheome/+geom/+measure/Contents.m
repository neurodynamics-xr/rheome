% MEASURE  Surface geometry that returns numbers: distances on the surface, the intrinsic mean, a field pooled over tiles.
%
% One of rheome.geom's three jobs (measure, select, draw: the Cortical Flow Editor's and nxr's words).
%
% Functions:
%   rheome.geom.measure.heatDistance      - single-/multi-source geodesic distance by the heat method (geometry-
%                        central's robust solver: intrinsic Delaunay, mollified, lumped mass) + reusable solver
%   rheome.geom.measure.heatDistanceEigen - the heat method in the LBO eigenbasis: two diagonal steps, no solve;
%                        accuracy set by the number of modes
%   rheome.geom.measure.karcherMean       - intrinsic (geodesic Fréchet) mean of points on the surface: the vertex
%                        minimising the weighted sum of heat distances^power (stays on the sheet)
%   rheome.geom.measure.tileMean          - a vertex map pooled over each tile (area-weighted mean, or median)
%
% Formerly rheome.geom.geodesic, .heatdistance, .karcher_mean, .tilemean (those names forward and warn once).
%
% See also: rheome.geom, rheome.geom.select, rheome.geom.draw, rheome.geom.edgegraph
%
% Author: Diellor Basha, 2026
