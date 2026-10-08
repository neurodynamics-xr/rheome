% SHOW  Visualizations and plots for the atom-design demo.
%
% Rendering helpers built on base MATLAB graphics.
%
% Functions:
%   rheome.show.surface     - render a cortical surface, optionally colored by a per-vertex field
%   rheome.show.vectors     - render a per-vertex 3-vector field (magnitude + quiver arrows)
%   rheome.show.spectrum2d  - joint eigenmode-frequency map (power / gain over sqrt(lambda), Hz)
%   rheome.show.filmstrip   - cortex snapshots of a scalar space-time field at time frames
%   rheome.show.vectorstrip - cortex quiver snapshots of a space-time VECTOR field at time frames
%   rheome.show.gif         - animated GIF of a space-time field evolving on the cortex
%   rheome.show.resolution  - the spatial scales of the problem on one log ruler, with rows of cortex
%                      panels above it: the pyramid's rungs, a measured point-spread function
%                      and real atlas scouts, all at ONE camera so the sizes compare
%                      (resolution_scales_omega.m; measurement in rheome.inverse.resolution)
%
% Author: Diellor Basha, 2026
