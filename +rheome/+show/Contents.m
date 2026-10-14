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
%                      (measurement in rheome.inverse.resolution)
%   rheome.show.timetracks  - a pattern's parameters over time (position, scale, heading, chirality, div/curl,
%                      amplitude): the plant's analytic course and the detected values on one time axis
%   rheome.show.track       - an analytic track's course (keys, profile, oscillator): the Editor's tracks
%   rheome.show.phasor      - phase as a phasor glyph: a polar plot of phases or a phasor marker on a map
%   rheome.show.levelstack  - the dyadic levels as the grid over time (no colormap), a field painted on its tiles
%   rheome.show.zstack      - a z-stack indicator: which level of the constant-Q stack a panel shows
%
% Author: Diellor Basha, 2026
