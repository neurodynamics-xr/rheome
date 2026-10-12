% +MODES  Resting-state pattern atlas: data-driven discovery of large-scale (>= 100 mm) complex spatial modes.
%
% Which large-scale spatiotemporal patterns does resting cortex express, within what MEG resolves? The modes
% are complex spatial patterns (amplitude topography and phase layout) on the shared group tiles, discovered
% by phase-invariant clustering, and tested against nulls of the instrument (leakage, lagged sources), of
% anatomy (spin), of the background (aperiodic) and of second-order structure (phase randomisation).
%
% Per participant
%   operator  - sensors -> group-tile normal current, low-passed at a declared wavelength (one matrix)
%   fields    - complex band tile series z(t) + nulls N1 (sources, empty room), N4 (aperiodic), N5 (phase-randomised)
%   analytic  - one tile series -> periodic/aperiodic split, band-pass, analytic signal, decimated
%   jointfields - joint graph-time wavelet field W(v,s,t,f) on tile x scale per octave + null F (field spread),
%               N4, N5, and the scale transfer T(s); Canary=true plants the J5 patterns on F instead
%   scaletransfer - T(s): the share of a planted graph-scale pattern the instrument returns at that scale, the edge
%   fit       - back-fit frozen centroids: occupancy, dwell, transitions, HMM, split-half, physical features, N2, N3, N5
% Group
%   discover  - pooled phase-invariant k-means, K by split-half reproducibility, DMD and cPCA sensitivity, frozen + hashed
%   jointcheck - J5: joint discovery on the canary plants, recovery |r| per plant, F against Level-1 (i)
%   reduce    - the keep rule (Level 1, labels, Level 2) per mode x band, Holm within cohort; the catalogue table
% Building blocks
%   octave    - the joint field's time factor: a jointfilterbank band member as an analytic, decimated signal
%   cluster   - phase-invariant complex k-means (centroid = principal eigenvector of the cluster scatter)
%   match     - Hungarian matching of two mode sets by |r|
%   hash      - SHA-256 of frozen centroids
%
% Design: the resting-state pattern atlas study protocol, sections 4-7; joint atoms (fields, discover, fit and
% reduce with M.kind "joint"): its section 14.
%
% Author: Diellor Basha, 2026
