% RHEOME  The multiscale geometry of human cortical dynamics: graph wavelets, phase and flow.
% Version 1.0.1 (R2023b) 09-Oct-2026
%
% Spatiotemporal dynamics on graphs -- cortical meshes, connectomes and sensor arrays: graph wavelet
% frames, the joint (lambda, omega) plane, phase geometry and flow. Everything lives in ONE namespace,
% rheome.*, so no name collides with MATLAB's or another toolbox's. Each package has its own
% Contents.m: help rheome.flow, help rheome.filters, ...
%
% Classes (the API of record)
%   rheome.graphfilterbank   - Spectral filterbank on a self-adjoint operator.
%   rheome.jointfilterbank   - Filterbank on the joint time-vertex plane (lambda, omega).
%   rheome.timefilterbank    - A designed constant-Q tight frame on the frequency axis.
%   rheome.flowpage          - One band, one page: the LBO coefficients of the cortical vorticity.
%   rheome.flowfeatures      - Windowed, band- and scale-stratified tabulation of the cortical flow tensor.
%   rheome.flowbrowser       - Interactive browser for the band- and scale-resolved cortical flow.
%   rheome.pagedrecording    - A long sensor recording read one time page at a time.
%   rheome.recordingbrowser  - Look at a stored recording through its pyramid, before loading any of it.
%
% Packages
%   rheome.operators         - Discrete differential operators on a triangle surface mesh.
%   rheome.eigen             - Spectral basis of the surface operators.
%   rheome.graphtransform    - Field <-> coefficient transforms for graphfilterbank.
%   rheome.filters           - Spectral filter kernels g(lambda[,omega/t]) and their application.
%   rheome.jtv               - Joint time-vertex filterbanks: construction, duals, analysis and synthesis.
%   rheome.differential      - Divergence, rotation, Helmholtz-Hodge and Poisson on the cortex.
%   rheome.dynamics          - DMD, dispersion, PDE fits and optical flow of activity maps.
%   rheome.spectral          - Periodic / aperiodic decomposition of spectra and of flow fields.
%   rheome.detect            - Feature detection on source vector fields.
%   rheome.flow              - Fused flow imaging kernels (sensors -> cortical flow maps).
%   rheome.forward           - Forward-model spectral transforms.
%   rheome.inverse           - Source-mapping inverse solutions, and what they can resolve.
%   rheome.source            - End-to-end source-mapping orchestration.
%   rheome.sensors           - Sensor arrays as graphs: geometry in, a calibrated operator out.
%   rheome.connectome        - Structural-connectome inputs.
%   rheome.geom              - Mesh generation and geometry utilities.
%   rheome.ingest            - Constant-Q tile summaries of a recording.
%   rheome.select            - The query side of the tile store.
%   rheome.selection         - Tiles, windows, bands, patches, ROIs and wavelet members.
%   rheome.domain            - Where a field lives.
%   rheome.fieldtype         - What a field is.
%   rheome.pipeline          - Operators chained by their arrows, run against an environment.
%   rheome.scale             - The per-participant driver and the group reduction.
%   rheome.qc                - Quality control of raw MEG.
%   rheome.report            - A measured result written up as a mini-article, with its provenance.
%   rheome.io                - Read Brainstorm files from disk.
%   rheome.import            - Brainstorm -> data cache (write once).
%   rheome.load              - Data cache -> analyses, and where outputs go (rheome.load.root, rheome.load.outroot).
%   rheome.show              - Visualizations.
%   rheome.utils             - Shared helpers.
%   rheome.demos             - Runnable, self-checking demonstrations (start here).
%
% Quick start
%   out = rheome.demos.filterbank_graph;   % rheome.graphfilterbank on a sphere, prints PASS/FAIL
%   out.ok
%   (the user guide: Help browser > Supplemental Software > Rheome)
%
% Data and outputs: an installed toolbox caches under <userpath>/rheome/data and writes
% results under <userpath>/rheome/results; set RHEOME_DATA / RHEOME_OUT to move them.
%
% Author: Diellor Basha, 2026
