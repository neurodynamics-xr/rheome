% +DEMOS  Runnable demonstrations.
%
% FILTERBANK SUITE -- three demos in the SAME five-figure order, so the analogy between
% MATLAB's cwtfilterbank and our two classes reads by scrolling. Each plants features with
% an analytic answer and prints PASS/FAIL, so running one is also a validation.
%
%   rheome.demos.filterbank_cwt    - MATLAB's cwtfilterbank: a chirp and a burst
%                             (needs the Wavelet Toolbox; SKIPS cleanly without it)
%   rheome.demos.filterbank_graph  - graphfilterbank: Gaussian bumps of known width on a sphere
%   rheome.demos.filterbank_joint  - jointfilterbank: a rigidly rotating pattern at known speed
%
% Each takes an optional outDir and exports its five figures as PNGs when given one:
%   rheome.demos.filterbank_graph('_figures')
%
% ⚠ THE THREE CALIBRATIONS. rheome.demos.filterbank_graph recovers sqrt(2)*sigma, not sigma,
% because vertexSpectrum is an ENERGY marginal and the relevant integral is of the SQUARED
% gain. A linear matched response would give sigma; rheome.filters.frame's Gamma = 1/sqrt(2)
% applies to a VORTEX. Confusing them is a 2x error that still looks like a measurement.
%
% Surface and operator demos:
%   rheome.demos.wavelet           - localize spectral filters as atoms on a cortex
%   rheome.demos.laplace_beltrami  - the LBO and its eigenmodes
%   (both take a Brainstorm cortex surface .mat)
%
% Sensor-array demos (sensor space only -- no leadfield, no cortex, no inverse):
%   rheome.demos.sensor_graph      - coordinates to a calibrated operator, in five figures
%   rheome.demos.sensor_limits     - each array's measurement window, predicted then swept until it breaks
%   rheome.demos.sensor_dynamics - what a moving pattern looks like on an array, and in its spectrum
%   rheome.demos.sensor_phase    - instantaneous phase across an array, and what its gradient measures
%   rheome.demos.sensor_topology - divergence and curl of the phase gradient, and what they classify
%   rheome.demos.sensor_tensor   - the joint wavelet tensor: aperture, rate and dispersion as controls
%   rheome.demos.sensor_wavelets - scale, rate and speed PLANTED then RECOVERED, generation and
%                           analysis through the same three factor lists
%
% Author: Diellor Basha, 2026
