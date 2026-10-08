% +SENSORS  Sensor arrays as graphs: geometry in, a calibrated operator out.
%
% ⚠ SENSOR SPACE ONLY. Every quantity here is a property of the signal ON THE SENSOR GRAPH.
% Speeds are metres per second ACROSS THE ARRAY; wavelengths are millimetres of SENSOR
% SEPARATION. What any of it implies about sources is a question this package does not ask,
% and no function here takes a leadfield, a cortex, or an inverse.
%
% The signal is SCALAR -- one number per sensor per sample -- so the +flow vector machinery
% (curl, divergence, helicity) does not apply. The phase-based readouts do.
%
% Geometry (all return the same array struct: .Name .Kind .Pos .Dim .Pitch .Aperture .Labels .nCh):
%   rheome.sensors.linear    - a laminar probe; the one Dim 1 array, and the only one with no faces
%   rheome.sensors.grid      - a planar lattice; 'utah' and 'ecog' presets
%   rheome.sensors.meg       - MEG helmet positions from a cached study or channel file
%   rheome.sensors.eeg       - a spherical-cap idealisation, or positions from a channel file
%   rheome.sensors.positions - the general case: an array from coordinates already in hand, so a
%                       geometry that survives only in a published figure can still be
%                       analysed on exactly the coordinates it was drawn on
%
% Operator:
%   rheome.sensors.graph     - Gaussian-weighted kNN W and L, returned in the SURFACE shape so the
%                       existing mesh readouts accept it unchanged
%   rheome.sensors.modes     - the spectrum, by exact eigenbasis or by Chebyshev
%
% Generation:
%   rheome.sensors.generate  - a space-time pattern from a dynamical family: seed a delta, let the
%                       family evolve it. An ADAPTER over rheome.filters.impulse, not a method --
%                       the family sets the dynamics, the seed sets where, the graph sets
%                       what the pattern can be. Args are in FAMILY units, not array units.
%
% Measurement:
%   rheome.sensors.sample    - evaluate an analytic field at the sensor coordinates
%   rheome.sensors.calibrate - alpha in lambda ~ alpha k^2, by closed form AND regression, so
%                       lambda carries metres; also the curve whose departure is the
%                       resolution floor
%   rheome.sensors.window    - which wavelengths this array admits, from coordinates alone
%
% ⭐ NO SPEED WITHOUT rheome.sensors.calibrate. The graph Laplacian's lambda is dimensionless;
% rheome.dynamics.dispersion fits omega = c_graph*sqrt(lambda), and the physical speed is
% c_graph*sqrt(alpha). Reporting c_graph as a speed is an error of alpha.
%
% Demos: rheome.demos.sensor_graph, rheome.demos.sensor_limits
%
%   rheome.sensors.tree      - recursive spectral bisection: the position pyramid (hemispheres,
%                       quadrants, patches ... sensors) for the tile store's group rows
%
% Author: Diellor Basha, 2026
