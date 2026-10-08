% DYNAMICS  The kinematic route: DMD, dispersion, PDE fits and optical flow of activity maps.
%
% Functions:
%   rheome.dynamics.dispersion         - Empirical dispersion relation omega(lambda) from a DMD result.
%   rheome.dynamics.dmd                - Reduced-rank Dynamic Mode Decomposition of a coefficient time series.
%   rheome.dynamics.flow_readout       - GCF readouts of a propagation velocity field.
%   rheome.dynamics.mode_coefficients  - Reduced-order coefficient time series for the dynamics fits.
%   rheome.dynamics.opticalflow_scalar - Manifold Horn-Schunck optical flow of a scalar activity map.
%   rheome.dynamics.opticalflow_vector - Vector optical flow of an unconstrained current field.
%   rheome.dynamics.pde_fit            - Fit the per-mode temporal dynamics to candidate PDEs (the model-based law).
%   rheome.dynamics.tangent_basis      - Per-vertex orthonormal tangent frame from the vertex normals.
%
% Author: Diellor Basha, 2026
