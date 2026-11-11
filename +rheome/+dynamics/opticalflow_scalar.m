function velocityField = opticalflow_scalar(activityMap, surface, opts)
% DYNAMICS.OPTICALFLOW_SCALAR  Manifold Horn-Schunck optical flow of a scalar activity map.
%
%   velocityField = rheome.dynamics.opticalflow_scalar(activityMap, surface [,opts])
%
% Estimates the propagation velocity field v(x,t) of a SCALAR activity map (e.g. |J|, the source magnitude) by the manifold
% Horn-Schunck variational method: per frame it minimises the brightness-constancy data term plus an isotropic smoothness prior,
%   min_v  sum_v M_v (dI/dt + v.gradI)^2  +  alpha * v'Cv,
% C the CONNECTION LAPLACIAN of the intrinsic Delaunay mesh: the Cortical Flow Editor's method and nxr's, frame-free (v'Cv
% compares a vector with its neighbour's parallel transport). This is a call of rheome.dynamics.opticalflow_connection; see it
% for the method, the data term, the gauge and the Editor parity.
%
% ⚠ THE EARLIER 'frame' SMOOTHNESS WAS REMOVED (PI ruling, 2026-11-11: "we should not be keeping any code that implements math
%   and analysis wrongly"). It penalised the scalar Laplacian of the COEFFICIENTS of v in per-vertex (1,1,1)xn frames, which
%   depends on the frame: neighbouring frames turned against each other make a parallel field look rough. On a real 1 s MEG
%   alpha window its divergence/curl maps correlated 0.05/0.00 with this flow (nxr-neuro audit 2026-11-11-opticalflow-smoothness).
%   Any result made with it (flow.apparent / measure_* before Rheome's 2026-11-11 releases) must be recomputed.
% ⚠ alpha: default 0.01 (the Editor's). v'Cv is not zero on a rigid rotation of a sphere, so the prior pulls the speed low as
%   alpha grows (sphere bump: 95 % of truth at 0.01, 33 % at 1). Report the alpha used.
%
% INPUTS:
%   activityMap  [nV x nT] non-negative scalar activity per vertex per frame
%   surface      struct with .Vertices .Faces .nV
%   opts.alpha   smoothness weight (default 0.01);  opts.mollify, opts.flipEps  as rheome.operators.intrinsic_delaunay
% OUTPUT:
%   velocityField [3nV x (nT-1)] ambient tangent velocity, rows [x1;y1;z1;x2;...] per frame
%
% See also: rheome.dynamics.opticalflow_connection, rheome.dynamics.flow_readout, rheome.flow.apparent
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    velocityField = rheome.dynamics.opticalflow_connection(activityMap, surface, opts);
end
