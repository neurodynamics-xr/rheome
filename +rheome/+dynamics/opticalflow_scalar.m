function velocityField = opticalflow_scalar(activityMap, surface, opts)
% DYNAMICS.OPTICALFLOW_SCALAR  Manifold Horn-Schunck optical flow of a scalar activity map.
%
%   velocityField = rheome.dynamics.opticalflow_scalar(activityMap, surface [,opts])
%   velocityField = rheome.dynamics.opticalflow_scalar(A, S, struct('regulariser','frame','alpha',1))   % the pre-2026-11-11 method
%
% Estimates the propagation velocity field v(x,t) of a SCALAR activity map (e.g. |J|, the source
% magnitude -- what Brainstorm's bst_opticalflow uses) by the manifold Horn-Schunck variational
% method: per frame it minimises the brightness-constancy data term plus an isotropic smoothness
% prior,  min_v  sum_v M_v (dI/dt + v.gradI)^2  +  alpha * R(v).  The regulariser R is a choice
% (PI decision, 2026-11-11):
%
%   'connection' (DEFAULT)  R = v'Cv, C the CONNECTION LAPLACIAN of the intrinsic Delaunay mesh -- the Cortical Flow Editor's
%       method and nxr's (rheome.dynamics.opticalflow_connection: the Editor's data term, signpost gauge, alpha default 0.01).
%       Frame-free: v'Cv compares a vector with its neighbour's parallel transport.
%   'frame'                 R = a'Ka + b'Kb, the scalar Laplacian on the COEFFICIENTS of v = a*e1 + b*e2 in the (1,1,1)xn
%       frames of rheome.dynamics.tangent_basis (alpha default 1). Kept only to reproduce old results.
%
% ⚠ THE 'frame' REGULARISER DEPENDS ON THE FRAME. Two neighbouring frames turned against each other make a parallel
%   field look rough, and the (1,1,1)xn frames turn across every fold. On a real 1 s MEG alpha window its divergence/curl
%   maps correlate 0.05/0.00 with the connection-Laplacian flow and its speed is nearly uniform (an artefact of the
%   prior, not a property of the data) -- nxr-neuro audit 2026-11-11-opticalflow-smoothness. ⚠ Results made with the old
%   default (every flow.apparent / measure_* result before 2026-11-11) should be re-checked with 'connection'.
% ⚠ alpha IS NOT COMPARABLE ACROSS REGULARISERS: C and a'Ka are different penalties (0.01 for 'connection' is not 0.01 for
%   'frame'). Each regulariser has its own default; record both.
%
% The 'frame' form reuses rheome.operators.face_gradient (gradient), rheome.operators.laplace_beltrami (smoothness K,
% mass M), rheome.dynamics.tangent_basis. This is the scalar BASELINE for comparison against the vector
% optical flow (rheome.dynamics.opticalflow_vector), which is aperture-free.
%
% INPUTS:
%   activityMap  [nV x nT] scalar activity per vertex per frame
%   surface      struct with .Vertices .Faces .nV (and .VertNormals for 'frame')
%   opts.regulariser  'connection' (default) | 'frame'
%   opts.alpha   smoothness weight (default 0.01 for 'connection', 1 for 'frame')
%   opts.mollify, opts.flipEps  ('connection') as rheome.operators.intrinsic_delaunay
% OUTPUT:
%   velocityField [3nV x (nT-1)] ambient tangent velocity, rows [x1;y1;z1;x2;...] per frame
%
% See also: rheome.dynamics.opticalflow_connection, rheome.dynamics.opticalflow_vector, rheome.dynamics.flow_readout, rheome.dynamics.tangent_basis
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'regulariser') || isempty(opts.regulariser), opts.regulariser = 'connection'; end
    reg = char(lower(string(opts.regulariser)));
    if strcmp(reg, 'connection')
        velocityField = rheome.dynamics.opticalflow_connection(activityMap, surface, opts);
        return
    elseif ~strcmp(reg, 'frame')
        error('dynamics:opticalflow_scalar:regulariser', 'regulariser must be ''connection'' or ''frame'', got ''%s''.', reg);
    end
    if ~isfield(opts, 'alpha') || isempty(opts.alpha), opts.alpha = 1; end

    % Normalise the activity map to unit amplitude (as Brainstorm's bst_opticalflow does): the true
    % velocity is invariant to scaling the brightness, but the DATA-vs-SMOOTHNESS balance is not, so
    % this makes the smoothness weight alpha interpretable and comparable across datasets.
    activityMap = activityMap / max(abs(activityMap(:)) + eps);

    fg = rheome.operators.face_gradient(surface.Vertices, surface.Faces);
    [stiffness, massMatrix] = rheome.operators.laplace_beltrami(surface.Vertices, surface.Faces, mesh='embedded', mass='galerkin');
    tb = rheome.dynamics.tangent_basis(surface);
    massVec = full(sum(massMatrix, 2));
    nV = surface.nV;  nT = size(activityMap, 2);
    velocityField = zeros(3*nV, nT-1);
    % ⭐ The smoothness block and the sparsity pattern are the same every frame: build them once and add
    % the four data diagonals with one sparse(). Same entries as the spdiags/concatenation form (each is
    % one addition, which IEEE makes order-free), so the solve sees the same matrix and returns the same
    % bits -- checked on a reference subject: assembly 1.8 -> 0.7 ms per frame against a 36 ms solve.
    aK = opts.alpha*stiffness;  Z = sparse(nV, nV);  smooth = [aK Z; Z aK];
    iD = [1:nV, 1:nV, nV+1:2*nV, nV+1:2*nV]';  jD = [1:nV, nV+1:2*nV, 1:nV, nV+1:2*nV]';

    for t = 1:nT-1
        activity = activityMap(:, t);
        activityDot = activityMap(:, t+1) - activityMap(:, t);
        vertexGradient = [fg.W*(fg.Gx*activity), fg.W*(fg.Gy*activity), fg.W*(fg.Gz*activity)];  % [nV x 3]
        g1 = sum(vertexGradient .* tb.e1, 2);
        g2 = sum(vertexGradient .* tb.e2, 2);
        g12 = massVec.*g1.*g2;
        A = smooth + sparse(iD, jD, [massVec.*g1.^2; g12; g12; massVec.*g2.^2], 2*nV, 2*nV);
        rhs = -[massVec.*g1.*activityDot; massVec.*g2.*activityDot];
        ab = A \ rhs;
        velocity = ab(1:nV).*tb.e1 + ab(nV+1:end).*tb.e2;                % [nV x 3] tangent
        velocityField(:, t) = reshape(velocity', [], 1);
    end
end

% Author: Diellor Basha, 2026
