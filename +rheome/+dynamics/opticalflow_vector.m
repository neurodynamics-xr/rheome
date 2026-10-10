function velocityField = opticalflow_vector(currentField, surface, opts)
% DYNAMICS.OPTICALFLOW_VECTOR  Vector optical flow of an unconstrained current field.  ⚠ RETIRED
%
% ⚠ RETIRED 2026-10-08: on a NORMAL current J = I*n -- what an
% inverse returns -- it fails noise-free on the FreeSurfer sphere (direction 25-69 deg, speed 35-89%, div/curl
% 27-64% of truth at alpha <= 1), because n turns along the flow: the static geometry of n is texture
% that says "no motion". It is exact only when its model holds (three independently advected channels:
% 1-5 deg, 96-104%), and there it is worse than the scalar under noise, since its data term carries no
% mass weight and alpha below ~1e3 does not regularise. The 'covariant' transport was never implemented.
% Use rheome.dynamics.opticalflow_scalar on the envelope. Kept callable; warns on every call.
%
%   velocityField = rheome.dynamics.opticalflow_vector(currentField, surface [,opts])
%
% Estimates the propagation velocity v(x,t) by advecting the WHOLE current vector field J (three
% ambient components) rather than a single scalar. Because the three components give three
% brightness-constancy constraints for the two tangent-velocity DOF, the per-vertex problem is
% OVER-DETERMINED -- aperture-free where at least two components have independent spatial gradients,
% unlike the scalar Horn-Schunck baseline. Per vertex, with p_c = e1.gradJ_c, q_c = e2.gradJ_c, solve
% the structure-tensor least squares  [sum p_c^2, sum p_c q_c; sum p_c q_c, sum q_c^2][a;b] =
% -[sum p_c dJ_c; sum q_c dJ_c], plus a light global smoothness prior alpha*K.
%
% INPUTS:
%   currentField [3nV x nT] ambient per-vertex 3-vector current, rows [x1;y1;z1;x2;...]
%   surface      struct .Vertices .Faces .VertNormals .nV
%   opts.alpha       smoothness weight (default 0.1)
%   opts.transport   'ambient' (default; multi-channel) | 'covariant' (connection-transported)
% OUTPUT:
%   velocityField [3nV x (nT-1)] ambient tangent velocity per frame
%
% See also: rheome.dynamics.opticalflow_scalar, rheome.dynamics.flow_readout, rheome.operators.connection_laplacian
%
% Author: Diellor Basha, 2026

    warning('dynamics:opticalflow_vector:retired', ...
        'rheome.dynamics.opticalflow_vector is retired (fails on a normal current; see its help). Use rheome.dynamics.opticalflow_scalar.');
    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'alpha')     || isempty(opts.alpha),     opts.alpha = 0.1;        end
    if ~isfield(opts, 'transport') || isempty(opts.transport), opts.transport = 'ambient'; end

    % scale-normalise the current amplitude (same rationale as the scalar OF)
    currentField = currentField / max(abs(currentField(:)) + eps);

    fg = rheome.operators.face_gradient(surface.Vertices, surface.Faces);
    [stiffness, ~] = rheome.operators.laplace_beltrami(surface.Vertices, surface.Faces, 'galerkin');
    tb = rheome.dynamics.tangent_basis(surface);
    nV = surface.nV;  nT = size(currentField, 2);
    velocityField = zeros(3*nV, nT-1);
    isCovariant = strcmpi(opts.transport, 'covariant');
    if isCovariant, conn = rheome.operators.connection_laplacian(surface.Vertices, surface.Faces); end

    for t = 1:nT-1
        currentNow = currentField(:, t);
        currentDot = currentField(:, t+1) - currentField(:, t);
        Spp = zeros(nV,1); Sqq = Spp; Spq = Spp; bp = Spp; bq = Spp;
        for c = 1:3
            comp    = currentNow(c:3:end);                       % component field [nV x 1]
            compDot = currentDot(c:3:end);
            grad = [fg.W*(fg.Gx*comp), fg.W*(fg.Gy*comp), fg.W*(fg.Gz*comp)];   % [nV x 3]
            p = sum(grad .* tb.e1, 2);
            q = sum(grad .* tb.e2, 2);
            Spp = Spp + p.^2;   Sqq = Sqq + q.^2;   Spq = Spq + p.*q;
            bp  = bp  + p.*compDot;  bq = bq + q.*compDot;
        end
        if isCovariant
            [Spp, Sqq, Spq, bp, bq] = i_covariant_correction(Spp, Sqq, Spq, bp, bq, ...
                currentNow, currentDot, tb, conn, fg);
        end
        A = [spdiags(Spp,0,nV,nV) + opts.alpha*stiffness,  spdiags(Spq,0,nV,nV); ...
             spdiags(Spq,0,nV,nV),  spdiags(Sqq,0,nV,nV) + opts.alpha*stiffness];
        ab = A \ (-[bp; bq]);
        velocity = ab(1:nV).*tb.e1 + ab(nV+1:end).*tb.e2;
        velocityField(:, t) = reshape(velocity', [], 1);
    end
end

% ---- Task-3 deliverable: connection-transported directional derivative. No-op stub for 'ambient'. ----
function [Spp, Sqq, Spq, bp, bq] = i_covariant_correction(Spp, Sqq, Spq, bp, bq, currentNow, currentDot, tb, conn, fg) %#ok<INUSD>
    % Implemented in Task 3 (uses conn.Rt to transport the complex tangent field before differencing).
end

% Author: Diellor Basha, 2026
