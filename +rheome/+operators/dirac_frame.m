function [A, B, scales] = dirac_frame(V, F, tau, N, normalize)
% OPERATORS.DIRAC_FRAME  Relative Dirac operator (Brainstorm production).
%
%   [A, B]         = rheome.operators.dirac_frame(V, F, tau)
%   [A, B]         = rheome.operators.dirac_frame(V, F, tau, N)              % supply the Gauss map
%   [A, B]         = rheome.operators.dirac_frame(V, F, tau, N, normalize)   % scale convention
%   [A, B, scales] = rheome.operators.dirac_frame(...)
%
% Faithful standalone port of Brainstorm's tess_operators 'Dirac' variant:
%
%   A = (1-tau)*(D_int^2 / sL)  +  tau*(E / sE) ,     B = kron(Mass, I4)
%
% D_int^2 (rheome.operators.dirac_intrinsic_sq) is the INTRINSIC spin-connection / tangent-
% frame Dirac squared (from the immersion / edge vectors); E (rheome.operators.dirac_extrinsic)
% is the EXTRINSIC Gauss-map (shape-operator) Dirac squared. Each block is co-normalized
% to unit largest generalized eigenvalue vs B (sL, sE) so tau is a DIMENSIONLESS dial
% (tau=0.5 == equal intrinsic/extrinsic frame weight), portable across mesh size/units.
% BOTH blocks couple the quaternion components (the rotating cortical frame) -- neither
% presupposes a tangent/normal split. Build per connected component (hemisphere).
%
% NORMALIZE (default true) -- the block co-normalization is only a MAGNITUDE balance (both
% blocks are already 1/length^2), so dropping it is a pure SCALAR sL on the operator:
%   true  : A = (1-tau)(D^2/sL) + tau(E/sE)         dimensionless, spectrum in [0,1] (Brainstorm)
%   false : A = (1-tau)D^2 + tau(sL/sE)E  ( = sL*A) PHYSICAL units 1/length^2 -> sqrt(lambda) is a
%           wavenumber, calibrated by the mesh's real size. This is EXACTLY sL times the
%           normalized operator, so the EIGENMODES are identical and only the eigenVALUES
%           rescale (lambda_phys = sL*lambda_norm). Use it to read the Dirac spectrum in real mm.
%
% INPUTS:  V,F surface;  tau in [0,1] (0=pure intrinsic, 1=pure extrinsic, 0.5=balanced);
%          N (optional) [V x 3] unit vertex normals for the Gauss map;
%          normalize (optional) true (default) | false (physical eigenvalues, same eigenmodes).
% OUTPUT:  A [4V x 4V] symmetric;  B [4V x 4V] SPD mass;  scales = [sL, sE].
%
% See also: rheome.operators.dirac_intrinsic_sq, rheome.operators.dirac_extrinsic, rheome.eigen.dirac_frame
%
% Author: Diellor Basha, 2026 (port of Brainstorm tess_operators 'Dirac')

    if nargin < 4, N = []; end
    if nargin < 5 || isempty(normalize), normalize = true; end
    if ~(isscalar(tau) && tau >= 0 && tau <= 1)
        error('operators:dirac_frame:tau', 'tau must be a scalar in [0,1].');
    end
    L4 = rheome.operators.dirac_intrinsic_sq(V, F);              % intrinsic Dirac^2 (immersion f)
    E  = rheome.operators.dirac_extrinsic(V, F, N);             % extrinsic Dirac^2 (Gauss map N)
    M  = rheome.operators.mass(V, F, 'galerkin');
    B  = kron(M, speye(4));
    sL = i_lambda_max(L4, B);
    sE = i_lambda_max(E,  B);
    if normalize
        A = (1 - tau) * (L4 / sL) + tau * (E / sE);      % dimensionless (spectrum in [0,1])
    else
        A = (1 - tau) * L4 + tau * (sL / sE) * E;        % = sL*(normalized): PHYSICAL, same eigenmodes
    end
    A  = (A + A') / 2;
    scales = [sL, sE];
end

% Largest generalized eigenvalue of a symmetric PSD pencil (A,B), B SPD -- the block
% co-normalization scale. Factorizes only the well-conditioned mass B; coarse tol is
% enough (only the scale matters).  (Port of tess_operators local_lambda_max.)
function lmax = i_lambda_max(A, B)
    A = (A + A')/2;  B = (B + B')/2;
    opts = struct('tol', 1e-4, 'maxit', 300, 'disp', 0);
    lmax = abs(eigs(A, B, 1, 'largestabs', opts));
    if ~isfinite(lmax) || lmax <= 0
        error('operators:dirac_frame:badScale', 'Could not estimate a positive block-normalization scale.');
    end
end

% Author: Diellor Basha, 2026
