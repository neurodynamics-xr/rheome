function sol = hodge_dirac(H, nModes)
% EIGEN.HODGE_DIRAC  The Hodge-Dirac eigenmodes of a closed surface, solved split into its exact and coexact ladders and cut at one lambda.
%
%   sol = rheome.eigen.hodge_dirac(H, nModes)       % H = rheome.operators.hodge_dirac(V, F)
%       sol.exactLambda [k0 x 1], sol.exactModes [nV x k0]      psi: (d0' star1 d0) psi = lambda star0 psi, star0-orthonormal
%       sol.coexactLambda [k2 x 1], sol.coexactModes [nF x k2]  chi: (star2 d1 star1^-1 d1' star2) chi = lambda star2 chi, star2-orthonormal
%       sol.mu [2(k0+k2) x 1]  the Dirac eigenvalues +-sqrt(lambda), ascending;  sol.mass  H.M
%
% On a closed surface every nonzero Hodge-Dirac pair (K, M) is exact or coexact: an exact mode is (psi, d0 psi / mu, 0)/sqrt(2),
% a coexact one (0, delta chi / mu, chi)/sqrt(2) with delta = star1^-1 d1' star2, mu = sqrt(lambda). So the pencil is solved as two
% Laplacian ladders: about nModes/2 + 1 pairs each (half and half by Weyl's law), the constants (lambda <= 1e-6 max) out, both cut at
% ONE lambda (the smaller ladder's top) so the band is complete. Shift-invert just below 0 (sigma = -1e-3), as nxr-convert.
%
% Port of nxr-convert eigen/hodge_dirac.py solve_hodge (the editor's make-hodge-package split); parity in tests/tHodgeDiracSolveParity.
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);
%   sol = rheome.eigen.hodge_dirac(rheome.operators.hodge_dirac(V, F), 30);
%
% See also: rheome.operators.hodge_dirac, rheome.geom.glue
%
% Author: Diellor Basha, 2026 (port of nxr-convert solve_hodge)
arguments, H (1,1) struct, nModes (1,1) double {mustBePositive, mustBeInteger}, end
s0 = H.star0; s1 = H.star1; s2 = H.star2; nV = numel(s0); nF = numel(s2);
want = floor(nModes/2) + 1;
L0 = H.d0' * spdiags(s1, 0, numel(s1), numel(s1)) * H.d0;
[lam0, psi] = i_ladder(L0, s0, min(want + 1, nV - 2));
B2 = spdiags(s2, 0, nF, nF) * H.d1;
S2 = B2 * spdiags(1 ./ s1, 0, numel(s1), numel(s1)) * B2';
[lam2, chi] = i_ladder(S2, s2, min(want + 1, nF - 2));
k0 = lam0 > 1e-6 * max(lam0);  k2 = lam2 > 1e-6 * max(lam2);    % the constants (mu = 0) out
top = min(max(lam0(k0)), max(lam2(k2)));                         % ONE lambda for both ladders: the band is complete
k0 = k0 & lam0 <= top;  k2 = k2 & lam2 <= top;
mu = sqrt([lam0(k0); lam2(k2)]);
sol = struct('exactLambda', lam0(k0), 'exactModes', psi(:,k0), 'coexactLambda', lam2(k2), 'coexactModes', chi(:,k2), ...
    'mu', sort([-mu; mu]), 'mass', H.M);
end

function [lam, X] = i_ladder(S, Mv, k)
% the k smallest eigenpairs of (S, diag Mv), M-orthonormal, ascending
n = numel(Mv); S = (S + S') / 2;
[X, D] = eigs(S, spdiags(Mv, 0, n, n), k, -1e-3, 'IsSymmetricDefinite', true);
[lam, o] = sort(diag(D)); lam = max(lam, 0); X = X(:,o);
end
