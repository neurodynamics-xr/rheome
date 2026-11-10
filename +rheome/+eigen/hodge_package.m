function pkg = hodge_package(sol, H)
% EIGEN.HODGE_PACKAGE  Hodge-Dirac modes as one vertex potential each, merged and ordered by mu: the basis the editor reads.
%
%   pkg = rheome.eigen.hodge_package(rheome.eigen.hodge_dirac(H, nModes), H)      % H = rheome.operators.hodge_dirac(V, F)
%       pkg.mu [K x 1]    the positive Dirac eigenvalues sqrt(lambda), ascending (the two ladders merged, stable)
%       pkg.kind [K x 1]  0 exact, 2 coexact
%       pkg.P [V x K]     single: the vertex potential of each mode, A = psi/sqrt(2) (exact) or S = <chi/sqrt(2), omega_v> (coexact:
%                         the star2-weighted mean of chi over the intrinsic faces round v); the flow of mode k is (grad A + n x grad S)/mu
%       pkg.mass [V x 1]  the 0-form (vertex) mass star0, the inner product the editor reads the modes under
%
% The potentials are rounded to single as nxr-convert stores them (package: float32) and the editor holds them on the device.
% Port of nxr-convert hodge.py _package_split (its W and the 1/sqrt(2) of a unit-normalised ladder).
%
% See also: rheome.eigen.hodge_dirac, rheome.inverse.hodge_weights
%
% Author: Diellor Basha, 2026 (port of nxr-convert hodge.package)
arguments, sol (1,1) struct, H (1,1) struct, end
lam = [sol.exactLambda; sol.coexactLambda];
mu = sqrt(max(lam, 0));
use = find(mu > 1e-6 * max(mu));
[~, o] = sort(mu(use)); o = use(o);                              % stable, as numpy's argsort(kind='stable')
n0 = numel(sol.exactLambda); nV = numel(H.star0); nF = numel(H.star2);
iF = reshape(H.mesh.F, 3, [])';                                  % the intrinsic faces (after the flips)
inc = sparse(iF(:), repmat((1:nF)', 3, 1), 1, nV, nF);
Wf = spdiags(1 ./ max(full(sum(inc, 2)), 1), 0, nV, nV) * inc * spdiags(H.star2, 0, nF, nF);
P = zeros(nV, numel(o));
ex = o <= n0;
P(:, ex) = sol.exactModes(:, o(ex)) / sqrt(2);
P(:, ~ex) = Wf * sol.coexactModes(:, o(~ex) - n0) / sqrt(2);
pkg = struct('mu', mu(o), 'kind', 2 * double(~ex), 'P', single(P), 'mass', H.star0);
end
