function [L, M] = laplace_beltrami(V, F, massType)
% OPERATORS.LAPLACE_BELTRAMI  Cotangent Laplace-Beltrami stiffness (and mass).
%
%   [L, M] = rheome.operators.laplace_beltrami(V, F)
%   [L, M] = rheome.operators.laplace_beltrami(V, F, massType)   % 'galerkin' (default) | 'lumped'
%
% Builds the discrete Laplace-Beltrami operator using the classic COTANGENT weights
% -- the same "cotan Laplacian" Brainstorm computes (tess_operators ->
% nxr_compute('operators', h, 'laplacian', 'cotan')), reproduced here in pure MATLAB
% so the demo needs no compiled plugin.
%
% CONVENTION (positive-semidefinite "stiffness" / Poisson form):
%   For an interior edge (i,j) with the two opposite triangle angles alpha, beta,
%   the edge weight is  w_ij = (cot(alpha) + cot(beta)) / 2.
%   The stiffness is    L = D - W,  where W is the symmetric weight matrix and
%   D = diag(sum(W,2)). Hence L is symmetric PSD, each row sums to zero, the
%   constant vector is the null mode (lambda_0 = 0), and the smooth eigenmodes
%   solve the generalized problem  L phi = lambda M phi.
%
% Sign note: this returns the PSD stiffness (positive diagonal, negative
% off-diagonal), matching Brainstorm's cotanL used as the operator A in
% eigs(A, M, k, 'smallestabs'). This is the negative of the "geometric Laplacian"
% div-grad; the eigenvalues lambda >= 0 are the ones the atom filters act on.
%
% INPUTS:
%   V        [nV x 3] vertex coordinates
%   F        [nF x 3] triangle vertex indices (1-based)
%   massType mass matrix pairing for the second output: 'galerkin' (default) | 'lumped'
%
% OUTPUTS:
%   L  [nV x nV] sparse symmetric PSD cotan stiffness
%   M  [nV x nV] sparse mass matrix (companion for the eigenproblem L phi = lambda M phi)
%
% Example:
%   [V, F] = rheome.geom.icosphere(4);
%   [L, M] = rheome.operators.laplace_beltrami(V, F);
%   z = V(:, 3);                                       % an l = 1 harmonic: eigenvalue l(l+1) = 2
%   assert(norm(L * ones(size(z))) < 1e-10 && abs((z' * L * z) / (z' * M * z) - 2) < 0.01)
%
% See also: rheome.operators.mass, rheome.eigen.modes
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(massType), massType = 'galerkin'; end
    nV = size(V, 1);

    % --- per-corner cotangents ---
    % For the corner at vertex c the two emanating edges are u, v; with triangle
    % area A, cot(angle at c) = dot(u,v) / (2A). Each corner's cotangent is the
    % weight contribution to the edge OPPOSITE that corner.
    i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);
    p1 = V(i1,:); p2 = V(i2,:); p3 = V(i3,:);

    % edge vectors used per corner
    u1 = p2 - p1;  v1 = p3 - p1;   % corner 1 -> opposite edge (2,3)
    u2 = p3 - p2;  v2 = p1 - p2;   % corner 2 -> opposite edge (3,1)
    u3 = p1 - p3;  v3 = p2 - p3;   % corner 3 -> opposite edge (1,2)

    twoA = sqrt(sum(cross(u1, v1, 2).^2, 2));   % 2 * triangle area (same from any corner)
    twoA = max(twoA, eps);                       % guard degenerate faces

    cot1 = sum(u1 .* v1, 2) ./ twoA;
    cot2 = sum(u2 .* v2, 2) ./ twoA;
    cot3 = sum(u3 .* v3, 2) ./ twoA;

    % --- accumulate the symmetric edge-weight matrix W ---
    % corner 1 cot -> edge (i2,i3); corner 2 -> edge (i3,i1); corner 3 -> edge (i1,i2).
    % 0.5 factor: each interior edge is visited from both adjacent triangles, and the
    % (cot alpha + cot beta)/2 definition carries the 1/2.
    I = [i2; i3; i1];
    J = [i3; i1; i2];
    Wc = 0.5 * [cot1; cot2; cot3];
    W = sparse([I; J], [J; I], [Wc; Wc], nV, nV);   % symmetric off-diagonal weights

    d = full(sum(W, 2));
    L = spdiags(d, 0, nV, nV) - W;
    L = (L + L') / 2;   % enforce exact symmetry

    if nargout > 1
        M = rheome.operators.mass(V, F, massType);
    end
end

% Author: Diellor Basha, 2026
