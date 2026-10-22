function [mu, X] = spin_dirac(V, F, op, k)
% EIGEN.SPIN_DIRAC  Eigenpairs of a spin Dirac operator on vertex spinors (Ye et al. quotient).
%
%   [mu, X] = rheome.eigen.spin_dirac(V, F, op, k)     op: 'extrinsic' | 'intrinsicHY' | 'intrinsicYe'
%
% Ye et al.'s quotient (dirac_eigen.m): A' D Mf^-1 D A psi = mu A' D A psi on VERTEX spinors psi
% carried to faces by A (rheome.operators.vertex_to_face; covariant for the intrinsic operators,
% plain for the extrinsic), Mf = face area (x) I4. Solved as K2 psi = nu K1 psi, nu = 1/mu, for the
% k largest |nu| (the k eigenvalues of D nearest 0), K2 = A'DA, K1 = A'D Mf^-1 D A + eps A'Mf A.
% A zero mode (the extrinsic's constants) is 0/0 and drops out. X is normalised to psi' A'Mf A psi = 1.
% mu is sorted ascending. Exact on the unit sphere: intrinsic +-(j+1), extrinsic j and -(j+2).
% Port of cortical-flow solve-spin-diracs.py. Small meshes on the hub; large solves belong to nsp.
%
% Author: Diellor Basha, 2026
arguments, V (:,3) double, F (:,3) double, op (1,1) string, k (1,1) double = 128, end
S = rheome.operators.spin_dirac(V, F);
D = S.(op);
A = rheome.operators.vertex_to_face(V, F, op ~= "extrinsic");
a = repelem(vecnorm(cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:), 2), 2, 2)/2, 4);
M = spdiags(a, 0, numel(a), numel(a)); Mi = spdiags(1./a, 0, numel(a), numel(a));
sym = @(X) (X + X')/2;
K2 = sym(A'*D*A); K1 = sym(A'*D*Mi*D*A); G = sym(A'*M*A);
epsr = 1e-9 * sum(diag(K1)) / sum(diag(G));
[X, nu] = eigs(K2, K1 + epsr*G, k, 'largestabs'); nu = diag(nu);
mu = 1 ./ nu; [mu, o] = sort(mu); X = X(:,o);
X = X ./ sqrt(sum(X .* (G*X), 1));
end
