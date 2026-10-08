function basis = connectome(N, K)
% EIGEN.CONNECTOME  Low eigenmodes of a normalized connectome Laplacian (normalized-adjacency trick).
%
%   basis = rheome.eigen.connectome(N, K)     % N = normalized adjacency from rheome.operators.connectome_laplacian
%
% The smooth network harmonics of the symmetric normalized connectome Laplacian L = I - N. Rather
% than shift-invert the (near-singular) L for its smallest eigenvalues, take the LARGEST of the
% normalized adjacency N -- factorization-free, well-conditioned -- and map back:
%     [Phi, mu] = eigs(N, K, 'largestreal');   lambda_L = 1 - mu.
% Faithful port of tess_eigen / proto_fiber_connectome's connectome branch. (For the LB-Connectome
% variant use rheome.eigen.modes(A,B,K) directly -- a generalized symmetric pencil; the Dirac-Connectome
% basis is that LB-Connectome basis lifted by eigen.dirac.)
%
% INPUTS:
%   N   [nk x nk] normalized adjacency D^{-1/2} W D^{-1/2} (rheome.operators.connectome_laplacian)
%   K   number of low modes (default 250)
%
% OUTPUT (struct basis): .Phi [nk x K] (orthonormal) .Lambda [K x 1] (ascending) .Mass (I) .nV .K
%
% See also: rheome.operators.connectome_laplacian, rheome.operators.connectome, rheome.eigen.modes, eigen.dirac
%
% Author: Diellor Basha, 2026 (port of tess_eigen 'Connectome Laplacian')

    if nargin < 2 || isempty(K), K = 250; end
    N = (N + N') / 2;  nk = size(N, 1);
    opts = struct('maxit', 400);
    [Phi, Ln] = eigs(N, K, 'largestreal', opts);      % largest of N = lowest of L
    lam = 1 - real(diag(Ln));
    [lam, o] = sort(lam);  Phi = real(Phi(:, o));
    basis = struct('Phi', Phi, 'Lambda', lam, 'Mass', speye(nk), 'nV', nk, 'K', K);
end

% Author: Diellor Basha, 2026
