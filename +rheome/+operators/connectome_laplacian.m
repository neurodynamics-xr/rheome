function [A, B, N] = connectome_laplacian(W, keep)
% OPERATORS.CONNECTOME_LAPLACIAN  Symmetric normalized graph Laplacian of a connectome.
%
%   [A, B, N] = rheome.operators.connectome_laplacian(W, keep)
%
% The whole-brain connectome operator: on the largest connected component (keep) of the
% connectome weights W, the symmetric normalized graph Laplacian
%     A = I - D^{-1/2} W D^{-1/2},   B = I
% whose low eigenmodes are the network's smooth (long-range) harmonics. Faithful port of
% tess_operators / local_build_connectome_operator ('Connectome Laplacian'). Also returns the
% normalized adjacency N = D^{-1/2} W D^{-1/2}, so rheome.eigen.connectome can use the factorization-
% free normalized-adjacency trick (lambda_L = 1 - lambda_N).
%
% INPUTS:
%   W     [nV x nV] connectome weights (rheome.operators.connectome);  keep  its largest-CC indices
%
% OUTPUTS:
%   A     [nk x nk] symmetric normalized Laplacian (nk = numel(keep))
%   B     [nk x nk] identity (mass)
%   N     [nk x nk] normalized adjacency (A = I - N)
%
% See also: rheome.operators.connectome, rheome.operators.lb_connectome, rheome.eigen.connectome
%
% Author: Diellor Basha, 2026 (port of tess_operators 'Connectome Laplacian')

    Wk = W(keep, keep);  nk = numel(keep);
    d  = full(sum(Wk, 2));
    Dm12 = spdiags(1 ./ sqrt(max(d, eps)), 0, nk, nk);
    N = Dm12 * Wk * Dm12;  N = (N + N') / 2;
    A = speye(nk) - N;
    B = speye(nk);
end

% Author: Diellor Basha, 2026
