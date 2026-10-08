function [A, B, gamma] = lb_connectome(S, W, keep)
% OPERATORS.LB_CONNECTOME  Laplace-Beltrami bridged by the structural connectome (whole brain).
%
%   [A, B, gamma] = rheome.operators.lb_connectome(S, W, keep)
%
% The combined whole-brain operator: the per-hemisphere cotan Laplace-Beltrami stiffness K
% (local geometry) bridged across hemispheres by the connectome graph Laplacian D-W (long-range
% fiber connectivity),
%     A = K + gamma*(D - W),   B = M            (M = LBO Galerkin mass)
% with gamma chosen scale-portably so the connectome term sits at the K diagonal scale
%     gamma = mean|diag K| / mean|diag(D-W)|.
% Low eigenmodes (rheome.eigen.modes(A,B,K)) respect BOTH cortical geometry and the connectome; lifting
% them to quaternions (rheome.eigen.lift) gives the Dirac-Connectome vector basis. Faithful port of
% tess_operators / local_build_connectome_operator ('LB-Connectome').
%
% INPUTS:
%   S     surface struct (.Vertices .Faces .nV .Hemi) -- the hemisphere split builds block-diagonal K,M
%   W     [nV x nV] connectome weights (rheome.operators.connectome);  keep  its largest-CC indices
%
% OUTPUTS:
%   A     [nV x nV] combined stiffness K + gamma*(D-W)
%   B     [nV x nV] Galerkin mass M (block-diagonal per hemisphere)
%   gamma scalar connectome/LBO balance factor
%
% See also: rheome.operators.connectome, rheome.operators.laplace_beltrami, rheome.eigen.modes, eigen.dirac
%
% Author: Diellor Basha, 2026 (port of tess_operators 'LB-Connectome')

    nV = S.nV;
    if ~isfield(S,'Hemi') || isempty(S.Hemi)
        error('operators:lb_connectome:noSplit', 'Surface needs a hemisphere split (S.Hemi) for the block-diagonal LBO.');
    end
    nH = numel(S.Hemi);

    % per-hemisphere cotan LBO assembled block-diagonal over the whole brain
    K = sparse(nV, nV);  M = sparse(nV, nV);
    for h = 1:nH
        Sh = rheome.utils.hemisphere(S, [], h);  g = Sh.GlobalVertices;
        [Lh, Mh] = rheome.operators.laplace_beltrami(Sh.Vertices, Sh.Faces, 'galerkin');
        K(g, g) = Lh;  M(g, g) = Mh;   %#ok<SPRIX>
    end

    % connectome graph Laplacian D - W on the connected component
    Wk = W(keep, keep);  nk = numel(keep);  d = full(sum(Wk, 2));
    Lc = sparse(nV, nV);  Lc(keep, keep) = spdiags(d, 0, nk, nk) - Wk;

    % scale-portable balance: match mean nonzero diagonal magnitude
    dK = full(diag(K));  dL = full(diag(Lc));
    gamma = mean(abs(dK(dK ~= 0))) / max(mean(abs(dL(dL ~= 0))), eps);

    A = K + gamma * Lc;  A = (A + A') / 2;
    B = M;
end

% Author: Diellor Basha, 2026
