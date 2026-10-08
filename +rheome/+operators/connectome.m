function [W, keep] = connectome(V, F, endPts, opts)
% OPERATORS.CONNECTOME  Vertex-resolution structural connectome from fiber endpoints.
%
%   [W, keep] = rheome.operators.connectome(V, F, endPts)
%   [W, keep] = rheome.operators.connectome(V, F, endPts, opts)   % opts.SmoothHops, opts.EdgeThr
%
% Builds a vertex x vertex structural connectome on a cortical mesh from tractography fiber
% endpoints -- faithful port of Brainstorm's tess_connectome vertex branch. Each fiber's two
% endpoints are snapped to their nearest vertices; the raw endpoint graph leaves ~46% of
% vertices isolated (sulcal walls carry no endpoints), so a continuous-kernel MESH SMOOTHING
% W = S^p * Wraw * S^p' (S = row-normalized one-hop mesh averaging) bridges the gaps before the
% largest connected component is returned. The result is the weight matrix consumed by
% rheome.operators.connectome_laplacian and rheome.operators.lb_connectome.
%
% INPUTS:
%   V       [nV x 3] vertices;  F [nF x 3] triangles
%   endPts  [nFib x 2 x 3] the two endpoints of each fiber (rheome.io.read.fibers), in V's coordinates
%   opts    (optional) struct: .SmoothHops (default 3, ~10mm) .EdgeThr (default 1e-4 of max)
%
% OUTPUTS:
%   W       [nV x nV] sparse symmetric connectome weights (zero diagonal)
%   keep    indices of the largest connected component of W>0
%
% See also: rheome.io.read.fibers, rheome.operators.connectome_laplacian, rheome.operators.lb_connectome, rheome.eigen.connectome, rheome.eigen.lift
%
% Author: Diellor Basha, 2026 (port of tess_connectome / local_vertex_connectome)

    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'SmoothHops'), opts.SmoothHops = 3;    end
    if ~isfield(opts,'EdgeThr'),    opts.EdgeThr    = 1e-4; end
    nV = size(V, 1);  F = double(F);

    % raw endpoint connectome: nearest vertex of each endpoint, symmetric fiber counts
    v1 = dsearchn(V, squeeze(endPts(:,1,:)));
    v2 = dsearchn(V, squeeze(endPts(:,2,:)));
    Wraw = sparse([v1;v2], [v2;v1], 1, nV, nV);
    Wraw = Wraw - spdiags(diag(Wraw), 0, nV, nV);

    % continuous-kernel mesh smoothing S^p (row-normalized one-hop averaging), bridges sulci
    ii = [F(:,1);F(:,2);F(:,3)];  jj = [F(:,2);F(:,3);F(:,1)];
    A  = double(sparse([ii;jj], [jj;ii], 1, nV, nV) > 0);
    Sm = A + speye(nV);  Sm = spdiags(1./sum(Sm,2), 0, nV, nV) * Sm;
    Sp = Sm ^ opts.SmoothHops;
    W  = Sp * Wraw * Sp';
    W  = (W + W')/2;  W = W - spdiags(diag(W), 0, nV, nV);
    W(W < opts.EdgeThr * max(W(:))) = 0;

    % largest connected component
    cc = conncomp(graph(W > 0));
    tb = tabulate(cc);  [~, big] = max(tb(:,2));
    keep = find(cc == big);
end

% Author: Diellor Basha, 2026
