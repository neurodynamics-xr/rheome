function file = connectome(name, Kc, opts)
% IMPORT.CONNECTOME  Build + cache the whole-brain structural-connectome operators & eigenbases.
%
%   file = rheome.import.connectome(name)              % Kc = 400, SmoothHops = 3, EdgeThr = 1e-4
%   file = rheome.import.connectome(name, Kc, opts)    % opts.SmoothHops, opts.EdgeThr
%
% Reads the cached (registered) fiber endpoints (rheome.io.read.fibers), builds the vertex connectome W
% (rheome.operators.connectome), assembles the two whole-brain connectome operators and their
% eigenbases, and caches everything to +data/<name>/connectome.mat:
%   'Connectome Laplacian'  A = I - D^-1/2 W D^-1/2, B = I   -- eigen = normalized-adjacency trick
%   'LB-Connectome'         A = K + gamma*(D-W),   B = M     -- eigen = generalized eigs(A,B)
% The Dirac-Connectome vector basis is NOT cached (it is the LB-Connectome basis lifted on demand
% by rheome.eigen.lift). See rheome.load.connectome.
%
% NUMERICAL NOTE: the LB-Connectome pencil spans ~14 orders of magnitude (near-null low mode),
% so it MUST be solved with the negative-sigma shift (rheome.eigen.smallest = bst_eigs_smallest) -- the
% naive eigs 'smallestabs' returns wrong low modes. Via rheome.eigen.modes it is now correct and
% reproducible (conn.lb.Reproducible = true), as is the Connectome-Laplacian basis.
%
% See also: rheome.io.read.fibers, rheome.operators.connectome, rheome.operators.connectome_laplacian,
%           rheome.operators.lb_connectome, rheome.eigen.connectome, rheome.eigen.modes, rheome.load.connectome
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(Kc), Kc = 400; end
    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'SmoothHops'), opts.SmoothHops = 3;    end
    if ~isfield(opts,'EdgeThr'),    opts.EdgeThr    = 1e-4; end

    dsdir = fullfile(rheome.load.root(), char(name));
    S = getfield(builtin('load', fullfile(dsdir,'surface.mat'), 'S'), 'S');
    endPts = rheome.io.read.fibers(name);
    fprintf('=== rheome.import.connectome[%s] (Kc=%d, %d fibers) ===\n', name, Kc, size(endPts,1));

    t = tic;
    [W, keep] = rheome.operators.connectome(S.Vertices, S.Faces, endPts, opts);
    fprintf('  connectome W: %d-vtx CC, %d edges (%.1fs)\n', numel(keep), nnz(triu(W)), toc(t));

    % --- Connectome Laplacian (well-conditioned, reproducible) ---
    t = tic;
    [A_cl, B_cl, N] = rheome.operators.connectome_laplacian(W, keep);
    bCL = rheome.eigen.connectome(N, Kc);
    fprintf('  Connectome-Laplacian + eigen(K=%d): lambda[%.4f..%.4f] (%.1fs)\n', Kc, bCL.Lambda(1), bCL.Lambda(end), toc(t));

    % --- LB-Connectome (ill-conditioned eigenbasis -- one sample) ---
    t = tic;
    [A_lb, B_lb, gamma] = rheome.operators.lb_connectome(S, W, keep);
    bLB = rheome.eigen.modes(A_lb, B_lb, Kc);   % rheome.eigen.smallest (negative-sigma shift) -> true, reproducible eigenpairs
    fprintf('  LB-Connectome (gamma=%.4g) + eigen(K=%d): lambda[%.3g..%.3g] (%.1fs)\n', gamma, Kc, bLB.Lambda(1), bLB.Lambda(end), toc(t));

    conn = struct();
    conn.meta = struct('Kc',Kc, 'SmoothHops',opts.SmoothHops, 'EdgeThr',opts.EdgeThr, ...
                       'gamma',gamma, 'builtFrom',char(name));
    conn.W    = W;
    conn.keep = keep;
    conn.laplacian = struct('A',A_cl, 'B',B_cl, 'N',N, 'Phi',bCL.Phi, 'Lambda',bCL.Lambda, 'Reproducible',true);
    conn.lb        = struct('A',A_lb, 'B',B_lb, 'gamma',gamma, 'Phi',bLB.Phi, 'Lambda',bLB.Lambda, 'Reproducible',true);  %#ok<STRNU> saved by name below

    file = fullfile(dsdir, 'connectome.mat');
    builtin('save', file, 'conn', '-v7.3');
    fprintf('rheome.import.connectome[%s]: W + operators + eigenbases -> %s\n', name, file);
end

% Author: Diellor Basha, 2026
