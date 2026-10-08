function file = bases(name, K, Kc, opts)
% IMPORT.BASES  Precompute + cache PER-HEMISPHERE operators & eigenbases for both hemispheres.
%
%   file = rheome.import.bases(name)                    % K = Kc = 400
%   file = rheome.import.bases(name, K, Kc)             % K scalar (LBO axis), Kc connection modes
%   file = rheome.import.bases(name, K, Kc, opts)
%
% Caches, for EACH hemisphere of the dataset cortex (L and R), the operator AND its eigenbasis
% for every operator family the module uses:
%   LBO         : cotan stiffness L + Galerkin mass M   (rheome.operators.laplace_beltrami)  + eigenmodes
%   Connection  : Levi-Civita connection Laplacian C     (rheome.operators.connection_laplacian) + eigenmodes
% The Dirac family is NOT duplicated here -- its whole-cortex eigenbasis is already cached by
% rheome.import.dirac (dirac__*.mat) and is sliced per hemisphere on load (rheome.load.bases attaches it), so a
% single copy serves both hemispheres.
%
% ⭐ THE CACHE IS KEYED BY K. Files are written as bases__K<K>__Kc<Kc>.mat, so a K = 400 and a
% K = 1000 basis COEXIST rather than one silently overwriting the other, and rheome.load.bases VERIFIES
% that the K it loaded is the K that was asked for. The failure mode otherwise is invisible:
% everything downstream still runs, on a different lambda axis, and nothing raises an error.
% (A legacy unkeyed bases.mat is still read if no keyed file exists.)
%
% ⭐ THE CONNECTION SOLVE IS REUSED BY DEFAULT. Kc has nothing to do with K -- the connection
% Laplacian does not depend on the LBO basis at all -- and it is the expensive half. When an
% existing cache already holds the requested Kc, its connection block is COPIED rather than
% re-solved, so raising the LBO K costs only the LBO eigensolve. opts.reuseConn = false forces it.
%
% INPUTS:
%   name  cached dataset name (rheome.import.surface first);  K  LBO modes per hemisphere (default 400);
%   Kc    connection eigenmodes per hemisphere (default = K)
%   opts  .reuseConn  true (default) -- copy the connection block from an existing cache with the
%                     same Kc instead of re-solving it
%        .overwrite   false (default) -- refuse to clobber an existing keyed file
%
% OUTPUT: writes +data/<name>/bases__K<K>__Kc<Kc>.mat holding a struct 'bases':
%   .meta  struct(K, Kc, builtFrom, lambdaMax, sigmaFloor)
%   .hemi  {1 x nH} hemisphere labels (from S.HemiLabel)
%   .(L|R) struct per hemisphere:
%       .gv    [nVh x 1] global vertex indices (local -> whole-cortex)
%       .S     slim hemisphere surface (.Vertices .Faces .VertNormals .nV .nF)
%       .lbo   struct(.L .M  .Phi .Lambda .Mass)      operator + eigenbasis (scalar)
%       .conn  struct(.C  .Psi .Lambda)               operator + eigenbasis (complex tangent)
%
% See also: rheome.load.bases, rheome.import.dirac, rheome.operators.connection_laplacian, rheome.detect.operator
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(K),  K  = 400; end
    if nargin < 3 || isempty(Kc), Kc = K;   end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'reuseConn') || isempty(opts.reuseConn), opts.reuseConn = true;  end
    if ~isfield(opts,'overwrite') || isempty(opts.overwrite), opts.overwrite = false; end

    dsdir = fullfile(rheome.load.root(), char(name));
    sf = fullfile(dsdir, 'surface.mat');
    if ~exist(sf, 'file')
        error('import:bases:noSurface', 'No surface cached for ''%s''. Run rheome.import.surface(''%s'', cortexFile) first.', name, name);
    end
    S = getfield(builtin('load', sf, 'S'), 'S');
    if ~isfield(S,'Hemi') || isempty(S.Hemi)
        error('import:bases:noSplit', 'Surface ''%s'' has no hemisphere split (needs the Structures atlas).', name);
    end
    nH = numel(S.Hemi);

    file = fullfile(dsdir, sprintf('bases__K%d__Kc%d.mat', K, Kc));
    if exist(file, 'file') && ~opts.overwrite
        error('import:bases:exists', ...
            ['%s already exists. Pass opts.overwrite = true to rebuild it, or just load it -- ' ...
             'not re-solving is the point of the cache.'], file);
    end

    % ---- can an existing connection block be reused? ----
    % The connection Laplacian is independent of K and its eigensolve is the expensive half, so
    % raising the LBO K should not pay for it twice.
    prior = [];
    if opts.reuseConn
        cand = dir(fullfile(dsdir, 'bases__K*__Kc*.mat'));
        if isempty(cand) && exist(fullfile(dsdir,'bases.mat'),'file')
            cand = dir(fullfile(dsdir, 'bases.mat'));
        end
        for c = 1:numel(cand)
            try
                p = getfield(builtin('load', fullfile(cand(c).folder, cand(c).name), 'bases'), 'bases');
                if isfield(p,'meta') && isfield(p.meta,'Kc') && p.meta.Kc == Kc
                    prior = p;
                    fprintf('rheome.import.bases[%s]: reusing the connection block (Kc=%d) from %s\n', ...
                        name, Kc, cand(c).name);
                    break;
                end
            catch, continue; end
        end
    end

    fprintf('=== rheome.import.bases[%s] (K=%d, Kc=%d, %d hemispheres) ===\n', name, K, Kc, nH);
    bases = struct('meta', struct('K',K, 'Kc',Kc, 'builtFrom',char(name)), 'hemi', {S.HemiLabel});
    lmaxAll = zeros(1, nH);
    for h = 1:nH
        lab = upper(regexp(S.HemiLabel{h}, '[LR]$', 'match', 'once'));   % 'Cortex L' -> 'L'
        if isempty(lab), lab = sprintf('H%d', h); end
        Sh = rheome.utils.hemisphere(S, [], h);  gv = Sh.GlobalVertices;
        fprintf('  [%s] %d vertices:', lab, Sh.nV);

        t = tic;
        [L, M] = rheome.operators.laplace_beltrami(Sh.Vertices, Sh.Faces, 'galerkin');
        if K >= Sh.nV
            error('import:bases:K', 'K = %d but hemisphere %s has only %d vertices.', K, lab, Sh.nV);
        end
        b = rheome.eigen.modes(L, M, K);
        lmaxAll(h) = max(b.Lambda);
        fprintf(' LBO(K=%d) %.1fs lmax=%.4g vpw=%.1f |', K, toc(t), lmaxAll(h), b.vertsPerHalfWave);

        if ~isempty(prior) && isfield(prior, lab) && isfield(prior.(lab), 'conn')
            conn = prior.(lab).conn;
            fprintf(' connection REUSED (Kc=%d)\n', Kc);
        else
            t = tic;
            C = rheome.operators.connection_laplacian(Sh.Vertices, Sh.Faces);
            [Psi, Lm]   = rheome.eigen.smallest(C.A, C.B, Kc);   % negative-sigma shift (not naive smallestabs)
            [Lambda, o] = sort(real(diag(Lm)));  Psi = Psi(:, o);
            conn = struct('C',C, 'Psi',Psi, 'Lambda',Lambda);
            fprintf(' connection(Kc=%d) %.1fs\n', Kc, toc(t));
        end

        Sslim = struct('Vertices',Sh.Vertices, 'Faces',Sh.Faces, 'VertNormals',Sh.VertNormals, 'nV',Sh.nV, 'nF',Sh.nF);
        bases.(lab) = struct('gv', gv, 'S', Sslim, ...
            'lbo',  struct('L',L, 'M',M, 'Phi',b.Phi, 'Lambda',b.Lambda, 'Mass',b.Mass), ...
            'conn', conn);
    end

    % What this K actually buys, recorded so it is not rederived downstream: the finest scale-space
    % member that responds properly is sigma_m >= 3.08/sqrt(lambda_max), since below that a mexhat
    % member loses over 5% of its mass past the end of the basis.
    bases.meta.lambdaMax  = max(lmaxAll);
    bases.meta.sigmaFloor = 3.08 / sqrt(max(lmaxAll));
    fprintf('lambda_max = %.4g -> frame usable floor sigma_m >= %.1f mm\n', ...
        bases.meta.lambdaMax, 1000*bases.meta.sigmaFloor);

    builtin('save', file, 'bases', '-v7.3');
    fprintf('rheome.import.bases[%s]: per-hemisphere operators + eigenbases -> %s\n', name, file);
end

% Author: Diellor Basha, 2026
