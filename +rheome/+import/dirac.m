function file = dirac(name, tau, K)
% IMPORT.DIRAC  Compute + cache the relative-Dirac eigenbasis for a dataset.
%
%   file = rheome.import.dirac(name)              % tau=0.5, K=400 (defaults)
%   file = rheome.import.dirac(name, tau, K)
%
% Builds the relative-Dirac eigenbasis (rheome.eigen.dirac_frame -- the expensive
% per-hemisphere eigensolve, ~50 s) from the dataset's cached surface and writes it to
% +data/<name>/dirac__tau<tau>__K<K>.mat. Requires rheome.import.surface(name, ...) first.
%
% See also: rheome.import.surface, rheome.import.dataset, rheome.eigen.dirac_frame, rheome.load.dirac
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(tau), tau = 0.5; end
    if nargin < 3 || isempty(K),   K   = 400; end

    dsdir = fullfile(rheome.load.root(), char(name));
    sf = fullfile(dsdir, 'surface.mat');
    if ~exist(sf, 'file')
        error('import:dirac:noSurface', ...
            'No surface cached for ''%s''. Run  rheome.import.surface(''%s'', cortexFile)  first.', name, name);
    end
    C = builtin('load', sf, 'S', 'L', 'M');  S = C.S;

    fprintf('rheome.import.dirac[%s]: computing relative Dirac (tau=%.2f, K=%d) -- this is the slow bit...\n', name, tau, K);
    t = tic;
    hemi = [];  if isfield(S,'Hemi'), hemi = S.Hemi; end                      % atlas L/R split (if cached)
    dbasis = rheome.eigen.dirac_frame(S.Vertices, S.Faces, tau, K, S.VertNormals, hemi);   %#ok<NASGU>
    basis  = rheome.eigen.modes(C.L, C.M, K);                    %#ok<NASGU> Laplace–Beltrami (for out.basis: eigenvalue axis)
    file = fullfile(dsdir, sprintf('dirac__tau%.2f__K%d.mat', tau, K));
    builtin('save', file, 'dbasis', 'basis', 'tau', 'K', '-v7.3');
    fprintf('  computed + cached in %.1fs -> %s\n', toc(t), file);
end

% Author: Diellor Basha, 2026
