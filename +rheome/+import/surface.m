function file = surface(name, cortexFile)
% IMPORT.SURFACE  Read a Brainstorm cortex + compute LBO operators, cache to +data.
%
%   file = rheome.import.surface(name, cortexFile)
%
% Reads the surface (rheome.io.read.surface), computes the cotan Laplace-Beltrami L and mass M
% (rheome.operators.laplace_beltrami), and writes +data/<name>/surface.mat. Cheap, but caches
% the parsed surface so later imports/loads never touch the Brainstorm file again.
%
% See also: rheome.import.dirac, rheome.import.study, rheome.import.dataset, rheome.load.surface
%
% Author: Diellor Basha, 2026

    dsdir = fullfile(rheome.load.root(), char(name));
    if ~exist(dsdir, 'dir'), mkdir(dsdir); end

    fprintf('rheome.import.surface[%s]: reading %s\n', name, cortexFile);
    S = rheome.io.read.surface(cortexFile);
    [L, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces, 'galerkin');   %#ok<ASGLU>

    file = fullfile(dsdir, 'surface.mat');
    builtin('save', file, 'S', 'L', 'M', 'cortexFile', '-v7.3');
    fprintf('  cached %d-vertex surface + LBO operators -> %s\n', S.nV, file);
end

% Author: Diellor Basha, 2026
