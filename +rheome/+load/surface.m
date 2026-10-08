function [S, L, M] = surface(name)
% LOAD.SURFACE  Load a cached cortex surface (and its LBO operators) from +data.
%
%   S            = rheome.load.surface(name)
%   [S, L, M]    = rheome.load.surface(name)     % also the cotan Laplace-Beltrami L and mass M
%
% Reads +data/<name>/surface.mat, written by rheome.import.surface / rheome.import.dataset. Errors
% (with a hint) if the dataset has not been imported.
%
% See also: rheome.import.surface, rheome.load.dirac, rheome.load.study
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'surface.mat');
    if ~exist(f, 'file')
        error('load:surface:missing', ...
            'No surface cached for ''%s''. Run  rheome.import.surface(''%s'', cortexFile)  first.', name, name);
    end
    C = builtin('load', f);
    S = C.S;
    if nargout > 1, L = C.L; end
    if nargout > 2, M = C.M; end
end

% Author: Diellor Basha, 2026
