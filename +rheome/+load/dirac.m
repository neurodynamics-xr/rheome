function [dbasis, basis] = dirac(name, tau, K)
% LOAD.DIRAC  Load a cached relative-Dirac eigenbasis (and Laplace–Beltrami basis) from +data.
%
%   dbasis          = rheome.load.dirac(name)              % tau=0.5, K=400 (defaults)
%   [dbasis, basis] = rheome.load.dirac(name, tau, K)      % also the Laplace–Beltrami basis (out.basis)
%
% Reads +data/<name>/dirac__tau<tau>__K<K>.mat, written by rheome.import.dirac / rheome.import.dataset.
% Errors (with a hint) if that (tau, K) has not been imported for the dataset.
%
% See also: rheome.import.dirac, rheome.eigen.dirac_frame, rheome.load.surface
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(tau), tau = 0.5; end
    if nargin < 3 || isempty(K),   K   = 400; end
    f = fullfile(rheome.load.root(), char(name), sprintf('dirac__tau%.2f__K%d.mat', tau, K));
    if ~exist(f, 'file')
        error('load:dirac:missing', ...
            'No Dirac eigenbasis cached for ''%s'' (tau=%.2f, K=%d). Run  rheome.import.dirac(''%s'', %.2f, %d)  first.', ...
            name, tau, K, name, tau, K);
    end
    C = builtin('load', f);
    dbasis = C.dbasis;
    if nargout > 1
        if isfield(C, 'basis'), basis = C.basis; else, basis = []; end   % older caches may lack it
    end
end

% Author: Diellor Basha, 2026
