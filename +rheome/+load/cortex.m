function cortex = cortex(name)
% LOAD.CORTEX  Load the bundled folded-cortex sandbox (operators + eigenmodes).
%
%   cortex = rheome.load.cortex(name)
%
% Returns the single struct that rheome.import.cortex cached to +data/<name>/cortex.mat -- the
% subject's folded cortical hemisphere with EVERYTHING on it:
%   .S .L .M .basis .dirac .conn(.C/.Psi/.Lambda) .globalVertices .hemi .meta
% The folded-metric twin of rheome.load.sphere (same mesh, real cortical embedding). Returns [] if
% not cached. Named rheome.load.cortex (not load.eigen) to avoid shadowing the rheome.eigen.* package.
%
% See also: rheome.import.cortex, rheome.load.sphere, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'cortex.mat');
    if ~exist(f, 'file')
        cortex = [];
        return;
    end
    cortex = getfield(builtin('load', f, 'cortex'), 'cortex');
end

% Author: Diellor Basha, 2026
