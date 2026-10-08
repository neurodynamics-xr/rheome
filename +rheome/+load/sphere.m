function sphere = sphere(name)
% LOAD.SPHERE  Load the bundled registration-sphere sandbox (operators + eigenmodes).
%
%   sphere = rheome.load.sphere(name)
%
% Returns the single struct that rheome.import.sphere cached to +data/<name>/sphere.mat -- the
% subject's FreeSurfer registration sphere (ico5 @ 100 mm) with EVERYTHING on it:
%   .S .L .M .basis .dirac .conn(.C/.Psi/.Lambda) .cortex .globalVertices .radius .hemi .meta
% One call gets the whole sandbox at cortical resolution. Returns [] if not cached.
%
% See also: rheome.import.sphere, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'sphere.mat');
    if ~exist(f, 'file')
        sphere = [];
        return;
    end
    sphere = getfield(builtin('load', f, 'sphere'), 'sphere');
end

% Author: Diellor Basha, 2026
