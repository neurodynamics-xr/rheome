function A = atlas(name, which)
% LOAD.ATLAS  Load the cached parcellations (the ROI axis) from +data.
%
%   A = rheome.load.atlas(name)                  % all atlases, [1 x nAtlas]
%   A = rheome.load.atlas(name, 'Desikan-Killiany')   % one, by name
%
% See rheome.io.read.atlas for the fields, and for the three traps the reader normalises away
% (row-vector vertices, Label-not-Name, partial coverage).
%
% See also: rheome.import.atlas, rheome.io.read.atlas, rheome.load.dataset
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'atlas.mat');
    if ~exist(f, 'file')
        error('load:atlas:missing', ...
            'No atlas cached for ''%s''. Run  rheome.import.atlas(''%s'')  first.', name, name);
    end
    C = builtin('load', f, 'atlas');
    A = C.atlas;

    if nargin < 2 || isempty(which), return; end

    k = find(strcmpi({A.Name}, char(which)), 1);
    if isempty(k)
        error('load:atlas:noSuchAtlas', ...
            'No atlas ''%s'' in ''%s''. Available: %s', ...
            char(which), name, strjoin({A.Name}, ', '));
    end
    A = A(k);
end

% Author: Diellor Basha, 2026
