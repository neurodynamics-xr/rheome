function arr = grid(varargin)
% SENSORS.GRID  A planar lattice: the Utah and ECoG geometries.
%
%   arr = rheome.sensors.grid('utah')                        % 10 x 10 at 400 um
%   arr = rheome.sensors.grid('ecog')                        % 8 x 8 at 10 mm
%   arr = rheome.sensors.grid('Size', [10 10], 'Pitch', 400e-6)
%
% The two presets are the same geometry class at pitches two orders of magnitude apart, and
% they therefore have disjoint measurement windows -- which is the cleanest demonstration in
% the suite that a window is a property of the instrument, not of the method.
%
% ⚠ APERTURE IS THE DIAGONAL. rheome.sensors.window reports 5.1 mm for the Utah preset, not the
% 3.6 mm side.
%
% INPUTS (name-value, or a preset name as the sole argument):
%   'Size'   [rows cols] (default [10 10])
%   'Pitch'  spacing, metres (default 400e-6)
%   'Name'   display name
% OUTPUT (the array contract, shared by every rheome.sensors.* constructor):
%   .Name .Kind .Pos [N x 3] metres .Dim .Pitch .Aperture .Labels .nCh
%
% See also: rheome.sensors.linear, rheome.sensors.meg, rheome.sensors.eeg, rheome.sensors.window
%
% Author: Diellor Basha, 2026

    presets = struct('utah', struct('Size', [10 10], 'Pitch', 400e-6, 'Name', 'Utah 10x10'), ...
                     'ecog', struct('Size', [ 8  8], 'Pitch',  10e-3, 'Name', 'ECoG 8x8'));

    if nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1})) && ...
            ~any(strcmpi(varargin{1}, {'Size', 'Pitch', 'Name'}))
        key = lower(char(varargin{1}));
        if ~isfield(presets, key)
            error('sensors:grid:preset', ...
                'Unknown preset ''%s''. Known: %s.', key, strjoin(fieldnames(presets).', ', '));
        end
        d = presets.(key);
        varargin = [{'Size', d.Size, 'Pitch', d.Pitch, 'Name', d.Name}, varargin(2:end)];
    end

    p = inputParser;
    p.addParameter('Size',  [10 10]);
    p.addParameter('Pitch', 400e-6);
    p.addParameter('Name',  'planar grid');
    p.parse(varargin{:});
    o = p.Results;

    sz = double(o.Size(:).');
    if numel(sz) ~= 2 || any(sz < 2) || any(mod(sz, 1) ~= 0)
        error('sensors:grid:size', 'Size must be two integers >= 2, got %s.', mat2str(sz));
    end
    if ~(o.Pitch > 0)
        error('sensors:grid:pitch', 'Pitch must be positive, got %g.', o.Pitch);
    end

    h = double(o.Pitch);
    [X, Y] = ndgrid((0:sz(1)-1)*h, (0:sz(2)-1)*h);      % rows vary fastest
    P = [X(:), Y(:), zeros(numel(X), 1)];
    [ri, ci] = ndgrid(1:sz(1), 1:sz(2));
    L = arrayfun(@(r, c) sprintf('r%02dc%02d', r, c), ri(:).', ci(:).', 'UniformOutput', false);
    arr = sen_array(o.Name, 'grid', P, 2, L);
end

% Author: Diellor Basha, 2026
