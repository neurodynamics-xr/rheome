function arr = linear(varargin)
% SENSORS.LINEAR  A laminar probe: N contacts on a line.
%
%   arr = rheome.sensors.linear()
%   arr = rheome.sensors.linear('NumSites', 32, 'Pitch', 100e-6)
%
% The one geometry with intrinsic dimension 1, and the reason that matters is topological
% rather than numerical: a winding number is summed around a face (a 2-cell), and a chain
% admits no faces, so a winding number cannot be defined on it at all. Rotating and spiral
% patterns are not hard to measure on a probe, they are not measurable -- see rheome.sensors.window,
% which reports this rather than leaving it implicit.
%
% A wave crossing the shank at angle theta also projects: the probe sees only k*cos(theta),
% so any speed inferred from it reads c/cos(theta), biased HIGH and uncorrectable from the
% probe alone.
%
% INPUTS (name-value):
%   'NumSites'  contacts (default 32)
%   'Pitch'     contact spacing, metres (default 100e-6)
%   'Name'      display name (default 'laminar probe')
% OUTPUT:
%   arr  see rheome.sensors.grid for the struct contract
%
% See also: rheome.sensors.grid, rheome.sensors.window, rheome.sensors.graph
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('NumSites', 32);
    p.addParameter('Pitch',    100e-6);
    p.addParameter('Name',     'laminar probe');
    p.parse(varargin{:});
    o = p.Results;

    n = double(o.NumSites);
    if ~isscalar(n) || n < 2 || mod(n, 1) ~= 0
        error('sensors:linear:numSites', 'NumSites must be an integer >= 2, got %s.', mat2str(n));
    end
    if ~(o.Pitch > 0)
        error('sensors:linear:pitch', 'Pitch must be positive, got %g.', o.Pitch);
    end

    z = (0:n-1).' * double(o.Pitch);
    P = [zeros(n, 2), z];
    L = arrayfun(@(i) sprintf('s%03d', i), 1:n, 'UniformOutput', false);
    arr = sen_array(o.Name, 'linear', P, 1, L);
end

% Author: Diellor Basha, 2026
