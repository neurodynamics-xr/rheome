function [U, info] = generate(G, family, varargin)
% SENSORS.GENERATE  A space-time pattern on a sensor array, from a dynamical family.
%
%   U = rheome.sensors.generate(G, 'wave')
%   U = rheome.sensors.generate(G, 'wave', 'Vertex', 45, 'Args', {30})
%   [U, info] = rheome.sensors.generate(G, 'travwave', 'SamplingFrequency', 200, 'SignalLength', 128)
%
% Seeds a delta at one sensor and evolves it with one of the joint dynamical families.
% The FAMILY sets the dynamics, the SEED sets where, and the array's own graph sets what
% the pattern can be. A line of contacts is not a special case here -- it is a graph with
% a different spectrum, and whatever the family produces on it is what that instrument
% can carry.
%
% ⭐ THIS IS AN ADAPTER, NOT A METHOD. The work is rheome.filters.impulse; this function builds
% the basis from a rheome.sensors.graph, evaluates the family on its NATIVE axis, and picks the
% domain so a caller does not have to know which families are propagators and which are
% joint kernels. Never reimplement the propagation here.
%
% FAMILIES, by native axis:
%   'ts'  propagators g(lambda, t):      diffusion  wave  dampedwave  kleingordon
%   'js'  joint kernels g(lambda, omega): travwave  resonator  gabor  stmatern
%
% ⚠ ARGS ARE IN FAMILY UNITS, NOT ARRAY UNITS. Each family's own defaults apply when Args
% is omitted, and they are usually wrong for a given array -- rheome.filters.diffusion's rate
% enters as exp(-tau*(lambda/lmax)*t), so a tau that suits one record leaves another
% looking like a delta that never moved. info reports what was used; check it rather than
% assuming a pattern evolved. rheome.sensors.calibrate is what carries metres.
%
% INPUTS:
%   G        a rheome.sensors.graph struct
%   family   one of the eight names above
%   'Vertex'             seed sensor (default: the one nearest the array centroid)
%   'SamplingFrequency'  Hz (default 200)
%   'SignalLength'       samples (default 128)
%   'Modes'              a precomputed rheome.sensors.modes struct (default: computed here)
%   'Args'               cell of trailing arguments for the family (default {})
% OUTPUT:
%   U     [nV x nT] the pattern
%   info  .Family .Domain .Vertex .Time .Args .Lmax
%
% See also: rheome.filters.impulse, rheome.sensors.graph, rheome.sensors.modes, rheome.sensors.calibrate
%
% Author: Diellor Basha, 2026

    TS = {'diffusion','wave','dampedwave','kleingordon'};
    JS = {'travwave','resonator','gabor','stmatern'};

    if nargin < 2 || ~(ischar(family) || isstring(family))
        error('sensors:generate:family', 'A family name is required.');
    end
    family = char(family);
    if any(strcmp(family, TS)),     domain = 'ts';
    elseif any(strcmp(family, JS)), domain = 'js';
    else
        error('sensors:generate:family', ...
            ['Unknown family ''%s''. Propagators (lambda x time): %s. Joint kernels ' ...
             '(lambda x omega): %s.'], family, strjoin(TS, ', '), strjoin(JS, ', '));
    end

    p = inputParser;
    p.addParameter('Vertex',            []);
    p.addParameter('SamplingFrequency', 200);
    p.addParameter('SignalLength',      128);
    p.addParameter('Modes',             []);
    p.addParameter('Args',              {});
    p.parse(varargin{:});
    o = p.Results;

    fs = o.SamplingFrequency;  N = o.SignalLength;
    args = o.Args;
    if ~iscell(args), args = {args}; end

    seed = o.Vertex;
    if isempty(seed)
        % The centre-most sensor: an interior seed is the one whose response is the
        % family's own, not the boundary's.
        c = mean(G.Vertices, 1);
        [~, seed] = min(vecnorm(G.Vertices - c, 2, 2));
    end
    if ~isscalar(seed) || seed < 1 || seed > G.nV || mod(seed,1) ~= 0
        error('sensors:generate:vertex', ...
            'Vertex must be an integer in 1..%d, got %s.', G.nV, mat2str(seed));
    end

    M = o.Modes;
    if isempty(M), M = rheome.sensors.modes(G); end
    if isempty(M.Phi)
        error('sensors:generate:modes', ...
            ['Generating a pattern needs the exact eigenbasis (M.Phi). The chebyshev ' ...
             'route returns none -- call rheome.sensors.modes(G) with the default ''eigen''.']);
    end

    t = (0:N-1) / fs;                       % s, the propagator axis
    f = (0:N-1) * (fs / N);                 % Hz, the joint-kernel axis

    switch domain
        case 'ts', Gk = feval("rheome.filters." + family, M.Lambda, t, args{:});
        case 'js', Gk = feval("rheome.filters." + family, M.Lambda, f, args{:});
    end

    basis = struct('Phi', M.Phi, 'Lambda', M.Lambda, ...
                   'Mass', speye(G.nV), 'nV', G.nV);
    if strcmp(domain, 'ts')
        U = rheome.filters.impulse(basis, seed, Gk, 'ts');
    else
        U = rheome.filters.impulse(basis, seed, Gk, 'js', N);
    end

    info = struct('Family', family, 'Domain', domain, 'Vertex', seed, ...
                  'Time', t, 'Args', {args}, 'Lmax', M.Lmax);
end

% Author: Diellor Basha, 2026
