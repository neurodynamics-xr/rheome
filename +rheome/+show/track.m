function v = track(spec, x)
% SHOW.TRACK  A parameter's course over its support, from an analytic track: the Editor's track, in Rheome.
%
%   v = rheome.show.track(spec, x)        x in [0, 1] of the support (any size); v the same size
%
% The same course the nxr-cortical-flow Editor gives a time atom's parameter (engine surfaceKit trackValue), so a
% plant drawn in a Rheome figure is the plant the Editor made. A TRACK is
%   KEYS     a struct array with .x (0..1), .v and .seg (the segment AFTER the key: 'hold' | 'linear' (default) |
%            'ease' | 'steps' (four dyadic quarters) | 'curve' (monotone cubic, flat at the ends and at a turn)),
%            or an [n x 2] matrix [x v] (linear segments); before the first key its value, after the last its value
%   or START -> TO along a PROFILE: 'ramp' | 'ease' | 'oscillate' (out and back, once) | 'step' (four steps)
%   then an OSCILLATOR on top, .Osc = struct('rate', cycles over the support, 'depth', +- amplitude,
%            'phase', rad (0), 'shape', 'sine' | 'triangle' | 'square')
% spec is a struct with Keys, or Start (default 0) with To and Profile, and Osc (each optional).
%
% Example: a heading turning one full turn, wobbling twice:
%   v = rheome.show.track(struct('Start', 0, 'To', 2*pi, 'Profile', 'ramp', 'Osc', struct('rate', 2, 'depth', 0.3)), linspace(0, 1, 200));
%
% Links between tracks (an offset, a ratio, nose-first, dispersion) are not resolved here: compute the follower's
% values and give them to rheome.show.timetracks as Values.
%
% See also: rheome.show.timetracks
%
% Author: Diellor Basha, 2026
    if ~isstruct(spec) || ~isscalar(spec), error('rheome:show:track', 'spec: a struct (Keys, or Start/To/Profile; Osc).'); end
    known = {'Keys', 'Start', 'To', 'Profile', 'Osc'};
    bad = setdiff(fieldnames(spec), known);
    if ~isempty(bad), error('rheome:show:track', '%s is not a track field (%s).', strjoin(bad, ', '), strjoin(known, ', ')); end
    x = double(x);
    if isfield(spec, 'Keys') && ~isempty(spec.Keys)
        v = arrayfun(@(u) keysAt(keyList(spec.Keys), u), x);
    else
        s = field(spec, 'Start', 0);  v = s + zeros(size(x));
        if isfield(spec, 'To') && ~isempty(spec.To)
            v = s + (spec.To - s) * profileAt(char(field(spec, 'Profile', 'ramp')), x, 0);
        end
    end
    if isfield(spec, 'Osc') && ~isempty(spec.Osc)
        o = spec.Osc;
        a = 2 * pi * o.rate * x + field(o, 'phase', 0);
        switch char(field(o, 'shape', 'sine'))
            case 'sine', w = sin(a);
            case 'triangle', w = (2 / pi) * asin(sin(a));
            case 'square', w = 2 * (sin(a) >= 0) - 1;
            otherwise, error('rheome:show:track', 'Osc.shape: sine, triangle or square.');
        end
        v = v + o.depth * w;
    end
end

function w = profileAt(p, x, phase)
    u = min(max(x, 0), 1);
    switch p
        case 'ramp', w = u;
        case 'ease', w = (1 - cos(pi * u)) / 2;
        case 'oscillate', w = (1 - cos(2 * pi * u + phase)) / 2;
        case 'step', w = min(3, floor(4 * u)) / 3;
        otherwise, error('rheome:show:track', 'Profile: ramp, ease, oscillate or step.');
    end
end

function K = keyList(K)
% keys as a struct array sorted by x, each with a segment
    if isnumeric(K), K = struct('x', num2cell(K(:, 1)), 'v', num2cell(K(:, 2))); end
    if ~isfield(K, 'seg'), [K.seg] = deal('linear'); end
    for i = 1:numel(K), if isempty(K(i).seg), K(i).seg = 'linear'; end, end
    [~, o] = sort([K.x]);  K = K(o);
end

function v = keysAt(K, x)
    if x <= K(1).x, v = K(1).v; return; end
    for q = 1:numel(K) - 1
        a = K(q);  b = K(q + 1);
        if x < b.x
            h = b.x - a.x;  u = 1;  if h > 0, u = (x - a.x) / h; end
            switch char(a.seg)
                case 'curve'   % cubic Hermite with Fritsch-Butland slopes: never overshoots
                    u2 = u * u;  u3 = u2 * u;
                    v = (2*u3 - 3*u2 + 1) * a.v + (u3 - 2*u2 + u) * h * slope(K, q) + (-2*u3 + 3*u2) * b.v + (u3 - u2) * h * slope(K, q + 1);
                    return
                case 'hold', w = 0;
                case 'ease', w = (1 - cos(pi * u)) / 2;
                case 'steps', w = min(3, floor(4 * u)) / 3;
                case 'linear', w = u;
                otherwise, error('rheome:show:track', 'Keys.seg: hold, linear, ease, steps or curve.');
            end
            v = a.v + (b.v - a.v) * w;  return
        end
    end
    v = K(end).v;
end

function m = slope(K, q)
% the weighted harmonic mean of the two secants at key q; zero at the ends and at a turn
    m = 0;
    if q <= 1 || q >= numel(K), return; end
    h0 = K(q).x - K(q-1).x;  h1 = K(q+1).x - K(q).x;
    if h0 <= 0 || h1 <= 0, return; end
    d0 = (K(q).v - K(q-1).v) / h0;  d1 = (K(q+1).v - K(q).v) / h1;
    if d0 * d1 <= 0, return; end
    m = 3 * (h0 + h1) / ((2*h1 + h0) / d0 + (h1 + 2*h0) / d1);
end

function v = field(s, n, d)
    if isfield(s, n) && ~isempty(s.(n)), v = s.(n); else, v = d; end
end
