function R = risefall(a, fs, t0, t1, opts)
% DETECT.RISEFALL  Where an amplitude episode switches on and off: a double-sigmoid fit.
%
%   R = rheome.detect.risefall(a, fs, t0, t1)                  % a candidate episode [t0, t1] s
%   R = rheome.detect.risefall(a, fs, t0, t1, Margin=1, Smooth=0.1)
%
%       a(t) = b + h * S((t - tOn)/tauR) * S(-(t - tOff)/tauF),    S(x) = 1/(1 + exp(-x))
%
% fitted to the amplitude a over [t0 - Margin, t1 + Margin]. tOn and tOff are the HALF-rise and
% HALF-fall times -- where the episode is half built and half gone -- and tauR, tauF how fast it
% switches. The fit splits an episode into three phases the rest of the analysis reads separately:
%     build-up  [tOn - 2 tauR, tOn + 2 tauR]    the amplitude is growing (12% -> 88% of h)
%     plateau   [tOn + 2 tauR, tOff - 2 tauF]   built, before it starts to go
%     decay     [tOff - 2 tauF, tOff + 2 tauF]
%
% ⭐ WHY A FIT AND NOT THE THRESHOLD CROSSINGS. A threshold detector's start and end are where the
% envelope crosses a fixed level, so they move with the episode's height and with the noise at the
% flanks (rheome.detect.spindle: a loud 0.2 s burst reported as 0.8 s). The sigmoid centres are properties
% of the SHAPE: a taller episode with the same timing gets the same tOn and tOff.
%
% INPUTS
%   a        [1 x nT] amplitude (e.g. the mean alpha envelope over a region)
%   fs       samples per second          t0, t1   the candidate episode (s), e.g. from rheome.detect.spindle
%   Margin   seconds either side included in the fit (1)
%   Smooth   moving mean applied first, symmetric (0.1 s)
%
% OUTPUT (struct R)
%   .tOn .tOff .tauR .tauF .b .h (seconds and amplitude units)   .r2   .ok
%   .buildup .plateau .decay   [start end] s; plateau is [] when the fit leaves none
%
% See also: rheome.detect.spindle, rheome.detect.tilepath
%
% Author: Diellor Basha, 2026

    arguments
        a (1,:) double
        fs (1,1) double {mustBePositive}
        t0 (1,1) double
        t1 (1,1) double
        opts.Margin (1,1) double {mustBeNonnegative} = 1
        opts.Smooth (1,1) double {mustBeNonnegative} = 0.1
    end
    R = struct('tOn',NaN,'tOff',NaN,'tauR',NaN,'tauF',NaN,'b',NaN,'h',NaN,'r2',NaN,'ok',false, ...
        'buildup',[],'plateau',[],'decay',[]);
    i0 = max(1, round((t0 - opts.Margin) * fs) + 1);  i1 = min(numel(a), round((t1 + opts.Margin) * fs) + 1);
    if i1 - i0 < 10, return; end
    y = movmean(a(i0:i1), max(1, round(opts.Smooth * fs)));  t = ((i0:i1) - 1) / fs;
    % ⚠ FIT ON A UNIT SCALE. A source envelope is ~1e-11 A.m; the optimiser's finite differences and
    % tolerances do not reach that scale and every fit stayed at its starting point (tau = 0.100 exactly).
    sc = max(abs(y));  if ~(sc > 0), return; end
    y = y / sc;
    S = @(x) 1 ./ (1 + exp(-x));
    f = @(p, t) p(1) + p(2) * S((t - p(3)) / p(4)) .* S(-(t - p(5)) / p(6));
    b0 = prctile(y, 5);  h0 = max(y) - b0;
    p0 = [b0 h0 t0 0.1 t1 0.1];
    lb = [min(y) - h0, 0,       t(1), 0.01, t0,   0.01];
    ub = [max(y),      3 * h0,  t1,   1.0,  t(end), 1.0];
    p0 = min(max(p0, lb), ub);
    o = optimoptions('lsqcurvefit', 'Display', 'off');
    try
        p = lsqcurvefit(f, p0, t, y, lb, ub, o);
    catch
        return
    end
    if p(5) <= p(3), return; end
    r = y - f(p, t);                                             % r2 is scale-free
    R.b = p(1) * sc;  R.h = p(2) * sc;  R.tOn = p(3);  R.tauR = p(4);  R.tOff = p(5);  R.tauF = p(6);
    R.r2 = 1 - sum(r.^2) / max(sum((y - mean(y)).^2), realmin);
    R.buildup = [R.tOn - 2*R.tauR, R.tOn + 2*R.tauR];
    R.decay = [R.tOff - 2*R.tauF, R.tOff + 2*R.tauF];
    if R.tOff - 2*R.tauF > R.tOn + 2*R.tauR
        R.plateau = [R.tOn + 2*R.tauR, R.tOff - 2*R.tauF];
    end
    R.ok = true;
end

% Author: Diellor Basha, 2026
