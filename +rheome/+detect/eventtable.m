function T = eventtable(Ev, Y, A, fs, opts)
% DETECT.EVENTTABLE  Shape-based boundaries of each event: a double-sigmoid on the event's own envelope.
%
%   T = rheome.detect.eventtable(Ev, Y, A, fs)
%   T = rheome.detect.eventtable(Ev, Y, A, fs, Margin=1, Smooth=0.1, TauMax=1, R2Min=0.7, T0=0)
%
% Generic over the event: one row of Ev per event, whatever detected it (the tracker's path, an atlas
% cell event). Its envelope is the AREA-WEIGHTED MEAN of the tile envelopes over the event's space tiles,
% at the envelope's native rate (not quantised to time tiles); rheome.detect.risefall is fitted to it over
% the event's first-to-last time ± Margin. tOn / tOff are the half-rise / half-fall times (MS1 SPEC1 §8 A).
%
% INPUTS
%   Ev   table, one row per event: hemi (string), tiles (cell of tile indices), t0, t1 (first and last
%        time of the event, s, on Y's clock)
%   Y    struct, Y.(hemi) [nTile x nT] tile envelope (e.g. tile means of rheome.flow.envelopemodes)
%   A    struct, A.(hemi) [nTile x 1] tile area (m^2)
%   fs   Y's sample rate          T0  time of Y's first sample (s)
%   R2Min  a fit counts if r2 >= R2Min (0.7, as rheome.flow.labelalpha)
%   Model  risefall's model: "linear" or "quadrature" (the episode adds to the background in power)
%
% OUTPUT T = Ev with, per event:
%   tOn tOff tauR tauF r2 (NaN without a fit)   fit  (fit made and r2 >= R2Min)   atBound (tau at TauMax)
%   durS = tOff - tOn   buildupS = 4 tauR   decayS = 4 tauF   plateauS (0 when the fit leaves none)
%   height (h, envelope units)   peakAmp (max of the envelope in [tOn, tOff])   nTiles   areaMM2
% Without an accepted fit the shape columns stay as fitted (NaN when none) and `fit` is false: the
% unfitted share is reported, never hidden.
%
% See also: rheome.detect.risefall, rheome.detect.eventlocked, rheome.scale.measure_events
%
% Author: Diellor Basha, 2026

    arguments
        Ev table
        Y struct
        A struct
        fs (1,1) double {mustBePositive}
        opts.Margin (1,1) double = 1
        opts.Smooth (1,1) double = 0.1
        opts.TauMax (1,1) double = 1
        opts.R2Min (1,1) double = 0.7
        opts.T0 (1,1) double = 0
        opts.Model (1,1) string = "linear"
    end
    n = height(Ev);  v = nan(n, 15);
    for i = 1:n
        h = char(Ev.hemi(i));  k = unique(Ev.tiles{i}(:));  w = A.(h)(k) / sum(A.(h)(k));
        a = w' * double(Y.(h)(k, :));
        R = rheome.detect.risefall(a, fs, Ev.t0(i) - opts.T0, Ev.t1(i) - opts.T0, ...
                                   Margin=opts.Margin, Smooth=opts.Smooth, TauMax=opts.TauMax, Model=opts.Model);
        pk = NaN;  pl = NaN;
        if R.ok
            s = max(1, round(R.tOn*fs) + 1):min(numel(a), round(R.tOff*fs) + 1);
            if ~isempty(s), pk = max(a(s)); end
            pl = 0;  if ~isempty(R.plateau), pl = diff(R.plateau); end
        end
        v(i, :) = [R.tOn + opts.T0, R.tOff + opts.T0, R.tauR, R.tauF, R.r2, R.ok && R.r2 >= opts.R2Min, ...
                   R.atBound, R.tOff - R.tOn, 4*R.tauR, 4*R.tauF, pl, R.h, pk, numel(k), 1e6 * sum(A.(h)(k))];
    end
    X = array2table(v, 'VariableNames', {'tOn','tOff','tauR','tauF','r2','fit','atBound','durS', ...
        'buildupS','decayS','plateauS','height','peakAmp','nTiles','areaMM2'});
    X.fit = X.fit == 1;  X.atBound = X.atBound == 1;
    T = [Ev, X];
end

% Author: Diellor Basha, 2026
