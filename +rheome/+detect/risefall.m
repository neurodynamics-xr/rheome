function R = risefall(a, fs, t0, t1, opts)
% DETECT.RISEFALL  Where an amplitude episode switches on and off: a double-sigmoid fit.
%
%   R = rheome.detect.risefall(a, fs, t0, t1)                  % a candidate episode [t0, t1] s
%   R = rheome.detect.risefall(a, fs, t0, t1, Margin=1, Smooth=0.1)
%
%       a(t) = b + h * S((t - tOn)/tauR) * S(-(t - tOff)/tauF),    S(x) = 1/(1 + exp(-x))      (Model "linear")
%       a(t) = sqrt(b^2 + (h * S((t - tOn)/tauR) * S(-(t - tOff)/tauF))^2)                    (Model "quadrature")
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
%   Model    "linear" (default) or "quadrature". ⭐ An episode on an ongoing background adds to it in POWER,
%            not in amplitude (incoherent sources): |p + n|^2 ~ p^2 + n^2. The linear model then reads the
%            flanks of sqrt(b^2 + h^2 S^2) as sharper than the episode's own (planted through MEG at 10 dB:
%            tau x 0.65). "quadrature" fits the episode's own amplitude S under that sum.
%            "posterior": the quadrature model, read by its POSTERIOR rather than its least-squares point.
%            The double sigmoid in power, P = B + H g, B and H solved in closed form (weights 1/P), is profiled
%            over a grid of (tOn, tOff, tauR, tauF) (flat prior on risefall's box, log-spaced tau); the
%            likelihood is tempered by the residuals' integrated autocorrelation (an envelope's samples are
%            not independent); tOn and tOff are the marginal posterior MEDIANS (the estimate that minimises the
%            absolute timing error), then tauR, tauF are refitted with tOn, tOff held there. ⭐ Why: planted
%            through MEG at 10 dB in own rest, the timing error is set by the beat between the episode and the
%            ongoing alpha of the same band (|p + n|^2 = p^2 + n^2 + 2 p n cos(dphi), a slow multiplicative
%            noise), not by the fit's start, bounds, margin or smoothing; the least-squares point then trades
%            timing error into tau (x0.7-0.8). Held-out plants, 4 participants at 10 dB, vs "quadrature": median
%            |tOn| 0.13-0.14 s (0.17-0.19), |tOff| 0.11-0.18 s (0.12-0.22), tau x0.73-0.98 (x0.62-0.87), fit
%            share 0.90-0.99 (0.81-0.88); at 20 dB tau x0.97-1.02. ⚠ At 10 dB timing sits near one 165 ms tile
%            and tau can stay compressed: the beat, not the fit, sets that limit. ~0.6 s per event; bounded past a 10-s window (GridMax, FramesMax: ~10 s at 120 s).
%   GridMax  "posterior" only: at most this many tOn (and tOff) points per pass (400: the 25-ms grid up to a 10-s window);
%            a longer window is read coarse, then at 25 ms over its posterior's support
%   FramesMax  "posterior" only: at most this many frames (1200: 30 Hz up to a 40-s window, slower past it)
%   TauMax   upper bound of tauR and tauF (1 s). On two OMEGA participants tau reached 1 s at the 95th
%            percentile (labelling design section 7); a fit AT the bound is flagged (.atBound), not trusted
%
% OUTPUT (struct R)
%   .tOn .tOff .tauR .tauF .b .h (seconds and amplitude units)   .r2   .ok
%   .atBound  tauR or tauF within 1% of TauMax
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
        opts.TauMax (1,1) double {mustBePositive} = 1
        opts.Model (1,1) string {mustBeMember(opts.Model, ["linear" "quadrature" "posterior"])} = "linear"
        opts.GridMax (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(opts.GridMax, 2)} = 400
        opts.FramesMax (1,1) double {mustBeInteger, mustBePositive} = 1200
    end
    R = struct('tOn',NaN,'tOff',NaN,'tauR',NaN,'tauF',NaN,'b',NaN,'h',NaN,'r2',NaN,'ok',false,'atBound',false, ...
        'buildup',[],'plateau',[],'decay',[]);
    i0 = max(1, round((t0 - opts.Margin) * fs) + 1);  i1 = min(numel(a), round((t1 + opts.Margin) * fs) + 1);
    if i1 - i0 < 10, return; end
    y = movmean(a(i0:i1), max(1, round(opts.Smooth * fs)));  t = ((i0:i1) - 1) / fs;
    % ⚠ FIT ON A UNIT SCALE. A source envelope is ~1e-11 A.m; the optimiser's finite differences and
    % tolerances do not reach that scale and every fit stayed at its starting point (tau = 0.100 exactly).
    sc = max(abs(y));  if ~(sc > 0), return; end
    y = y / sc;
    S = @(x) 1 ./ (1 + exp(-x));
    if opts.Model == "posterior"
        R = i_posterior(R, y, t, t0, t1, sc, S, opts);  return
    end
    f = @(p, t) p(1) + p(2) * S((t - p(3)) / p(4)) .* S(-(t - p(5)) / p(6));
    b0 = prctile(y, 5);  h0 = max(y) - b0;
    p0 = [b0 h0 t0 0.1 t1 0.1];
    lb = [min(y) - h0, 0,       t(1), 0.01, t0,   0.01];
    if opts.Model == "quadrature"
        f = @(p, t) sqrt(p(1)^2 + (p(2) * S((t - p(3)) / p(4)) .* S(-(t - p(5)) / p(6))).^2);
        h0 = sqrt(max(max(y)^2 - b0^2, 0));  p0(2) = h0;  lb(1) = 0;
    end
    ub = [max(y),      3 * h0,  t1,   opts.TauMax,  t(end), opts.TauMax];
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
    R.atBound = max(R.tauR, R.tauF) >= 0.99 * opts.TauMax;
    R.ok = true;
end


% Model "posterior": profile posterior of the power double sigmoid, marginal medians of tOn, tOff, then tau given them.
function R = i_posterior(R, y, t, t0, t1, sc, S, opts)
    STEP = 0.025;  NTAU = 16;  RATE = 30;  TAUMIN = 0.02;                       % grid: 25 ms, 16 log-spaced tau
    % ⚠ TIME: cost is frames x tOn points x tOff points x 256, i.e. the window CUBED at a fixed step (cf-events: >30 min
    % on one event). Declared bound: <= FramesMax frames (the rate drops below 30 Hz past FramesMax/30 s) and <= GridMax
    % points per axis. A capped axis is read coarse first, then at 25 ms (or GridMax points) over its posterior's
    % support (mass outside < 1e-6). Below the caps the result is the uncapped grid's, to rounding.
    span = max(t(end) - t(1), eps);
    dec = max([1, round(numel(t)/span/RATE), ceil(numel(t)/opts.FramesMax)]);
    td = t(1:dec:end)';  Pd = (y(1:dec:end).^2)';
    w = 1 ./ max(movmean(Pd, 3), 1e-3);                                        % amplitude noise ~ constant: var(P) ~ P
    tauG = exp(linspace(log(TAUMIN), log(opts.TauMax), NTAU));
    onG = i_axis(t(1), t1, STEP, opts.GridMax);  offG = i_axis(t0, t(end), STEP, opts.GridMax);
    if isempty(onG) || isempty(offG) || ~any(offG(:) > onG, 'all'), return; end
    [p, emin] = i_post(i_sse(td, Pd, w, onG, offG, tauG, S), td, Pd, w, onG, offG, tauG, S);
    if ~isfinite(emin), return; end
    if numel(onG) > 1 && onG(2) - onG(1) > STEP*(1 + 1e-9) || numel(offG) > 1 && offG(2) - offG(1) > STEP*(1 + 1e-9)
        onG = i_refine(onG, sum(p, [2 3 4]), STEP, opts.GridMax);  offG = i_refine(offG, sum(p, [1 3 4]), STEP, opts.GridMax);
        if ~any(offG(:) > onG, 'all'), return; end
        [p, emin] = i_post(i_sse(td, Pd, w, onG, offG, tauG, S), td, Pd, w, onG, offG, tauG, S);
        if ~isfinite(emin), return; end
    end
    q = [i_median(onG, sum(p, [2 3 4])), i_median(offG, sum(p, [1 3 4]))];
    if q(2) <= q(1), return; end
    tq = exp([i_median(log(tauG), sum(p, [1 2 4])), i_median(log(tauG), sum(p, [1 2 3]))]);
    f = @(x, t) sqrt(x(1)^2 + (x(2)*S((t - q(1))/x(3)) .* S(-(t - q(2))/x(4))).^2);
    b0 = prctile(y, 5);  x0 = [b0, sqrt(max(max(y)^2 - b0^2, 1e-6)), tq];
    lb = [0 0 0.01 0.01];  ub = [max(y), 3*x0(2), opts.TauMax, opts.TauMax];  x0 = min(max(x0, lb), ub);
    try
        x = lsqcurvefit(f, x0, t, y, lb, ub, optimoptions('lsqcurvefit', 'Display', 'off'));
    catch
        return
    end
    r = y - f(x, t);
    R.b = x(1)*sc;  R.h = x(2)*sc;  R.tOn = q(1);  R.tauR = x(3);  R.tOff = q(2);  R.tauF = x(4);
    R.r2 = 1 - sum(r.^2) / max(sum((y - mean(y)).^2), realmin);
    R.buildup = [R.tOn - 2*R.tauR, R.tOn + 2*R.tauR];  R.decay = [R.tOff - 2*R.tauF, R.tOff + 2*R.tauF];
    if R.tOff - 2*R.tauF > R.tOn + 2*R.tauR, R.plateau = [R.tOn + 2*R.tauR, R.tOff - 2*R.tauF]; end
    R.atBound = max(R.tauR, R.tauF) >= 0.99 * opts.TauMax;  R.ok = true;
end

% Grid a:step:b, or nMax points over [a, b] when that is more.
function g = i_axis(a, b, step, nMax)
    g = a:step:b;
    if numel(g) > nMax, g = linspace(a, b, nMax); end
end

% The fine grid over a coarse axis' posterior support (cumulative 1e-6 .. 1 - 1e-6, one coarse point either side).
function g = i_refine(g, pm, step, nMax)
    c = cumsum(pm(:)) / sum(pm);  lo = max(1, find(c >= 1e-6, 1) - 1);  hi = min(numel(g), find(c >= 1 - 1e-6, 1) + 1);
    g = i_axis(g(lo), g(hi), step, nMax);
end

% Weighted SSE of P = B + H g at the profiled B, H for every (tOn, tOff, tauR, tauF): [nOn x nOff x NTAU x NTAU], Inf where
% tOff <= tOn or H <= 0. Its sums over time are matrix products: g = SR^2 SF^2, so w'g = SR2' (w.SF2) for all pairs at once.
function L = i_sse(td, Pd, w, onG, offG, tauG, S)
    n = numel(td);  nOn = numel(onG);  nOff = numel(offG);  NTAU = numel(tauG);
    sw = sum(w);  swp = w'*Pd;  swpp = w'*Pd.^2;
    SF2 = zeros(n, nOff, NTAU);
    for a2 = 1:NTAU, SF2(:, :, a2) = S(-(td - offG)/tauG(a2)).^2; end
    SF2 = reshape(SF2, n, []);  wSF2 = w.*SF2;  wpSF2 = (w.*Pd).*SF2;  wSF4 = wSF2.*SF2;  clear SF2
    bad = repmat(offG <= onG(:), 1, NTAU);  L = zeros(nOn, nOff, NTAU, NTAU);
    for a1 = 1:NTAU
        SR2 = S((td - onG)/tauG(a1)).^2;
        swg = SR2'*wSF2;  swgg = (SR2.^2)'*wSF4;  swgp = SR2'*wpSF2;            % [nOn x nOff*NTAU]
        H = (sw*swgp - swg*swp) ./ max(sw*swgg - swg.^2, realmin);  B = (swp - H.*swg)/sw;
        e = swpp - 2*B*swp - 2*H.*swgp + B.^2*sw + 2*B.*H.*swg + H.^2.*swgg;      % weighted SSE at the profiled B, H
        z = B < 0;                                                               % B clipped at 0
        if any(z, 'all'), H(z) = swgp(z)./max(swgg(z), realmin);  e(z) = swpp - H(z).^2.*swgg(z); end
        e(H <= 0 | bad) = Inf;  L(:, :, a1, :) = reshape(e, nOn, nOff, 1, NTAU);
    end
end

% The tempered posterior over L (in place), its minimum and mode: the temperature from the residuals' integrated ACF.
function [p, emin, kb] = i_post(p, td, Pd, w, onG, offG, tauG, S)
    n = numel(td);  [emin, kb] = min(p(:));  if ~isfinite(emin), return; end
    [i, j, a1, a2] = ind2sub(size(p), kb);
    X = [ones(n, 1), (S((td - onG(i))/tauG(a1)) .* S(-(td - offG(j))/tauG(a2))).^2];
    r = Pd - X*((X.*sqrt(w)) \ (Pd.*sqrt(w)));  r = r - mean(r);
    ac = xcorr(r, 'coeff');  ac = ac(n:end);  m = find(ac < 0, 1);  if isempty(m), m = n; end
    Tint = max(1, 1 + 2*sum(ac(2:m-1)));                                         % frames per independent sample
    T = 2*emin/n*Tint;  p = exp(-(p - emin) / T);  p = p / sum(p(:));
end

% The median of a discrete distribution p over x, interpolated on its cumulative.
function m = i_median(x, p)
    c = cumsum(p(:)) / sum(p);  [c, u] = unique(c);  x = x(:);  x = x(u);
    if numel(c) < 2, m = x(1); else, m = interp1(c, x, 0.5, 'linear', 'extrap'); end
end

% Author: Diellor Basha, 2026
