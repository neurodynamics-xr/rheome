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
%            and tau can stay compressed: the beat, not the fit, sets that limit. ~0.3 s per event.
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
    dec = max(1, round(numel(t)/max(t(end) - t(1), eps)/RATE));
    td = t(1:dec:end)';  Pd = (y(1:dec:end).^2)';  n = numel(td);
    w = 1 ./ max(movmean(Pd, 3), 1e-3);                                        % amplitude noise ~ constant: var(P) ~ P
    onG = t(1):STEP:t1;  offG = t0:STEP:t(end);  tauG = exp(linspace(log(TAUMIN), log(opts.TauMax), NTAU));
    if isempty(onG) || isempty(offG), return; end
    [I, J] = ndgrid(1:numel(onG), 1:numel(offG));  k = offG(J) > onG(I);  I = I(k);  J = J(k);
    if isempty(I), return; end
    sw = sum(w);  swp = w'*Pd;  swpp = w'*Pd.^2;  L = zeros(numel(I), NTAU, NTAU);
    for a1 = 1:NTAU
        SR = S((td - onG)/tauG(a1));
        for a2 = 1:NTAU
            G = (SR(:, I) .* S(-(td - offG(J))/tauG(a2))).^2;
            swg = w'*G;  swgg = w'*G.^2;  swgp = (w.*Pd)'*G;
            H = (sw*swgp - swg*swp) ./ max(sw*swgg - swg.^2, realmin);  B = (swp - H.*swg)/sw;
            e = swpp - 2*B*swp - 2*H.*swgp + B.^2*sw + 2*B.*H.*swg + H.^2.*swgg;  % weighted SSE at the profiled B, H
            z = B < 0;                                                           % B clipped at 0
            if any(z), H(z) = swgp(z)./max(swgg(z), realmin);  e(z) = swpp - H(z).^2.*swgg(z); end
            e(H <= 0) = Inf;  L(:, a1, a2) = e;
        end
    end
    [emin, kb] = min(L(:));  if ~isfinite(emin), return; end
    [ib, a1, a2] = ind2sub(size(L), kb);                                        % the mode: for the residuals' ACF
    X = [ones(n, 1), (S((td - onG(I(ib)))/tauG(a1)) .* S(-(td - offG(J(ib)))/tauG(a2))).^2];
    r = Pd - X*((X.*sqrt(w)) \ (Pd.*sqrt(w)));  r = r - mean(r);
    ac = xcorr(r, 'coeff');  ac = ac(n:end);  m = find(ac < 0, 1);  if isempty(m), m = n; end
    Tint = max(1, 1 + 2*sum(ac(2:m-1)));                                         % frames per independent sample
    p = exp(-(L - emin) / (2*emin/n*Tint));  p = p / sum(p(:));
    pt = sum(sum(p, 3), 2);
    q = [i_median(onG, accumarray(I, pt, [numel(onG) 1])), i_median(offG, accumarray(J, pt, [numel(offG) 1]))];
    if q(2) <= q(1), return; end
    tq = exp([i_median(log(tauG), squeeze(sum(sum(p, 3), 1))), i_median(log(tauG), squeeze(sum(sum(p, 2), 1)))]);
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

% The median of a discrete distribution p over x, interpolated on its cumulative.
function m = i_median(x, p)
    c = cumsum(p(:)) / sum(p);  [c, u] = unique(c);  x = x(:);  x = x(u);
    if numel(c) < 2, m = x(1); else, m = interp1(c, x, 0.5, 'linear', 'extrap'); end
end

% Author: Diellor Basha, 2026
