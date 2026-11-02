function R = burstfit(a, fps, t0, t1, opts)
% CORES.BURSTFIT  Where an amplitude episode switches on and off: the quadrature double-sigmoid fit, nxr's version.
%
%   R = rheome.cores.burstfit(a, fps, t0, t1, Margin=, Smooth=, TauMax=)
%
%   a(t) = sqrt(b^2 + (h S((t - tOn)/tauR) S(-(t - tOff)/tauF))^2),   S(x) = 1/(1 + e^-x), S clipped at |x| <= 50
%
% fitted to the amplitude a [1 x m] (sampled at fps) over [t0 - Margin, t1 + Margin] seconds, after a moving mean of Smooth
% seconds (rheome.cores.boxsmooth: zero-padded) and a division by max|y| (b, h are multiplied back; r2 is scale-free). Index
% window (0-based) i0 = max(0, round((t0-Margin)*fps)), i1 = min(m-1, round((t1+Margin)*fps)) inclusive, round half to even; fewer
% than 10 samples between them: not ok. Smoothing length max(1, round(Smooth*fps)), half to even.
% Start p0 = [b0, h0, t0, 0.1, t1, 0.1], b0 = 5th percentile of y by numpy's linear rule (rheome.cores.linpercentile),
% h0 = sqrt(max(max(y)^2 - b0^2, 0)); bounds lb = [0 0 t(1) 0.01 t0 0.01], ub = [max(y) 3 max(h0,1e-6) t1 TauMax t(end) TauMax],
% p0 clipped to [lb, ub - 1e-12]. Optimiser: lsqcurvefit (trust-region-reflective; nxr uses scipy least_squares, trf) with
% tolerances tightened to 1e-10: both are local least squares from the same start, so parameters agree to ~1e-3 s, not to the last
% digit (see tests/tVortexSignal for the achieved agreement). Fails (ok = false) when tOff <= tOn. r2 = 1 - SSE/max(SST, 1e-300).
% atBound = max(tauR, tauF) >= 0.99 TauMax.
%
% OUTPUT R (struct): .tOn .tOff .tauR .tauF (s) .b .h .r2 .ok .atBound; .buildup [tOn +- 2 tauR], .decay [tOff +- 2 tauF], .plateau
% [tOn + 2 tauR, tOff - 2 tauF] ([] when none), each [start end] s. The same fields as rheome.detect.risefall.
%
% Differs from rheome.detect.risefall (which is unchanged, for the sleep-spindle detectors): only the quadrature model; start
% clip ub - 1e-12, the 1e-6 floor in ub, the exp clip and the 1e-300 floor of r2 are nxr's; margin, smoothing and TauMax in seconds
% are passed by rheome.cores.burststages as 10, 1 and 10 cycles. Ports nxr-analysis src/nxr/vortex.py risefall (line 801, commit 30b9622).
%
% See also: rheome.cores.burststages, rheome.cores.burstenvelope, rheome.detect.risefall
%
% Author: Diellor Basha, 2026
    arguments
        a (1,:) double
        fps (1,1) double {mustBePositive}
        t0 (1,1) double
        t1 (1,1) double
        opts.Margin (1,1) double {mustBeNonnegative} = 1
        opts.Smooth (1,1) double {mustBeNonnegative} = 0.1
        opts.TauMax (1,1) double {mustBePositive} = 1
    end
    R = struct('tOn',NaN,'tOff',NaN,'tauR',NaN,'tauF',NaN,'b',NaN,'h',NaN,'r2',NaN,'ok',false,'atBound',false, ...
        'buildup',[],'plateau',[],'decay',[]);
    i0 = max(0, evenround((t0 - opts.Margin) * fps));
    i1 = min(numel(a) - 1, evenround((t1 + opts.Margin) * fps));
    if i1 - i0 < 10, return; end
    k = max(1, evenround(opts.Smooth * fps));
    y = rheome.cores.boxsmooth(a(i0+1:i1+1), k);
    t = (i0:i1) / fps;
    sc = max(abs(y));  if ~(sc > 0), return; end
    y = y / sc;
    S = @(x) 1 ./ (1 + exp(-min(max(x, -50), 50)));
    f = @(p, tt) sqrt(p(1)^2 + (p(2) * S((tt - p(3)) / p(4)) .* S(-(tt - p(5)) / p(6))).^2);
    b0 = rheome.cores.linpercentile(y, 5);  h0 = sqrt(max(max(y)^2 - b0^2, 0));
    lb = [0, 0, t(1), 0.01, t0, 0.01];
    ub = [max(y), 3 * max(h0, 1e-6), t1, opts.TauMax, t(end), opts.TauMax];
    p0 = min(max([b0, h0, t0, 0.1, t1, 0.1], lb), ub - 1e-12);
    o = optimoptions('lsqcurvefit', 'Display', 'off', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-10, 'MaxIterations', 1000, 'MaxFunctionEvaluations', 20000);
    try
        p = lsqcurvefit(f, p0, t, y, lb, ub, o);
    catch
        return
    end
    if p(5) <= p(3), return; end
    r = y - f(p, t);
    R.b = p(1) * sc;  R.h = p(2) * sc;  R.tOn = p(3);  R.tauR = p(4);  R.tOff = p(5);  R.tauF = p(6);
    R.r2 = 1 - sum(r.^2) / max(sum((y - mean(y)).^2), 1e-300);
    R.buildup = [R.tOn - 2*R.tauR, R.tOn + 2*R.tauR];
    R.decay = [R.tOff - 2*R.tauF, R.tOff + 2*R.tauF];
    if R.tOff - 2*R.tauF > R.tOn + 2*R.tauR
        R.plateau = [R.tOn + 2*R.tauR, R.tOff - 2*R.tauF];
    end
    R.atBound = max(R.tauR, R.tauF) >= 0.99 * opts.TauMax;
    R.ok = true;
end
