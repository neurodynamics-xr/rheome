function [Ys, info] = laggedrest(Y, Gp, front, fs, opts)
% FORWARD.LAGGEDREST  The lagged-source null of a resting recording: rheome.forward.laggedpair window by window.
%
%   [Ys, info] = rheome.forward.laggedrest(Y, Gp, front, fs)
%   [Ys, info] = rheome.forward.laggedrest(Y, Gp, front, fs, WindowS=4, MaxLagS=0.5, Group=hemi, Mask=bad)
%
% Why: a resting direction reading (rheome.flow.upstreamshare) is a claim about travel only if it beats what
% generators that do NOT travel, at the participant's own lags and amplitude ratios, produce through the same
% head and the same readout (gate N of the AGE-1 and WAVE-OR-LEAKAGE test plans). An event has one window
% (rheome.forward.laggedpair); a resting span has no events, so it is cut into windows and each is replaced
% by its own lagged-pair null: generators matched to that stretch of the participant's sources, with no travel.
%
% 1. WINDOWS. Overlapping windows of WindowS seconds, hop WindowS - 2 MaxLagS; each window's central part (its
%    MaxLagS margins dropped, because laggedpair's delay is circular over the window) is kept, and the kept
%    parts tile the span edge to edge. The first and last windows keep their outer margin too.
% 2. NULL. In each window, rheome.forward.laggedpair: the frontal and the posterior generator (one candidate
%    patch per Group each, the best by least squares), the back one replaced by the front one's waveform at
%    the window's ratio, polarity and lag, plus the residual with each channel's phases randomised.
% 3. MASK. A window whose samples touch Mask (artefact) is not fitted: its kept part is left at zero and
%    info.read is false there; the caller must not read those samples.
%
% ⚠ THE SPAN IS THE SIGNAL THE READOUT TAKES, already band-limited (AGE-1: the periodic 1-8 Hz part, real, at
%   the readout's frame rate). The surrogate is then made analytic and read exactly as the data are.
% ⚠ The window and lag limit are a declared choice, not a fit: MaxLagS = 0.5 s is half a cycle at 1 Hz, the
%   largest lag the readout can tell from its mirror (rheome.flow.upstreamshare, aliasing).
% ⚠ The random phases use the global stream: seed it (rng) before calling for a reproducible null.
% ⚠ NOT VALIDATED (2026-10-11, rheome.scale.measure_upstreamnull in two MEG heads):
%   lagged 20-40 mm discs that do not travel leave |U*| up to 0.24, as large as planted travel at 100-1000 mm,
%   0.5-2 m/s; travel power 0. Tile candidates fit compact generators at r2 ~0.7, and the per-channel phase-
%   randomised residual loses the background's zero-lag structure. Do not read U* as travel until redesigned.
%
% INPUTS
%   Y      [nChan x nT] real sensor signal    Gp [nChan x nPatch] candidate generator topographies
%   front  [1 x nPatch] logical, frontal candidates     fs  sampling rate (Hz)
% NAME-VALUE
%   WindowS (4)   MaxLagS (0.5)   Group [1 x nPatch] (one group)   Mask [1 x nT] logical (none)
% OUTPUT
%   Ys     [nChan x nT] the surrogate
%   info   .read [1 x nT] logical (samples carrying a surrogate)   .windows table: first, last, lagS, ratio,
%          r2, polarity, fitted (one row per window)
%
% See also: rheome.forward.laggedpair, rheome.scale.measure_upstream, rheome.flow.upstreamshare
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        Gp double
        front (1,:) logical
        fs (1,1) double {mustBePositive}
        opts.WindowS (1,1) double {mustBePositive} = 4
        opts.MaxLagS (1,1) double {mustBeNonnegative} = 0.5
        opts.Group (1,:) = ones(1, size(Gp, 2))
        opts.Mask (1,:) logical = false(1, size(Y, 2))
    end
    nT = size(Y, 2);  w = round(opts.WindowS * fs);  m = round(opts.MaxLagS * fs);  hop = w - 2*m;
    if hop < 1, error('rheome:laggedrest:window', 'WindowS (%g s) must exceed 2 MaxLagS (%g s)', opts.WindowS, 2*opts.MaxLagS); end
    if nT < w, error('rheome:laggedrest:short', 'span of %d samples is shorter than one %d-sample window', nT, w); end
    st = unique([1:hop:nT-w+1, nT-w+1]);  nW = numel(st);
    Ys = zeros(size(Y));  read = false(1, nT);
    W = table(zeros(nW,1), zeros(nW,1), nan(nW,1), nan(nW,1), nan(nW,1), nan(nW,1), false(nW,1), ...
              'VariableNames', {'first','last','lagS','ratio','r2','polarity','fitted'});
    kept = 0;                                        % last sample already written
    for i = 1:nW
        s = st(i):st(i)+w-1;
        k1 = max(st(i) + m, kept + 1);  if i == 1, k1 = 1; end
        k2 = st(i) + w - 1 - m;  if i == nW, k2 = nT; end
        W.first(i) = k1;  W.last(i) = k2;  kept = k2;
        if any(opts.Mask(s)), continue; end
        [Yw, fit] = rheome.forward.laggedpair(Y(:, s), Gp, front, fs, MaxLagS=opts.MaxLagS, Group=opts.Group);
        Ys(:, k1:k2) = Yw(:, k1 - st(i) + 1 : k2 - st(i) + 1);  read(k1:k2) = true;
        W{i, {'lagS','ratio','r2','polarity'}} = [fit.lagS fit.ratio fit.r2 fit.polarity];  W.fitted(i) = true;
    end
    info = struct('read', read, 'windows', W);
end

% Author: Diellor Basha, 2026
