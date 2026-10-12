function [Ys, fit] = laggedpair(Y, Gp, front, fs, opts)
% FORWARD.LAGGEDPAIR  The lagged-source null of one event: two fixed generators fitted to it, with its amplitude ratio and lag and no travel, back through the same head with its own noise.
%
%   [Ys, fit] = rheome.forward.laggedpair(Y, Gp, front, fs)
%   [Ys, fit] = rheome.forward.laggedpair(Y, Gp, front, fs, MaxLagS=0.3, Group=hemi)   % bilateral generators
%
% Why: a frontal and a posterior generator that only differ in lag make a field that every reader of direction can
% call a wave (MS1 G16b: the framework declares propagation in about half of two-lagged-generator nulls). A
% direction claim is a claim about travel only if it beats this null, matched to the participant's own event.
%
% 1. FIT. Every candidate is a patch with a fixed sensor topography (a column of Gp: e.g. a tile of normal current
%    through the participant's own gain). A generator is one patch per Group (e.g. per hemisphere: slow
%    oscillations are bilateral, and one patch for both hemispheres misfits them -- synthetic, lag 0.15 s planted
%    in both hemispheres: one pair reads 0.02 s, one patch per hemisphere 0.15 s). Over every combination of a
%    frontal patch (front true) and a posterior one (front false) in each group, the combination whose time courses
%    explain most of Y by least squares is kept; each generator's patches are then reduced to their first singular
%    pair (one topography, one time course: s1 front, s2 back), the rest joining the residual.
% 2. LAG AND AMPLITUDE. The lag tau is the peak of |xcorr(s2, s1)| within +-MaxLagS (tau > 0: the back lags the
%    front), the amplitude ratio rms(s2)/rms(s1), the polarity the sign of that peak.
% 3. NULL. The back generator is replaced by the front one's waveform, scaled by the ratio and delayed by tau, so
%    the two generators are fixed in place and share one waveform: nothing travels. Ys = the pair's topographies x
%    [s1; ratio * polarity * s1(t - tau)] + the residual Y - fit with each channel's Fourier phases randomised
%    (matched noise: each channel keeps its spectrum of everything the pair does not explain; any travel left in
%    the residual loses its cross-channel lags).
% OPTIONS (default = the null above):
%    Keep=true        the two generators keep their OWN fitted time courses s1 and s2 (still two fixed topographies:
%                     nothing travels), instead of the back one becoming the front one's lagged copy.
%    Residual="zerolag"  the residual's Fourier phases are randomised by ONE phase per frequency for every channel,
%                     applied to the real part of that frequency's spatial pattern: it keeps the residual's zero-lag
%                     cross-spectrum (field spread) and drops its lags. "perchan" (default) randomises each channel.
%    MaxCombos (1e4)  above this many candidate combinations the best one is found by coordinate ascent (each
%                     slot in turn takes its best candidate, the others fixed), not by trying every one.
% ⚠ THE DELAY IS CIRCULAR over the window. Give the window a margin of at least MaxLagS beyond what is read.
% ⚠ THE PAIR IS THE BEST OF THE CANDIDATES, not a free dipole fit: positions are as fine as the patches.
% ⚠ The random phases use the global stream: seed it (rng) before calling for a reproducible null.
%
% INPUTS
%   Y      [nChan x nT] the event window, in the sensor space and reference Gp is in
%   Gp     [nChan x nPatch] candidate generator topographies     front  [1 x nPatch] logical, frontal candidates
%   fs     sampling rate (Hz)            Group  [1 x nPatch] labels, one patch per label per generator (default: one)
% OUTPUT
%   Ys     [nChan x nT] the surrogate
%   fit    .front .back (patch indices, one per group, in unique(Group) order) .wFront .wBack (each patch's
%          share of its generator, energy^(1/2)) .ratio (back over front, source units) .lagS .polarity .r2
%          (fraction of Y the two generators explain)
%
% See also: rheome.scale.measure_slowosc, rheome.forward.leadfield
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        Gp double
        front (1,:) logical
        fs (1,1) double
        opts.MaxLagS (1,1) double = 0.3
        opts.Group (1,:) = ones(1, size(Gp, 2))
        opts.Keep (1,1) logical = false
        opts.Residual (1,1) string {mustBeMember(opts.Residual, ["perchan" "zerolag"])} = "perchan"
        opts.MaxCombos (1,1) double {mustBePositive} = 1e4
    end
    K = Gp' * Gp;  B = Gp' * Y;  C = B * B';                      % explained energy of a set c: trace(K_cc \ C_cc)
    gr = unique(opts.Group, 'stable');  nG = numel(gr);  cand = cell(1, 2 * nG);
    for g = 1:nG
        cand{g} = find(front & opts.Group == gr(g));  cand{nG + g} = find(~front & opts.Group == gr(g));
    end
    cand = cand(~cellfun(@isempty, cand));  nG = numel(cand) / 2;
    assert(nG == round(nG) && nG >= 1, 'rheome:laggedpair:candidates', 'every group needs a frontal and a posterior candidate');
    if prod(cellfun(@numel, cand)) <= opts.MaxCombos                    % every combination
        I = cell(1, 2 * nG);  [I{:}] = ndgrid(cand{:});  I = cell2mat(cellfun(@(x) x(:), I, 'UniformOutput', false));
        e = zeros(size(I, 1), 1);
        for r = 1:size(I, 1), c = I(r, :);  e(r) = i_E(K, C, c); end
        [~, r] = max(e);  c = I(r, :);
    else, c = i_ascent(K, C, cand);                                    % too many: coordinate ascent
    end
    fit.front = c(1:nG);  fit.back = c(nG+1:end);
    s = Gp(:, c) \ Y;
    [AF, s1, fit.wFront] = i_rank1(Gp(:, fit.front), s(1:nG, :));
    [AB, s2, fit.wBack] = i_rank1(Gp(:, fit.back), s(nG+1:end, :));
    R = Y - AF * s1 - AB * s2;
    fit.r2 = 1 - sum(R(:).^2) / sum(Y(:).^2);
    m = round(opts.MaxLagS * fs);  [xc, lg] = xcorr(s2, s1, m);
    [~, k] = max(abs(xc));  fit.lagS = lg(k) / fs;  fit.polarity = sign(xc(k));
    fit.ratio = rms(s2) / rms(s1);
    if opts.Residual == "zerolag", R = i_zerolag(R); else, R = i_phaserand(R); end
    if opts.Keep, Ys = AF * s1 + AB * s2 + R;
    else, Ys = AF * s1 + AB * (fit.ratio * fit.polarity * circshift(s1, lg(k))) + R; end
end

function e = i_E(K, C, c)
% the energy of Y that the candidate set c explains by least squares (a repeated candidate explains nothing new)
    if numel(unique(c)) < numel(c), e = -Inf; else, e = trace(K(c, c) \ C(c, c)); end
end

function c = i_ascent(K, C, cand)
% each slot (a group's front or back) starts at its best single candidate; then each slot in turn takes the
% candidate that, with the others fixed, explains most, until a sweep changes nothing (at most 10 sweeps)
    ns = numel(cand);  c = zeros(1, ns);
    for s = 1:ns, [~, j] = max(arrayfun(@(k) i_E(K, C, k), cand{s}));  c(s) = cand{s}(j); end
    for sweep = 1:10
        moved = false;
        for s = 1:ns
            e = arrayfun(@(k) i_E(K, C, [c(1:s-1) k c(s+1:end)]), cand{s});  [~, j] = max(e);
            if cand{s}(j) ~= c(s), c(s) = cand{s}(j);  moved = true; end
        end
        if ~moved, break; end
    end
end

function [a, s, w] = i_rank1(A, S)
% a generator of several patches with one time course: the first singular pair of their summed SENSOR contribution
% (not of their time courses: near-collinear patches get large cancelling time courses that the sum removes);
% w = each patch's share of that contribution's energy. s is in the patches' units (a unit-norm vector of patch
% currents), so the ratio of two generators is a ratio of source amplitudes
    [U, D, V] = svd(A * S, 'econ');  u = A \ U(:, 1);                 % the patch currents making that topography
    k = norm(u) * sign(sum(u) + (sum(u) == 0));  a = U(:, 1) / k;  s = k * D(1) * V(:, 1)';   % s in source units
    e = zeros(size(A, 2), 1);  for i = 1:size(A, 2), e(i) = norm(A(:, i) * S(i, :), 'fro')^2; end
    w = sqrt(e / sum(e));
end

function y = i_phaserand(y)
% each channel with its own random Fourier phases: its amplitude spectrum kept, its lags to other channels destroyed
    n = size(y, 2);  h = floor((n - 1) / 2);  Fy = fft(y, [], 2);
    ph = exp(2i * pi * rand(size(y, 1), h));
    Fy(:, 2:h+1) = Fy(:, 2:h+1) .* ph;  Fy(:, n:-1:n-h+1) = Fy(:, n:-1:n-h+1) .* conj(ph);
    y = real(ifft(Fy, [], 2));
end

function y = i_zerolag(y)
% one random phase per frequency, common to every channel, applied to the REAL part of that frequency's spatial
% pattern: the zero-lag cross-spectrum (field spread, in-phase sources) is kept in expectation, every lag between
% channels (its imaginary part) is dropped
    n = size(y, 2);  h = floor((n - 1) / 2);  Fy = fft(y, [], 2);
    X = sqrt(2) * real(Fy(:, 2:h+1) .* exp(2i * pi * rand(1, h))) .* exp(2i * pi * rand(1, h));
    Fy(:, 2:h+1) = X;  Fy(:, n:-1:n-h+1) = conj(X);  y = real(ifft(Fy, [], 2));
end

% Author: Diellor Basha, 2026
