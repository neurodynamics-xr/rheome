function [T, X] = measure_patterns(name, S, item, opts)
% SCALE.MEASURE_PATTERNS  The catalogue's detectors on this participant's resting MEG against its own nulls (MS1 G5, G13).
%
%   [T, X] = rheome.scale.measure_patterns(name, S, "patterns")       % G5, and the G13 singularity count
%   [T, X] = rheome.scale.measure_patterns(name, S, "patternnulls")   % G13: speed sweep, dispersion, two-flow r
%
% Port of the single-recording catalogue analysis (MS1 Fig. 15C: one recording, 120 s, left hemisphere) to every
% participant, both hemispheres (Hemis), the whole imported record (rheome.scale.run DurationS, 300 s), the
% participant's own plain minimum norm (S.Res) read along vertex normals smoothed over NormalSigmaMM.
%
% patterns      Per frame (every FrameS of the analytic field, Band 8-12 Hz, and IAF +- IAFHalfWidth Hz):
%                 planarity     phase-gradient directionality |sum a k^| / sum a over faces
%                 source sink saddle vortex   rheome.detect.criticalPoints counts of the amplitude-weighted
%                               phase-velocity direction
%                 standing      |sum a z^2| / sum a |z|^2
%                 speed         median phase speed on faces above median amplitude (m/s)
%                 singularities phase singularities on faces inside the top 20% of the amplitude smoothed
%                               at 62 mm (the chirality analysis's mask: a singularity IS an amplitude zero,
%                               so masking on its own face removes it) -- G13's count
%               and per PacWindowS window, band "thetagamma":
%                 pac           median over vertices of the Theta-phase x Gamma-amplitude modulation index
%               NULLS. Primary: NSurr phase-randomised surrogates, one random phase per frequency bin shared
%               by every channel, so every auto- and cross-spectrum is kept. Secondary: the empty room through
%               the same kernel columns (channels by name); NaN when the participant has none.
%               ⭐ THE SURROGATES ARE EVALUATED, NOT REBUILT. Band-pass, analytic signal and phase
%               randomisation are all diagonal in frequency, so a surrogate's analytic field at any time is
%               (band bins .* random phases) * exp(2 pi i f t): only the SurrStrideS grid frames (and
%               PacSurrWindows windows) are computed, never a 300 s record per surrogate.
%               Per detector: thr = 95th percentile of all surrogate frames; f = fraction of the participant's
%               frames above thr (all frames); e = f - 0.05; p = (1 + #{f_s >= f_grid}) / (1 + NSurr), f_s and
%               f_grid on the SAME grid frames, so the p-value carries the autocorrelation the frame count
%               ignores; f_er against the empty room's 95th percentile; z_median, p_median: the participant's
%               grid median against the surrogates' medians (G13 for "singularities").
%               Halves 1 and 2 (frames in the first / second half of the record) for the ICC; half 0 = all.
%
% patternnulls  MS1 section 6.3, per participant and hemisphere, each against a null that can move it:
%                 peakedness    travelling-wave speed sweep (dynamics_spectral.m): max / mean over 0.05-1 m/s
%                               of the joint-energy share on rheome.filters.travwave's ridge (width 1.5 Hz), LBO
%                               coefficients of the normal current, NumWindows strongest WindowS alpha windows;
%                               median over windows. Single-subject reference: 1.32 (planted travelling 3.98, standing 2.07).
%                 dispersion    top / bottom |lambda| quintile ratio of each Dirac mode's periodic peak
%                               frequency (Welch 8 s, aperiodic knee fit 1-45 Hz, peak in 6-16 Hz), the
%                               hemisphere's Dirac modes applied to the MNE current. Single-subject reference: 1.01, a wave needs 3.1.
%                 rotation_r    corr of |curl v| (rheome.flow.apparent of the 8-16 Hz envelope, 300 Hz) with
%                               |curl J| over vertices, per 2 s tile, NumTiles tiles; median. Single-subject reference: r = -0.11.
%               ⚠⚠ THE PHASE-RANDOMISED SURROGATE CANNOT BE THE NULL FOR THE FIRST TWO. Both are functions of
%               each mode's power spectrum, and a shared random phase per bin leaves every mode's spectrum
%               EXACTLY unchanged (|Pk * Y_f * e^{i theta_f}|^2 = |Pk * Y_f|^2): z would be 0/0. Their null is
%               NPerm permutations of the mode labels (which lambda each spectrum sits at) -- H0: temporal
%               spectrum independent of spatial scale, i.e. no dispersion. rotation_r's null is the
%               tile-mismatched pairs (|curl v| of one tile against |curl J| of another).
%               z = (value - null mean) / null sd; p one-sided upper for peakedness and dispersion (travelling
%               and dispersive read high), two-sided for rotation_r. Halves as above (windows, tiles, Welch span).
%
% X is the per-unit table (rheome.scale.run writes it as <item>.csv); T the rheome.scale.rows.
%
% ⚠ The MNE resolution floor (median r50 ~52 mm on one participant) bounds every count: a mesoscale pattern cannot
%   appear here as itself (the catalogue validation's 'meso' rows; G6 measures that loss).
%
% See also: rheome.scale.run, rheome.flow.phasegradient, rheome.detect.criticalPoints,
%           rheome.detect.phasesingularity, rheome.filters.travwave, rheome.flow.apparent
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S
        item (1,1) string {mustBeMember(item, ["patterns" "patternnulls"])}
        opts.Hemis string = ["L" "R"]
        opts.Band (1,2) double = [8 12]
        opts.IAFHalfWidth (1,1) double = 2          % 0 = no IAF sensitivity band
        opts.Rate (1,1) double = 100                % frame clock (Hz): the second sample of each frame
        opts.FrameS (1,1) double = 0.2
        opts.NSurr (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.SurrStrideS (1,1) double = 4           % surrogate grid: 75 frames per 300 s surrogate
        opts.Theta (1,2) double = [4 8]
        opts.Gamma (1,2) double = [30 45]
        opts.PacWindowS (1,1) double = 2
        opts.PacSurrWindows (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.NormalSigmaMM (1,1) double = 20
        opts.NumWindows (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.WindowS (1,1) double = 3
        opts.NumTiles (1,1) double {mustBeInteger, mustBePositive} = 8
        opts.NPerm (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.Seed (1,1) double = 43
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    rng(opts.Seed);
    if item == "patterns", [T, X] = i_patterns(name, S, F, fs, opts);
    else, [T, X] = i_nulls(name, S, F, fs, opts);
    end
end

%% ---------- G5: detectors per frame against the surrogates and the empty room
function [T, X] = i_patterns(name, S, F, fs, o)
    Fy = fft(F, [], 2);  N = size(F, 2);  dur = N / fs;
    iaf = i_iaf(F, fs);
    bands = struct('name', "alpha", 'lim', o.Band);
    if o.IAFHalfWidth > 0 && isfinite(iaf), bands(2) = struct('name', "iaf", 'lim', iaf + [-1 1]*o.IAFHalfWidth); end
    E = [];  iE = [];  FyE = [];
    try, [E, iE] = i_emptyroom(name, S);  FyE = fft(E.F, [], 2);
    catch e, fprintf('[patterns %s] no empty-room arm: %s\n', name, e.message);
    end
    tA = 0 : o.FrameS : dur - 2/o.Rate;                    % every frame
    tG = 0 : o.SurrStrideS : dur - 2/o.Rate;               % the surrogate grid, a subset of tA
    [~, gi] = min(abs(tA(:) - tG(:)'), [], 1);
    X = table();
    for h = o.Hemis
        C = i_hemi(S, h, o);
        for b = bands
            B = i_band(Fy, N, fs, b.lim);
            xa = i_frames(C, B, tA, 1, o.Rate);                % [nA x nDet]
            xs = zeros(o.NSurr, numel(tG), size(xa, 2));
            for s = 1:o.NSurr
                rng(o.Seed + s);                               % the same surrogates in both hemispheres
                xs(s, :, :) = reshape(i_frames(C, B, tG, exp(2i*pi*rand(1, numel(B.f))), o.Rate), 1, numel(tG), []);
            end
            xe = NaN(1, size(xa, 2));
            if ~isempty(E)
                tE = 0 : o.FrameS : size(E.F, 2)/E.fs - 2/o.Rate;
                xe = i_frames(i_kcols(C, iE), i_band(FyE, size(E.F, 2), E.fs, b.lim), tE, 1, o.Rate);
            end
            X = [X; i_summarise(h, b.name, i_detnames(), xa, gi, xs, xe, tA, tG, dur)]; %#ok<AGROW>
            fprintf('[patterns %s] %s %s %.1f-%.1f Hz: %d frames, %d surrogates x %d\n', name, h, b.name, b.lim, numel(tA), o.NSurr, numel(tG));
        end
        X = [X; i_pac(C, Fy, N, fs, E, iE, FyE, dur, h, o)]; %#ok<AGROW>
    end
    X.iafHz = repmat(iaf, height(X), 1);
    k = X.half == 0;  lab = X.hemi(k) + "_" + X.band(k) + "_" + X.detector(k);
    T = [rheome.scale.rows("patterns", "iaf_hz", iaf, "Hz");
         rheome.scale.rows("patterns", repmat("excess_e", sum(k), 1), X.e(k), "fraction", lab);
         rheome.scale.rows("patterns", repmat("p_surrogate", sum(k), 1), X.p(k), "p", lab);
         rheome.scale.rows("patterns", repmat("frac_above_emptyroom", sum(k), 1), X.fER(k), "fraction", lab);
         rheome.scale.rows("patterns", repmat("z_median", sum(k), 1), X.zMedian(k), "z", lab)];
end

function n = i_detnames()
    n = ["planarity" "source" "sink" "saddle" "vortex" "standing" "speed" "singularities"];
end

function x = i_frames(C, B, t, r, rate)
% the detectors on frames at times t (s) of the band's analytic field, bins rotated by r
    Z = i_at(B, [t(:); t(:) + 1/rate]', r);               % [nCh x 2nt] -> each frame and its next sample
    nt = numel(t);  Zv = C.Kn * Z;  x = zeros(nt, numel(i_detnames()));
    for j = 1:nt, x(j, :) = i_frame(Zv(:, [j j+nt]), C, rate); end
end

function v = i_frame(z, C, rate)
    unitv = @(V) V ./ max(vecnorm(V, 2, 2), eps);
    Fc = double(C.Sh.Faces);
    pg = rheome.flow.phasegradient(z, C.Sh, 'Rate', rate);
    a = abs(z(:,1));  aF = mean(a(Fc), 2);  hi = aF > median(aF);
    kh = unitv(pg.k(:,:,1));  vel = pg.velocity(:,:,1);
    pgd = norm(sum(aF .* kh, 1)) / sum(aF);
    spd = median(pg.speed(hi, 1), 'omitnan');
    sidx = abs(sum(C.av .* z(:,1).^2)) / sum(C.av .* a.^2);
    U = (C.Afv * (aF .* unitv(vel))) ./ C.afvsum;  U(~isfinite(U)) = 0;
    ty = string(rheome.detect.criticalPoints(reshape(U.', [], 1), C.Sh, 'all', C.op).type(:));
    q = rheome.detect.phasesingularity(z(:,1), C.Sh).perFace(:);
    As = C.lbo.Phi * (C.heat62 .* (C.lbo.Phi' * (C.lbo.Mass * a)));  fc = mean(As(Fc), 2);
    v = [pgd nnz(ty == "source") nnz(ty == "sink") nnz(ty == "saddle") nnz(ty == "vortex") sidx spd ...
         nnz(q ~= 0 & fc >= prctile(fc, 80))];
end

function X = i_pac(C, Fy, N, fs, E, iE, FyE, dur, h, o)
% theta-phase x gamma-amplitude modulation index per window, median over vertices
    Bt = i_band(Fy, N, fs, o.Theta);  Bg = i_band(Fy, N, fs, o.Gamma);
    nw = round(o.PacWindowS * o.Rate);  w0 = 0 : o.PacWindowS : dur - o.PacWindowS;
    sub = unique(round(linspace(1, numel(w0), min(o.PacSurrWindows, numel(w0)))));
    mi = @(Ct, Bt_, Bg_, t0, rt, rg) median(i_mi(Ct.Kn * i_at(Bt_, t0 + (0:nw-1)/o.Rate, rt), ...
                                                  Ct.Kn * i_at(Bg_, t0 + (0:nw-1)/o.Rate, rg)));
    xa = arrayfun(@(t0) mi(C, Bt, Bg, t0, 1, 1), w0)';
    xs = zeros(o.NSurr, numel(sub));
    for s = 1:o.NSurr
        rng(o.Seed + s);  rt = exp(2i*pi*rand(1, numel(Bt.f)));  rg = exp(2i*pi*rand(1, numel(Bg.f)));
        xs(s, :) = arrayfun(@(t0) mi(C, Bt, Bg, t0, rt, rg), w0(sub));
    end
    xe = NaN;
    if ~isempty(E)
        Ce = i_kcols(C, iE);  nE = size(E.F, 2);
        Bte = i_band(FyE, nE, E.fs, o.Theta);  Bge = i_band(FyE, nE, E.fs, o.Gamma);
        xe = arrayfun(@(t0) mi(Ce, Bte, Bge, t0, 1, 1), 0 : o.PacWindowS : nE/E.fs - o.PacWindowS)';
    end
    X = i_summarise(h, "thetagamma", "pac", xa, sub, xs, xe, w0 + o.PacWindowS/2, w0(sub) + o.PacWindowS/2, dur);
end

function m = i_mi(Zt, Zg)
    ag = abs(Zg);  m = abs(mean(ag .* exp(1i*angle(Zt)), 2)) ./ mean(ag, 2);
end

function X = i_summarise(h, band, dets, xa, gi, xs, xe, tA, tG, dur)
% one row per detector x half (0 = all, 1, 2): f, e, p, f_er, z of the median
    X = table();
    for d = 1:numel(dets)
        for half = 0:2
            ka = true(size(tA(:)));  kg = true(size(tG(:)));
            if half, ka = (tA(:) < dur/2) == (half == 1);  kg = (tG(:) < dur/2) == (half == 1); end
            a = xa(ka, d);  g = xa(gi(kg), d);  s = xs(:, kg, d);
            thr = prctile(s(isfinite(s)), 95);
            f = i_frac(a, thr);  fg = i_frac(g, thr);
            fs_ = sum(s > thr, 2) ./ max(sum(isfinite(s), 2), 1);
            p = (1 + sum(fs_ >= fg)) / (1 + numel(fs_));
            ms = median(s, 2, 'omitnan');  mg = median(g, 'omitnan');
            z = (mg - mean(ms)) / std(ms);  pm = (1 + sum(ms >= mg)) / (1 + numel(ms));
            e = xe(:, d);  e = e(isfinite(e));  fER = NaN;  if ~isempty(e), fER = i_frac(a, prctile(e, 95)); end
            X = [X; table(h, band, dets(d), half, f, f - 0.05, p, fg, fER, thr, median(a, 'omitnan'), ...
                 median(s(:), 'omitnan'), median(e), z, pm, nnz(ka), size(s, 1), nnz(kg), ...
                 'VariableNames', {'hemi','band','detector','half','f','e','p','fGrid','fER','thrSurr', ...
                 'subjectMedian','surrMedian','emptyroomMedian','zMedian','pMedian','nFrames','nSurr','nGridFrames'})]; %#ok<AGROW>
        end
    end
end

function f = i_frac(x, thr)
    x = x(isfinite(x));  f = mean(x > thr);
end

function iaf = i_iaf(F, fs)
% individual alpha frequency: the periodic peak (aperiodic fit 2-40 Hz) of the mean sensor PSD in 7-14 Hz
    nf = 2^nextpow2(4*fs);  [P, f] = pwelch(F.', hann(nf), nf/2, nf, fs);  P = mean(P, 2);
    m = f >= 2 & f <= 40;  ap = rheome.spectral.aperiodic(P(m), f(m));
    per = P(m) - ap.Pap;  fm = f(m);  k = fm >= 7 & fm <= 14;
    [pk, i] = max(per(k));  fk = fm(k);  iaf = NaN;  if pk > 0, iaf = fk(i); end
end

%% ---------- G13: speed sweep, dispersion, the two flows
function [T, X] = i_nulls(name, S, F, fs, o)
    N = size(F, 2);  dur = N / fs;  X = table();
    [b, a] = butter(3, o.Band / (fs/2), 'bandpass');
    env = movmean(mean(abs(hilbert(filtfilt(b, a, F.'))), 2)', round(o.WindowS * fs));
    w = round(o.WindowS * fs);  env([1:w, end-w+1:end]) = 0;
    [~, wc] = findpeaks(env, 'MinPeakDistance', w, 'SortStr', 'descend', 'NPeaks', o.NumWindows);
    d = i_dirac(name);
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';  dec = max(1, round(fs/300));  rate = fs/dec;  Lt = round(2*fs);
    tc = round(linspace(Lt, N - Lt, o.NumTiles));
    for h = o.Hemis
        C = i_hemi(S, h, o);
        X = [X; i_peakedness(C, F, fs, wc, w, dur, h, o)]; %#ok<AGROW>
        X = [X; i_dispersion(C, d, S, h, F, fs, o)]; %#ok<AGROW>
        X = [X; i_twoflows(C, Fa, tc, Lt, dec, rate, fs, dur, h)]; %#ok<AGROW>
        fprintf('[patternnulls %s] %s done\n', name, h);
    end
    k = X.half == 0;  lab = X.hemi(k) + "_" + X.statistic(k);
    T = [rheome.scale.rows("patternnulls", repmat("value", sum(k), 1), X.value(k), "", lab);
         rheome.scale.rows("patternnulls", repmat("z", sum(k), 1), X.z(k), "z", lab);
         rheome.scale.rows("patternnulls", repmat("p", sum(k), 1), X.p(k), "p", lab)];
end

function X = i_peakedness(C, F, fs, wc, w, dur, h, o)
    Pk = C.lbo.Phi' * (C.lbo.Mass * C.Kn);  lam = C.lbo.Lambda(:);  sp = linspace(0.05, 1, 60);
    P = cell(1, numel(wc));  fr = [];
    for i = 1:numel(wc)
        idx = wc(i) - floor(w/2) + (0:w-1);
        [ch, fr] = rheome.filters.jspectrum(Pk * F(:, idx), fs);
        P{i} = abs(ch).^2;
    end
    kf = fr > 0 & fr <= 45;  P = cellfun(@(p) p(:, kf), P, 'UniformOutput', false);
    G = arrayfun(@(c) rheome.filters.travwave(lam, fr(kf), c, 1.5), sp, 'UniformOutput', false);
    pk = @(p) i_peak(p, G);
    obs = cellfun(pk, P);  nul = zeros(o.NPerm, numel(wc));
    for q = 1:o.NPerm
        pm = randperm(numel(lam));  nul(q, :) = cellfun(@(p) pk(p(pm, :)), P);
    end
    X = i_nullrows(h, "peakedness", obs, nul, wc/fs, dur, "upper");
end

function v = i_peak(p, G)
    e = cellfun(@(g) sum(p .* g, 'all'), G) / sum(p, 'all');  v = max(e) / mean(e);
end

function X = i_dispersion(C, d, S, h, F, fs, o)
    gv = double(S.B.(char(h)).gv(:))';  rr = reshape((gv - 1)*4 + (2:4)', [], 1);
    cols = find(d.Hemisphere == double(h == "R") + 1);
    Dk = full(d.Phi(rr, cols))' * (d.Mass(rr, rr) * C.K);   % Dirac coefficients of the MNE current [M x nCh]
    lam = abs(d.Lambda(cols));  N = size(F, 2);  X = table();
    nf = 2^nextpow2(8*fs);
    for half = 0:2
        span = 1:N;  if half, span = (1:floor(N/2)) + (half - 1)*floor(N/2); end
        [P, f] = pwelch((Dk * F(:, span)).', hann(nf), nf/2, nf, fs);
        m = f >= 1 & f <= 45;  for hz = [60 120], m = m & ~(f > hz-2 & f < hz+2); end
        ap = rheome.spectral.aperiodic(P(m, :), f(m), struct('knee', true));
        per = P(m, :) - ap.Pap;  fm = f(m);  ka = fm >= 6 & fm <= 16;  fa = fm(ka);
        [mx, i] = max(per(ka, :), [], 1);  ok = mx > 0;  pkf = fa(i(ok));  l = lam(ok);
        rq = @(l_) i_quint(l_, pkf);
        nul = arrayfun(@(~) rq(l(randperm(numel(l)))), 1:o.NPerm)';
        X = [X; i_nullrow(h, "dispersion", half, rq(l), nul, nnz(ok), "upper")]; %#ok<AGROW>
    end
end

function r = i_quint(l, pkf)
    e = prctile(l, [20 80]);  r = median(pkf(l >= e(2))) / median(pkf(l <= e(1)));
end

function X = i_twoflows(C, Fa, tc, Lt, dec, rate, fs, dur, h)
    n = numel(tc);  cv = cell(1, n);  cj = cell(1, n);
    for i = 1:n
        idx = tc(i) - Lt/2 + (0:dec:Lt-1);
        Jc = C.K * Fa(:, idx);
        A = rheome.flow.activation(Jc, Rate=rate, Band=[8 16]);
        Ap = rheome.flow.apparent(A, C.Sh, Rate=rate, Band=[8 16], Alpha=1);
        vJ = rheome.differential.curl(real(Jc), C.Sh, C.fg);
        cv{i} = abs(mean(Ap.vorticity, 2));  cj{i} = abs(mean(vJ(:, 1:end-1), 2));
    end
    R = zeros(n);  for i = 1:n, for j = 1:n, R(i, j) = corr(cv{i}, cj{j}); end, end
    X = table();
    for half = 0:2
        k = true(1, n);  if half, k = (tc/fs < dur/2) == (half == 1); end
        off = R(k, k);  off = off(~eye(nnz(k)));
        X = [X; i_nullrow(h, "rotation_r", half, median(diag(R(k, k))), off, nnz(k), "two")]; %#ok<AGROW>
    end
end

function X = i_nullrows(h, stat, obs, nul, t, dur, side)
% per-window statistic: median over the half's windows, against the same median of each null draw
    X = table();
    for half = 0:2
        k = true(size(t));  if half, k = (t < dur/2) == (half == 1); end
        X = [X; i_nullrow(h, stat, half, median(obs(k)), median(nul(:, k), 2), nnz(k), side)]; %#ok<AGROW>
    end
end

function X = i_nullrow(h, stat, half, v, nul, nU, side)
    nul = nul(isfinite(nul));  mu = mean(nul);  sd = std(nul);
    if side == "upper", p = (1 + sum(nul >= v)) / (1 + numel(nul));
    else, p = (1 + sum(abs(nul - mu) >= abs(v - mu))) / (1 + numel(nul));
    end
    X = table(h, stat, half, v, mu, sd, (v - mu)/sd, p, numel(nul), nU, ...
              'VariableNames', {'hemi','statistic','half','value','nullMean','nullSD','z','p','nNull','nUnits'});
end

%% ---------- shared
function C = i_hemi(S, h, o)
    H = S.B.(char(h));  C.Sh = H.S;  C.lbo = H.lbo;  C.nV = size(H.S.Vertices, 1);
    if ~isfield(C.Sh, 'nV'), C.Sh.nV = C.nV; end
    rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
    C.K = S.Res.ImagingKernel(rows, :);
    lbo = H.lbo;  sg = o.NormalSigmaMM * 1e-3;            % normals smoothed over sigma: the MNE cannot follow sulcal flips
    Nv = lbo.Phi * (exp(-lbo.Lambda(:) * sg^2/2) .* (lbo.Phi' * (lbo.Mass * H.S.VertNormals)));
    Nv = Nv ./ max(vecnorm(Nv, 2, 2), eps);
    C.Kn = Nv(:,1) .* C.K(1:3:end, :) + Nv(:,2) .* C.K(2:3:end, :) + Nv(:,3) .* C.K(3:3:end, :);
    C.heat62 = exp(-lbo.Lambda(:) * (0.062/(2*pi))^2);
    C.fg = rheome.operators.face_gradient(H.S.Vertices, double(H.S.Faces));  C.op = rheome.detect.operator(C.Sh);
    nF = size(H.S.Faces, 1);
    C.Afv = sparse(double(H.S.Faces(:)), repmat((1:nF)', 3, 1), repmat(C.fg.FaceArea, 3, 1), C.nV, nF);
    C.afvsum = full(sum(C.Afv, 2));  C.av = C.afvsum / 3;
end

function C = i_kcols(C, iE)
    C.K = C.K(:, iE);  C.Kn = C.Kn(:, iE);
end

function B = i_band(Fy, N, fs, lim)
% the positive-frequency bins of a band, scaled so (B.F .* r) * exp(2 pi i f t) is the analytic signal
    k = 1 : floor((N-1)/2);  f = k * fs / N;  m = f >= lim(1) & f <= lim(2);
    B = struct('F', 2 * Fy(:, k(m) + 1) / N, 'f', f(m));
end

function Z = i_at(B, t, r)
    Z = (B.F .* r) * exp(2i*pi * B.f(:) * t(:)');
end

function [E, iE] = i_emptyroom(name, S)
% the participant's own empty room on the channels it shares, by name, with S.iSel (as measure_noisefloor)
    N = builtin('load', fullfile(rheome.load.root(), name, 'noise.mat'), 'nrec');  nrec = N.nrec;
    st = rheome.load.study(name);  nm = string(st.chan.Name(S.iSel));  clear st
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;
    [~, iE, jn] = intersect(nm, nn(ok), 'stable');  jo = find(ok);  jn = jo(jn);
    if isempty(iE), error('scale:patterns:noise', 'the empty room shares no good channel with the recording'); end
    E = struct('F', double(nrec.F(jn, :)), 'fs', nrec.sfreq);
end

function d = i_dirac(name)
    try, d = rheome.load.dirac(name);
    catch, rheome.import.dirac(name);  d = rheome.load.dirac(name);
    end
end

% Author: Diellor Basha, 2026
