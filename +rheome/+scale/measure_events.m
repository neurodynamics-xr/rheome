function [T, X, Acf, TF] = measure_events(name, S, opts)
% SCALE.MEASURE_EVENTS  What a detected alpha event is: boundaries from its shape, its spectrum, its timing (MS1 SPEC1 §8).
%
%   [T, X, Acf, TF] = rheome.scale.measure_events(name)
%   [T, X, Acf, TF] = rheome.scale.measure_events(name, S, Depth=8, JLevel=4, NRandom=10)
%
% The events are the tracker's DETECTIONS: per hemisphere and 2 s window, the best Viterbi path that beats
% the mode-shift surrogate threshold (rheome.scale.measure_grouptrack, same seed and calibration, at one
% configuration: depth-8 tiles, 18.75 frames/s by default). Each is one row of an event table (hemisphere,
% tiles visited, first and last time), so the same code reads atlas-cell events once they exist.
% ⭐ A DETECTION IS A PIECE OF AN EPISODE. The tracker cuts at its 2 s windows: on OMEGA sub-0002, 92% of the
% detections spanned their whole window and 70% sat in consecutive windows, so a +-1 s fit on one piece cannot
% see the episode's start or end (18% of fits reached r2 0.7). With Merge (default), detections of one
% hemisphere less than MergeGapS apart whose tiles overlap or touch are one EVENT (tiles: the union).
% The surrogate's detections are merged by the same rule.
%   A. BOUNDARIES   rheome.detect.eventtable: risefall on the event's own envelope (area-weighted mean of
%                   its tiles' 8-16 Hz envelope at 300 frames/s, first-to-last time +- 1 s): tOn, tOff, tauR,
%                   tauF, build-up / plateau / decay, r2 (a fit counts at r2 >= 0.7), with risefall's Model
%                   ("posterior": the quadrature model -- the episode adds to the background in power -- read by its
%                   posterior medians; chosen by measure_eventplants).
%                   The tau bound is revisited: every event is fitted at TauMax (1 s, labelalpha's) and at
%                   TauMaxAlt (2 s). frac_fit_detections (unmerged) and frac_fit_linear are the comparisons.
%   B1. TIME-FREQUENCY  rheome.detect.eventlocked: power at the event's tile centres (minimum-norm kernel
%                   rows, area-weighted), -1..+1 s around tOn and around tOff, 1-45 Hz, each event against
%                   ITS OWN matched random-onset baseline (the same spatial filter at NRandom random times in
%                   the clean span, averaged over time), in dB. Aligned on the fit (fitted events) and on the
%                   tracker's window boundaries (all events: the fallback if the boundary validation fails).
%                   Peak frequency = argmax in Band of the onset-locked ratio over [0, min(dur, 1 s)].
%                   Cycles = duration x peak frequency; cyclesTiles = duration / TimeTileS (a time tile is one
%                   cycle of the band's lowest frequency, 165 ms for alpha) is the cross-check.
%   B2. TIMING      on the BEST boundaries (the fit where it counts, else the tracker's window; "_window":
%                   the window's for all): inter-event intervals (onset to onset, per hemisphere, pooled), CV; the
%                   autocorrelation of the event train (1 while an event is on, per time tile) and of the
%                   hemisphere's mean envelope (per time tile), against the same statistic on the
%                   mode-shift surrogate's detections / envelope and on NPoisson Poisson trains of the same
%                   rate (durations resampled from the real events). The train's rhythmicity is its ACF EXCESS
%                   over the Poisson median at lags >= 2 s (long events alone make the ACF high at short lags),
%                   against the Poisson draws' p95 of the same excess. The envelope's own ACF is kept by the
%                   mode-shift surrogate by construction (each mode keeps its time course): a check, not a test.
%
% OUTPUTS
%   T    rheome.scale.rows (analysis "events"), band "d<depth>/<fps>fps" (and "..._h1"/"_h2", the recording's
%        halves, for split-half ICC of the per-participant summaries)
%   X    events.csv: the event table, one row per real detection (tiles as text), + peakFreqHz cycles
%        cyclesTiles half and the TauMaxAlt refit (tauRAlt tauFAlt r2Alt fitAlt atBoundAlt)
%   Acf  eventacf.csv: lag_s, series, acf (train_fit (best boundaries), train_window, train_surrogate, train_poisson_median,
%        train_poisson_p95, env_real, env_surrogate)
%   TF   eventtf.mat: .f .t, dB maps .on_fit .off_fit .on_window .off_window [nF x nT], .on_fit_h1/_h2,
%        .n_fit .n_window
%
% ⚠ The tracker takes at most ONE detection per 2 s window and per hemisphere, so interval and train
%   statistics at lags below ~2 s are shaped by the detector, not only by the data. The surrogate's
%   detections share that constraint; the Poisson train does not.
% ⚠ These characterise the detector's output; they are not hypothesis tests (SPEC1 §8 B).
% ⚠ The envelope surrogate here is a fresh mode-shift draw (Seed), not grouptrack's own: the same
%   construction, a different realisation.
%
% See also: rheome.detect.eventtable, rheome.detect.eventlocked, rheome.scale.measure_grouptrack,
%           rheome.scale.measure_eventplants, rheome.scale.groupevents
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depth (1,1) double = 8
        opts.JLevel (1,1) double = 4
        opts.Hemis string = ["L" "R"]
        opts.Margin (1,1) double = 1
        opts.Smooth (1,1) double = 0.1
        opts.TauMax (1,1) double = 1
        opts.TauMaxAlt (1,1) double = 2
        opts.R2Min (1,1) double = 0.7
        opts.TimeTileS (1,1) double = 0.165
        opts.Band (1,2) double = [8 16]
        opts.Freqs (:,1) double = (1:45)'
        opts.Cycles (1,1) double = 5
        opts.LockS (1,1) double = 1
        opts.NRandom (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.MaxLagS (1,1) double = 10
        opts.NPoisson (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.Seed (1,1) double = 31
        opts.EdgeS (1,1) double = 10
        opts.Merge (1,1) logical = true
        opts.MergeGapS (1,1) double = 0.25
        opts.Model (1,1) string = "posterior"
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    [~, ~, ~, ~, Pa] = rheome.scale.measure_grouptrack(name, S, Depths=opts.Depth, JLevels=opts.JLevel, ...
                                                       Hemis=opts.Hemis, EdgeS=opts.EdgeS);
    rng(opts.Seed);
    FR = 300;  rate = FR / 2^opts.JLevel;  b = sprintf('d%d/%gfps', opts.Depth, rate);
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);  span = [(CS.first - 1)/fs, CS.last/fs];
    E = rheome.flow.envelopemodes(name, Band=[8 16], Rate=FR);
    Y = struct();  Ys = Y;  A = Y;  Gc = Y;  Ga = Y;
    for hh = opts.Hemis
        h = char(hh);  H = S.B.(h);  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  av = full(sum(Mm, 2));
        Tr = rheome.geom.tree(H.S, L=H.lbo.L, M=Mm, MaxDepth=opts.Depth);
        G = rheome.geom.tiles(Tr, H.S, opts.Depth, Ruler=rheome.geom.edgegraph(H.S));
        P = double(G.P);  Wd = (P'*(av.*Phi)) ./ (P'*av);  A.(h) = P'*av;
        C = double(E.(h).C);  Y.(h) = Wd*C;
        sh = randi([2*FR, size(C, 2) - 2*FR], size(C, 1), 1);
        for k = 1:size(C, 1), C(k, :) = circshift(C(k, :), sh(k)); end
        Ys.(h) = Wd*C;  clear C
        gv = double(H.gv(:));  Gc.(h) = gv(G.centre);  Ga.(h) = G.A > 0 | speye(size(G.A));
    end
    clear E
    R = i_events(Pa(Pa.kind == "real" & Pa.passed, :));  Su = i_events(Pa(Pa.kind == "surrogate" & Pa.passed, :));
    R0 = R;
    if opts.Merge, R = i_merge(R, Ga, opts.MergeGapS);  Su = i_merge(Su, Ga, opts.MergeGapS); end
    n = height(R);
    T = rheome.scale.rows("events", ["n_events" "n_detections" "n_surrogate_events" "span_s"], ...
                          [n height(R0) height(Su) diff(span)], ["events" "detections" "events" "s"], b);
    X = table();  Acf = table();  TF = struct();
    if n == 0, return, end

    % A. boundaries, at the bound and at the revisited bound
    fo = {'Margin', opts.Margin, 'Smooth', opts.Smooth, 'R2Min', opts.R2Min, 'Model', opts.Model};
    X = rheome.detect.eventtable(R, Y, A, FR, fo{:}, TauMax=opts.TauMax);
    X2 = rheome.detect.eventtable(R, Y, A, FR, fo{:}, TauMax=opts.TauMaxAlt);
    X0 = rheome.detect.eventtable(R0, Y, A, FR, fo{:}, TauMax=opts.TauMax);                % unmerged detections
    fo{end} = "linear";  XL = rheome.detect.eventtable(R, Y, A, FR, fo{:}, TauMax=opts.TauMax);
    X.tauRAlt = X2.tauR;  X.tauFAlt = X2.tauF;  X.r2Alt = X2.r2;  X.fitAlt = X2.fit;  X.atBoundAlt = X2.atBound;
    X.half = 1 + (X.t0 >= mean(span));

    % B1. event-locked time-frequency, each event against its own random-onset baseline
    filt = cell(n, 1);
    for i = 1:n
        h = char(X.hemi(i));  k = X.tiles{i}(:);  w = sqrt(A.(h)(k) / sum(A.(h)(k)));
        rows = reshape((Gc.(h)(k)' - 1)*3 + (1:3)', [], 1);
        filt{i} = repelem(w, 3, 1) .* S.Res.ImagingKernel(rows, :);     % (w, 3, 1): a one-tile event stays a column
    end
    lo = {'Window', [-1 1]*opts.LockS, 'Freqs', opts.Freqs, 'Cycles', opts.Cycles};
    pad = opts.LockS + 3*opts.Cycles/(2*pi*min(opts.Freqs)) + 0.1;
    ra = span(1) + pad + rand(n*opts.NRandom, 1) * (diff(span) - 2*pad);
    Pr = rheome.detect.eventlocked(F, fs, ra, repelem(filt, opts.NRandom), lo{:});
    nF = numel(Pr.f);  Bse = squeeze(mean(reshape(Pr.P, nF, [], opts.NRandom, n), [2 3], 'omitnan'));   % [nF x n]
    Bse = reshape(Bse, nF, 1, n);  clear Pr
    fit = X.fit;  on = X.tOn;  off = X.tOff;  on(~fit) = NaN;  off(~fit) = NaN;
    L = struct('on_fit', on, 'off_fit', off, 'on_window', X.t0, 'off_window', X.t1);
    TF.f = opts.Freqs;  TF.n_fit = nnz(fit);  TF.n_window = n;  Rat = struct();
    for a = string(fieldnames(L))'
        Pk = rheome.detect.eventlocked(F, fs, L.(a), filt, lo{:});  TF.t = Pk.t;
        Rat.(a) = Pk.P ./ Bse;  TF.(a) = 10*log10(mean(Rat.(a), 3, 'omitnan'));
    end
    for hf = 1:2, TF.("on_fit_h" + hf) = 10*log10(mean(Rat.on_fit(:, :, X.half == hf), 3, 'omitnan')); end
    fb = TF.f >= opts.Band(1) & TF.f <= opts.Band(2);  fbi = find(fb);
    X.peakFreqHz = nan(n, 1);
    for i = 1:n
        a = "on_window";  d = X.t1(i) - X.t0(i);  if fit(i), a = "on_fit"; d = X.durS(i); end
        r = mean(Rat.(a)(fb, TF.t >= 0 & TF.t <= min(max(d, 1/rate), opts.LockS), i), 2);
        [mx, j] = max(r);  if isfinite(mx), X.peakFreqHz(i) = TF.f(fbi(j)); end
    end
    X.cycles = X.durS .* X.peakFreqHz;  X.cyclesTiles = X.durS / opts.TimeTileS;

    % B2. timing: intervals, event train and envelope autocorrelation
    ed = span(1):opts.TimeTileS:span(2);  nb = numel(ed) - 1;  ml = floor(opts.MaxLagS / opts.TimeTileS);
    lag = (0:ml)' * opts.TimeTileS;
    bo = X.t0;  bf = X.t1;  bo(fit) = X.tOn(fit);  bf(fit) = X.tOff(fit);         % best boundaries: the fit, else the window
    sets = struct('fit', X(:, {'hemi'}), 'window', X(:, {'hemi'}), 'surrogate', Su(:, {'hemi'}));
    sets.fit.on = bo;  sets.fit.off = bf;  sets.window.on = X.t0;  sets.window.off = X.t1;
    sets.surrogate.on = Su.t0;  sets.surrogate.off = Su.t1;
    tr = struct();  cv = struct();  iei = struct();
    for s = string(fieldnames(sets))'
        [tr.(s), cv.(s), iei.(s)] = i_train(sets.(s), opts.Hemis, ed, ml);
    end
    du = bf - bo;  pc = zeros(opts.NPoisson, 1);  pa = zeros(ml + 1, opts.NPoisson);
    for q = 1:opts.NPoisson
        Q = table();
        for hh = opts.Hemis
            m = nnz(X.hemi == hh);  if m == 0, continue, end
            o = span(1) + cumsum(exprnd(diff(span)/m, 2*m + 10, 1));  o = o(o < span(2));
            Q = [Q; table(repmat(hh, numel(o), 1), o, o + du(randi(n, numel(o), 1)), 'VariableNames', {'hemi','on','off'})]; %#ok<AGROW>
        end
        [pa(:, q), pc(q)] = i_train(Q, opts.Hemis, ed, ml);
    end
    ea = zeros(ml + 1, 2);  ix = discretize((0:size(Y.(char(opts.Hemis(1))), 2) - 1)/FR, ed);  ok = ~isnan(ix);
    for hh = opts.Hemis
        h = char(hh);  w = A.(h) / sum(A.(h));
        ea(:, 1) = ea(:, 1) + i_acf(accumarray(ix(ok)', (w'*Y.(h)(:, ok))', [nb 1], @mean), ml) / numel(opts.Hemis);
        ea(:, 2) = ea(:, 2) + i_acf(accumarray(ix(ok)', (w'*Ys.(h)(:, ok))', [nb 1], @mean), ml) / numel(opts.Hemis);
    end
    sr = ["train_fit" "train_window" "train_surrogate" "train_poisson_median" "train_poisson_p95" "env_real" "env_surrogate"];
    V = [tr.fit tr.window tr.surrogate median(pa, 2) prctile(pa, 95, 2) ea];
    Acf = table(repmat(lag, numel(sr), 1), repelem(sr', ml + 1), V(:), 'VariableNames', {'lag_s','series','acf'});

    % metrics
    beyond = lag >= 2;  pm = median(pa, 2);  xs = @(c) max(c(beyond) - pm(beyond));      % ACF excess over Poisson, lags >= 2 s
    dec = @(c) lag(find(c < exp(-1), 1));  e1 = dec(ea(:, 1));  e2 = dec(ea(:, 2));
    if isempty(e1), e1 = NaN; end
    if isempty(e2), e2 = NaN; end
    m = ["frac_fit" "frac_unfitted" "frac_tau_at_bound" "frac_tau_at_bound_alt" "frac_fit_alt" ...
         "dur_s" "tauR_s" "tauF_s" "plateau_s" "frac_plateau0" "buildup_s" "decay_s" "peak_freq_hz" ...
         "cycles" "cycles_tiles" "n_tiles" "area_mm2" "height" "rate_per_min" ...
         "iei_median_s" "iei_cv" "iei_cv_window" "sur_iei_cv" "poisson_iei_cv_median" "poisson_iei_cv_p95" ...
         "train_acf_excess" "train_acf_excess_window" "sur_train_acf_excess" "poisson_train_acf_excess_p95" ...
         "env_acf_decay_s" "sur_env_acf_decay_s"];
    f = X(fit, :);  pq = zeros(opts.NPoisson, 1);  for q = 1:opts.NPoisson, pq(q) = xs(pa(:, q)); end
    v = [mean(fit) 1 - mean(fit) mean(X.atBound(X.r2 >= opts.R2Min)) mean(X.atBoundAlt(X.r2Alt >= opts.R2Min)) mean(X.fitAlt) ...
         median(f.durS) median(f.tauR) median(f.tauF) median(f.plateauS) mean(f.plateauS == 0) median(f.buildupS) median(f.decayS) ...
         median(X.peakFreqHz, 'omitnan') median(f.cycles, 'omitnan') median(f.cyclesTiles) median(X.nTiles) median(X.areaMM2) ...
         median(f.height) 60*n/diff(span) median(iei.fit) cv.fit cv.window cv.surrogate median(pc) prctile(pc, 95) ...
         xs(tr.fit) xs(tr.window) xs(tr.surrogate) prctile(pq, 95) e1 e2];
    u = ["fraction" "fraction" "fraction" "fraction" "fraction" "s" "s" "s" "s" "fraction" "s" "s" "Hz" ...
         "cycles" "tiles" "tiles" "mm2" "A.m" "1/min" "s" "ratio" "ratio" "ratio" "ratio" "ratio" ...
         "r" "r" "r" "r" "s" "s"];
    T = [T; rheome.scale.rows("events", [m "frac_fit_detections" "frac_fit_linear" "n_detections_per_event"], ...
         [v mean(X0.fit) mean(XL.fit) median(X.nDetections)], [u "fraction" "fraction" "detections"], b)];
    tm = ["on_band_rise_db" "on_broad_rise_db" "off_band_fall_db" "off_broad_fall_db" ...
          "win_on_band_rise_db" "win_on_broad_rise_db" "win_off_band_fall_db" "win_off_broad_fall_db"];
    T = [T; rheome.scale.rows("events", tm, [i_tfsum(TF.on_fit, TF, opts.Band, 1) i_tfsum(TF.off_fit, TF, opts.Band, -1) ...
         i_tfsum(TF.on_window, TF, opts.Band, 1) i_tfsum(TF.off_window, TF, opts.Band, -1)], "dB", b)];
    for hf = 1:2                                                  % split halves, for the ICC
        k = X.half == hf;  fk = k & fit;
        [~, ch] = i_train(sets.fit(k, :), opts.Hemis, ed, ml);
        T = [T; rheome.scale.rows("events", ["frac_fit" "dur_s" "tauR_s" "tauF_s" "peak_freq_hz" "cycles" "iei_cv" "on_band_rise_db"], ...
             [mean(fit(k)) median(X.durS(fk)) median(X.tauR(fk)) median(X.tauF(fk)) median(X.peakFreqHz(k), 'omitnan') ...
              median(X.cycles(fk), 'omitnan') ch i_tfsum(TF.("on_fit_h" + hf), TF, opts.Band, 1) * [1; 0]], ...
             ["fraction" "s" "s" "s" "Hz" "cycles" "ratio" "dB"], b + "_h" + hf)]; %#ok<AGROW>
    end
    X.tiles = string(cellfun(@(t) char(strjoin(string(t(:)'), " ")), X.tiles, 'uni', 0));
end

% Detections -> event table rows: hemi, tiles visited, first and last time (s), window.
function Ev = i_events(Q)
    t0 = zeros(height(Q), 1);  t1 = t0;
    for i = 1:height(Q)
        t0(i) = Q.start_s(i) + (Q.frames{i}(1) - 1) / Q.rate(i);  t1(i) = Q.start_s(i) + Q.frames{i}(end) / Q.rate(i);
    end
    Ev = table(Q.hemi, cellfun(@(x) unique(x(:)), Q.tiles, 'uni', 0), t0, t1, Q.window, ones(height(Q), 1), ...
               'VariableNames', {'hemi','tiles','t0','t1','window','nDetections'});
end

% Detections of one hemisphere less than gap s apart, tiles overlapping or touching (Ga), are one event.
function M = i_merge(Ev, Ga, gap)
    Ev = sortrows(Ev, {'hemi', 't0'});  keep = true(height(Ev), 1);  j = 1;
    for i = 2:height(Ev)
        h = char(Ev.hemi(i));
        if Ev.hemi(i) == Ev.hemi(j) && Ev.t0(i) - Ev.t1(j) <= gap && any(Ga.(h)(Ev.tiles{j}, Ev.tiles{i}), 'all')
            Ev.tiles{j} = union(Ev.tiles{j}, Ev.tiles{i});  Ev.t1(j) = max(Ev.t1(j), Ev.t1(i));
            Ev.nDetections(j) = Ev.nDetections(j) + 1;  keep(i) = false;
        else
            j = i;
        end
    end
    M = Ev(keep, :);
end

% Event train per hemisphere (1 while on, per time tile): mean ACF over hemispheres, IEI CV pooled.
function [ac, cv, iei] = i_train(Q, hemis, ed, ml)
    ac = zeros(ml + 1, 1);  iei = [];  nh = 0;  nb = numel(ed) - 1;  dt = ed(2) - ed(1);
    for hh = hemis
        k = find(Q.hemi == hh);
        if numel(k) < 2, continue, end
        x = zeros(nb, 1);
        for j = k(:)'
            x(max(1, floor((Q.on(j) - ed(1))/dt) + 1):min(nb, ceil((Q.off(j) - ed(1))/dt))) = 1;
        end
        ac = ac + i_acf(x, ml);  nh = nh + 1;  iei = [iei; diff(sort(Q.on(k)))]; %#ok<AGROW>
    end
    ac = ac / max(nh, 1);  if nh == 0, ac(:) = NaN; end
    cv = std(iei) / mean(iei);  if numel(iei) < 2, cv = NaN; end
end

function c = i_acf(x, ml)
    x = x(:) - mean(x);
    if ~any(x), c = nan(ml + 1, 1); return, end
    c = xcorr(x, ml, 'coeff');  c = c(ml + 1:end);
end

% Band and broadband (20-45 Hz) change: onset (sgn 1): [0, 0.5] minus [-1, -0.5] s; offset (sgn -1): [-0.5, 0] minus [0.5, 1].
function v = i_tfsum(M, TF, band, sgn)
    fb = TF.f >= band(1) & TF.f <= band(2);  fw = TF.f >= 20;
    a = TF.t >= 0 & TF.t <= 0.5;  p = TF.t >= -1 & TF.t <= -0.5;
    if sgn < 0, a = TF.t >= -0.5 & TF.t <= 0; p = TF.t >= 0.5 & TF.t <= 1; end
    d = @(r) mean(M(r, a), 'all', 'omitnan') - mean(M(r, p), 'all', 'omitnan');
    v = [d(fb) d(fw)];
end

% Author: Diellor Basha, 2026
