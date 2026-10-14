function [T, X] = measure_freqtrack(name, S, opts)
% SCALE.MEASURE_FREQTRACK  Planted chirps in the participant's own rest, tracked across frequency (MS1 2.5).
%
%   [T, X] = rheome.scale.measure_freqtrack(name)
%   [T, X] = rheome.scale.measure_freqtrack(name, S, Depth=8, JLevel=4, NPer=10, SNRs=[20 10])
%
% The through-MEG test of rheome.detect.jointpath, in the injection design of rheome.scale.measure_inject (same
% source, segment and readout), with frequency, not place, as the thing that moves:
%   bands     one octave wide, centred at 2.^LogF (half-octave steps, 4-32 Hz: 2.8-45 Hz); per band, 8th-order
%             band-pass and Hilbert at the sensors, minimum norm, |J| per vertex, LB modes, depth-Depth tile
%             means (rheome.flow.envelopemodes' chain), 300 frames/s block-averaged by 2^JLevel (the time tile)
%   source    a Laplace-Beltrami heat-kernel blob (SigmaMM 15) along the normal at a random vertex, forwarded
%             through the participant's leadfield; it stays put -- the place is not what is tested
%   cases     still (FcHz), hop (a CROSSFADE of two tones at 2.^HopLogF: both present mid-window; nothing drifts), and
%             exponential chirps at Rates octaves/s (log2 f = log2 FcHz + r (t - W/2), held at the window's ends);
%             NPer 2 s windows each from the clean span, present for the window +-1 s with 0.25 s ramps
%   amplitude per-band SNR (dB) = 20 log10 of the plant's full-height envelope on its tile in the band of its
%             instantaneous frequency (a unit tone read alone) over the null's median there (N.norm: the real
%             envelope's median over the calibration windows, circularly shifted); the chirp's
%             amplitude follows that gain in log f, so the SNR is the same at every frequency it passes
%   arms      noiseless (the plant alone at 20 dB: the instrument's own bias), each SNR in the own rest, and
%             control (the same rest segment without the plant)
%   null      each LB mode's envelope circularly shifted by ONE shift shared by every band (keeps each series'
%             spectrum and the coupling of overlapping bands; destroys alignment across places);
%             rheome.detect.jointpathnull, odd windows calibrate, even windows test: Alpha false paths/s, the
%             cross-band criterion (the 95th percentile of |netOct| over passing null paths, or CrossOct if
%             given) and the CROSS-BAND false rate at it
%   tracker   rheome.detect.jointpath, top 3 paths; FOUND = a passing path with >= half its frames within NearMM
%             of the plant and a MEAN |frequency centroid - truth| <= OctTol octaves over its frames (the truth
%             clamped to the band set's centres, which bound the readout). The threshold is calibrated on the
%             best path per window and the plants may use any of 3: control_found measures that specificity,
%             and hit_minus_control is the hit rate above it
%   real      every 2 s window of the clean span (rheome.scale.cleanspan: edges and bad blocks out, then a guard of
%             the lowest band's two-sided envelope response down GuardDB, guard_s; the null built on the same
%             frames) through the same tracker and threshold (case "real")
%
% X (freqtrack.csv): one row per hemisphere x case x plant x arm (and per real window): truth, readout (pass found
%   near octErr sweepOctS octSpan durS netMM netOct jumpOct sweepErr scoreRatio), the half of the recording, the seed, the
%   null's held-out false rates (per s and per window), its cross criterion, the number of passing calibration
%   paths that set it (nullCalN) and band occupancy, and the median adjacent-band envelope correlation of the real
%   field and of the null (realBandR, nullBandR: equal under the shared shift, the null's guard), the median null
%   calibration score over the threshold (nullScoreRatio) and the guard (guardS, s). sweepErr = the
%   read sweep - the LS slope of the (clamped) truth over the path's own frames: a perfect readout reads 0 in
%   every case (the hop's step truth has an LS slope of ~1.1 oct/s, not 0).
% T (rheome.scale.rows, analysis "freqtrack"), band "<arm>/<case>/d<depth>/<fps>fps":
%   hit_rate        fraction of plants FOUND                                          fraction
%   sweep_ratio     median read sweep / planted rate, found chirps ⭐                 ratio
%   sweep_oct_s     median read sweep of found plants (still: 0; hop: see sweep_err)  oct/s
%   sweep_err_med   median sweepErr of found plants (0 for a perfect readout)        oct/s
%   net_oct_med     median netOct of found plants                                     oct
%   jump_oct_med    median jumpOct of found plants (hop: one large step; chirp: small) oct
%   oct_err_med     median frequency-centroid error of found plants                   oct
%   hit_minus_control  hit_rate - control_found of the same case                      fraction
%   control_found   (arm control) fraction of plant-free segments with a FOUND path    fraction
% and band "d<depth>/<fps>fps": false_rate, cross_false_rate (1/s, held-out surrogate), false_frac_win,
% cross_frac_win (the same per 2 s window, the unit of real_frac_*), cross_oct (the criterion), cross_cal_n (passing
% calibration paths behind it), real_band_r and null_band_r (adjacent-band envelope correlation), real_frac_pass,
% real_frac_cross (passing real windows with |netOct| >= cross_oct), real_sweep_abs_med (oct/s),
% real_score_ratio_med and null_score_ratio_med (median best-path score / threshold of the real and the null
% calibration windows: equal if the null is a null of this recording; real below null = a conservative null),
% guard_s (s); per band
% "d<depth>/<fps>fps/<Hz>Hz": null_band_occ; hit_rate per half ("_h1" / "_h2"). All medians over hemispheres
% of ONE participant: the per-participant false rate.
%
% ⚠ The frequency readout is the centroid over the path's band and its neighbours, so it is bounded by the
%   band set: below 2^(LogF(1)) or above 2^(LogF(end)) a chirp reads clamped, and sweep_ratio < 1 there.
% ⚠ On noise the best path's centroid spans >= 1 octave in nearly every passing window: a path that crosses
%   bands is evidence only by |netOct| against cross_oct and cross_false_rate, never by its span.
% ⚠ The null's paths sit mostly in the low bands (null_band_occ): the threshold is set by them, and the
%   sensitivity differs by band.
% ⚠ The per-mode shift also removes the minimum norm's spatial point-spread correlation between places (as
%   measure_inject and grouptrack): the null is conservative about place, not about band.
%
% See also: rheome.detect.jointpath, rheome.detect.jointpathnull, rheome.scale.measure_inject
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depth (1,1) double = 8
        opts.JLevel (1,1) double = 4
        opts.LogF double = 2:0.5:5
        opts.BandOct (1,1) double = 1
        opts.FcHz (1,1) double = 10
        opts.Cases string = ["still" "hop" "up05" "up1" "up2" "down1"]
        opts.Rates double = [0 0 0.5 1 2 -1]
        opts.HopLogF (1,2) double = [2.5 4]
        opts.SNRs double = [20 10]
        opts.NPer (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.SigmaMM (1,1) double = 15
        opts.NearMM (1,1) double = 50
        opts.OctTol (1,1) double = 0.5
        opts.Alpha (1,1) double = 0.1
        opts.CrossOct double {mustBeScalarOrEmpty, mustBePositive} = []
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 53
        opts.EdgeS (1,1) double = 10
        opts.GuardDB (1,1) double {mustBePositive} = 60
    end
    if numel(opts.Rates) ~= numel(opts.Cases), error('scale:freqtrack:cases', 'Cases and Rates must pair up.'); end
    if isempty(S), S = rheome.scale.sensors(name); end
    rng(opts.Seed);
    FR = 300;  W = 2;  MARGIN = 1;  PAD = 2;  bs = 2^opts.JLevel;  nF = floor(W*FR/bs);  rate = FR/bs;
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    lf = opts.LogF(:)';  nB = numel(lf);  ed = 2.^(lf' + [-1 1]*opts.BandOct/2);
    if ed(end) >= 0.45*fs, error('scale:freqtrack:nyquist', 'Top band edge %.1f Hz is too close to fs/2 = %g.', ed(end), fs/2); end
    bp = cell(1, nB);  E = bp;
    for b = 1:nB
        bp{b} = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', ed(b, 1), 'HalfPowerFrequency2', ed(b, 2), 'SampleRate', fs);
        E{b} = rheome.flow.envelopemodes(name, Band=ed(b, :), Rate=FR);
    end
    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);
    % the envelope frames inside the clean span and off its bad blocks: the recording's onset transient (1-4 Hz at
    % ~200x) otherwise passes as the one real path and, scattered by the per-mode shift, sets the null's threshold.
    % The envelopes are filtered over the whole recording, so a guard of the lowest band's two-sided (filtfilt +
    % Hilbert) response, down GuardDB, keeps their ringing out too (cf-chirps calibration 2: with the clean span
    % alone, the first or last window of it was still the one passing real window in 57/80 hemispheres)
    gS = i_guard(bp{1}, fs, opts.GuardDB);
    cFirst = CS.first + gS;  cLast = CS.last - gS;  badG = movmax(CS.mask, 2*gS + 1);
    sFr = 1 + (0:size(E{1}.(char(opts.Hemis(1))).C, 2) - 1)*round(fs/FR);
    keepFr = sFr >= cFirst & sFr <= cLast & ~badG(min(sFr, numel(badG)));
    nSeg = (W + 2*MARGIN + 2*PAD)*fs;  tSeg = (0:nSeg-1)/fs - (MARGIN + PAD);
    env = i_env(tSeg, W, MARGIN, 0.25);  cIdx = (MARGIN + PAD)*fs + (1:round(fs/FR):W*fs);
    tFr = ((0:nF-1)*bs + bs/2)/FR;                                  % tracker frame centres in the window (s)
    uW = min(max(tSeg/W, 0), 1);
    X = table();
    for hh = opts.Hemis
        h = char(hh);  H = S.B.(h);  Sh = H.S;  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  lam = H.lbo.Lambda(:);
        av = full(sum(Mm, 2)) .* rheome.geom.cortexmask(H.S);  gv = double(H.gv(:));  nV = numel(gv);  nrm = double(Sh.VertNormals);
        rows = reshape((gv' - 1)*3 + (1:3)', [], 1);  K = S.Res.ImagingKernel(rows, :);  GL = S.G(:, rows);
        hk = exp(-lam*(opts.SigmaMM*1e-3)^2/2);
        ge = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=H.lbo.L, M=Mm, MaxDepth=opts.Depth);
        G = rheome.geom.tiles(Tr, Sh, opts.Depth, Ruler=ge);  P = double(G.P);  nTile = size(P, 2);
        Wd = (P'*(av.*Phi))./(P'*av);  Wd(~isfinite(Wd)) = 0;   % a tile wholly on the medial wall: no cortex, reads 0 (weight 0)  Dcv = distances(ge, G.centre);
        nFr = nnz(keepFr);  nWin = floor(nFr/(W*FR));
        Yr = zeros(nTile, nB, nF*nWin);  Ys = Yr;
        sh = randi([W*FR, nFr - W*FR], size(E{1}.(h).C, 1), 1);      % one shift per LB mode, shared by every band
        for b = 1:nB                                                % one band at a time: never 7 full-rate copies
            C = double(E{b}.(h).C(:, keepFr));
            Yr(:, b, :) = reshape(i_block(Wd*C, bs, W*FR, nWin), nTile, 1, []);
            for k = 1:size(C, 1), C(k, :) = circshift(C(k, :), sh(k)); end
            Ys(:, b, :) = reshape(i_block(Wd*C, bs, W*FR, nWin), nTile, 1, []);
        end
        clear C
        bR = [i_bandr(Yr) i_bandr(Ys)];                             % a per-band shift drives the null's towards 0
        N = rheome.detect.jointpathnull(Ys, G, nF, rate, LogF=lf, Alpha=opts.Alpha, WindowS=W, CrossOct=opts.CrossOct);
        clear Ys
        base = table(hh, opts.Depth, median(G.diameterMM), rate, opts.Seed, N.falseTestS, N.crossTestS, ...
                     N.falseTestW, N.crossTestW, N.crossOct, numel(N.calNet), N.bandOcc, bR(1), bR(2), ...
                     median(N.calScores)/N.threshold, gS/fs, ...
                     'VariableNames', {'hemi','depth','tileMM','rate','seed','nullFalseS','nullCrossS','nullFalseW', ...
                     'nullCrossW','nullCrossOct','nullCalN','nullOcc','realBandR','nullBandR','nullScoreRatio','guardS'});
        for w = 1:nWin
            r = i_path(Yr(:, :, (w-1)*nF + (1:nF)), G, N, lf, rate, [], [], opts);
            X = [X; base, i_row("real", NaN, w, 1 + (w > nWin/2), "real", r)]; %#ok<AGROW>
        end
        clear Yr
        rd = @(Fs, bb) i_read(Fs, bp, bb, cIdx, K, Phi, Mm, Wd, bs, nF);
        lims = [max(cFirst, 1), floor((cFirst + cLast - nSeg)/2); floor((cFirst + cLast - nSeg)/2) + 1, cLast - nSeg + 1];
        for ci = 1:numel(opts.Cases)
            cs = opts.Cases(ci);  r0 = opts.Rates(ci);
            for n = 1:opts.NPer
                hf = 1 + (n > opts.NPer/2);  s0 = randi(lims(hf, :)) + (0:nSeg-1);
                cvx = find(rheome.geom.cortexmask(Sh));  v0 = cvx(randi(numel(cvx)));  k0 = G.tileOf(v0);  g0 = GL*i_blob(v0, Phi, Mm, hk, nrm);   % on cortex
                gain = zeros(1, nB);                               % per band: real median / unit-tone envelope
                for b = 1:nB
                    Yu = rd(g0*(env.*cos(2*pi*2^lf(b)*tSeg)), b);
                    gain(b) = N.norm(k0, b)/max(Yu(k0, b, :));
                end
                ga = @(l) exp(interp1(lf, log(gain), min(max(l, lf(1)), lf(end))));
                if cs == "hop"
                    car = ga(opts.HopLogF(1))*cos(pi/2*uW).*cos(2*pi*2^opts.HopLogF(1)*tSeg + 2*pi*rand) + ...
                          ga(opts.HopLogF(2))*sin(pi/2*uW).*cos(2*pi*2^opts.HopLogF(2)*tSeg + 2*pi*rand);
                    lT = opts.HopLogF(1) + diff(opts.HopLogF)*(tFr/W > 0.5);
                else
                    l = log2(opts.FcHz) + r0*(min(max(tSeg, 0), W) - W/2);
                    car = ga(l).*cos(2*pi*cumsum(2.^l)/fs + 2*pi*rand);
                    lT = log2(opts.FcHz) + r0*(tFr - W/2);
                end
                bsig = g0*(env.*car);  tru = {Dcv(:, v0), lT};
                arms = ["control" "noiseless" "snr" + opts.SNRs];
                for a = arms
                    switch a
                        case "control",   Fs = F(:, s0);
                        case "noiseless", Fs = 10^(20/20)*bsig;
                        otherwise,        Fs = F(:, s0) + 10^(sscanf(a, 'snr%f')/20)*bsig;
                    end
                    r = i_path(rd(Fs, 1:nB), G, N, lf, rate, tru{:}, opts);
                    X = [X; base, i_row(cs, r0, n, hf, a, r)]; %#ok<AGROW>
                end
            end
            fprintf('[freqtrack %s] %s %-6s: %d plants\n', name, h, cs, opts.NPer);
        end
    end
    T = i_metrics(X, opts);
end

function r = i_row(cs, r0, n, hf, arm, v)
    r = [table(cs, r0, n, hf, arm, 'VariableNames', {'case','rate_oct_s','n','half','arm'}), ...
         array2table(v, 'VariableNames', {'pass','found','near','octErr','sweepOctS','octSpan','durS','netMM', ...
                                          'netOct','jumpOct','sweepErr','scoreRatio'})];
end

function T = i_metrics(X, opts)
    T = table();  r = @(m, v, u, b) rheome.scale.rows("freqtrack", m, v, u, b);
    b0 = sprintf('d%d/%gfps', opts.Depth, X.rate(1));
    q = X(X.case == "real", :);  p = q(q.pass == 1, :);
    T = [T; r(["false_rate" "cross_false_rate" "false_frac_win" "cross_frac_win" "cross_oct" "cross_cal_n" ...
               "real_band_r" "null_band_r" "real_frac_pass" "real_frac_cross" "real_sweep_abs_med" ...
               "real_score_ratio_med" "null_score_ratio_med" "guard_s"], ...
              [median(q.nullFalseS) median(q.nullCrossS) median(q.nullFalseW) median(q.nullCrossW) median(q.nullCrossOct) ...
               median(q.nullCalN) median(q.realBandR) median(q.nullBandR) ...
               mean(q.pass) mean(q.pass == 1 & abs(q.netOct) >= q.nullCrossOct) median(abs(p.sweepOctS), 'omitnan') ...
               median(q.scoreRatio) median(q.nullScoreRatio) median(q.guardS)], ...
              ["1/s" "1/s" "fraction" "fraction" "oct" "count" "r" "r" "fraction" "fraction" "oct/s" ...
               "ratio" "ratio" "s"], b0)];
    T = [T; r(repmat("null_band_occ", numel(opts.LogF), 1), median(q.nullOcc, 1), "fraction", b0 + "/" + compose("%.3g", 2.^opts.LogF(:)) + "Hz")];
    for a = unique(X.arm(X.case ~= "real"))'
        for c = opts.Cases
            q = X(X.arm == a & X.case == c, :);  f = q(q.found == 1, :);  bc = a + "/" + c + "/" + b0;
            if a == "control"
                T = [T; r("control_found", mean(q.found), "fraction", bc)]; %#ok<AGROW>
                continue
            end
            sr = f.sweepOctS ./ f.rate_oct_s;  sr(f.rate_oct_s == 0) = NaN;
            cf = mean(X.found(X.arm == "control" & X.case == c));
            T = [T; r(["hit_rate" "sweep_ratio" "sweep_oct_s" "sweep_err_med" "net_oct_med" "jump_oct_med" "oct_err_med" ...
                       "hit_minus_control"], ...
                      [mean(q.found) median(sr, 'omitnan') median(f.sweepOctS, 'omitnan') median(f.sweepErr, 'omitnan') ...
                       median(f.netOct, 'omitnan') median(f.jumpOct, 'omitnan') median(f.octErr, 'omitnan') mean(q.found) - cf], ...
                      ["fraction" "ratio" "oct/s" "oct/s" "oct" "oct" "oct" "fraction"], bc)]; %#ok<AGROW>
            for hf = 1:2
                T = [T; r("hit_rate", mean(q.found(q.half == hf)), "fraction", bc + "_h" + hf)]; %#ok<AGROW>
            end
        end
    end
end

% Top 3 disjoint paths (1 without a truth); FOUND = passing, near the plant for >= half its frames, and its
% frequency centroid within OctTol of the truth clamped to the band centres (mean over its frames). Returns
% [pass found near octErr sweepOctS octSpan durS netMM netOct jumpOct sweepErr scoreRatio] of the found path, else of
% the best one (scoreRatio = its score / the null's threshold).
function r = i_path(Y, G, N, lf, rate, dv, lT, opts)
    r = [0 0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0];
    p = rheome.detect.jointpath(Y, G, Null=N, LogF=lf, SampleRate=rate, NumPaths=3 - 2*isempty(dv));
    if p.nTracks == 0, return; end
    ok = [p.tracks.score] > N.threshold;  nr = nan(p.nTracks, 1);  oe = nr;  q = 1;
    if ~isempty(dv)
        lT = min(max(lT, lf(1)), lf(end));                         % the readout cannot leave the band centres
        for i = 1:p.nTracks
            x = p.tracks(i);  nr(i) = mean(dv(x.tiles) <= opts.NearMM*1e-3);
            oe(i) = mean(abs(x.logf(:) - lT(x.frames(:))'));
        end
        fnd = ok(:) & nr >= 0.5 & oe <= opts.OctTol;
        if any(fnd), [~, q] = max(nr.*fnd - 1e-3*oe.*fnd); end
    else
        fnd = false;
    end
    x = p.tracks(q);  one = p;  one.tracks = x;  one.nTracks = 1;  s = rheome.detect.trackstats(one, G, rate);
    se = NaN;
    if ~isempty(dv) && numel(x.frames) >= 3, pt = polyfit(x.frames(:)/rate, lT(x.frames(:))', 1);  se = x.sweepOctS - pt(1); end
    r = [ok(q) any(fnd) nr(q) oe(q) x.sweepOctS x.octSpan s.durationS s.netMM x.netOct x.jumpOct se x.score/N.threshold];
end

% The band envelopes of a sensor segment, tile means at the tracker's frames: [nTile x nB x nF] (bands not in bb
% stay 0). The chain of rheome.flow.envelopemodes on one segment, as measure_inject's readout.
function Y = i_read(Fs, bp, bb, cIdx, K, Phi, Mm, Wd, bs, nF)
    Y = zeros(size(Wd, 1), numel(bp), nF);
    for b = bb
        Fa = hilbert(filtfilt(bp{b}, Fs.')).';
        Cw = Phi'*(Mm*rheome.flow.activation(K*Fa(:, cIdx)));
        Y(:, b, :) = reshape(i_block(Wd*Cw, bs, numel(cIdx), 1), [], 1, nF);
    end
end

% Half-width (samples) of a band's two-sided envelope response (filtfilt, then Hilbert) above -db of its peak.
function g = i_guard(bp, fs, db)
    c = round(30*fs);  x = zeros(2*c - 1, 1);  x(c) = 1;
    a = abs(hilbert(filtfilt(bp, x)));
    g = max(abs(find(a > 10^(-db/20)*max(a)) - c));
end

% Median over tiles and adjacent band pairs of the envelope correlation over time, Y [nTile x nB x nT].
function m = i_bandr(Y)
    Z = Y - mean(Y, 3);  Z = Z./sqrt(sum(Z.^2, 3));
    m = median(sum(Z(:, 1:end-1, :).*Z(:, 2:end, :), 3), 'all', 'omitnan');
end

% Heat-kernel blob at vertex v along the normal, interleaved [x1;y1;z1;x2;...] (as measure_inject).
function j = i_blob(v, Phi, Mm, hk, nrm)
    d = zeros(size(Phi, 1), 1);  d(v) = 1;
    a = Phi*(hk.*(Phi'*(Mm*d)));  a = a/max(a);
    j = reshape((a.*nrm)', [], 1);
end

% 1 inside [-margin, W+margin] with cosine ramps of r seconds at both outer ends (as measure_inject).
function e = i_env(t, W, margin, r)
    e = ones(size(t));  a = -margin;  b = W + margin;
    e(t < a | t > b) = 0;
    up = t >= a & t < a + r;  e(up) = 0.5 - 0.5*cos(pi*(t(up) - a)/r);
    dn = t > b - r & t <= b;  e(dn) = 0.5 - 0.5*cos(pi*(b - t(dn))/r);
end

% nWin windows of nW frames, each block-averaged by bs: [nRow x floor(nW/bs)*nWin].
function Y = i_block(Y0, bs, nW, nWin)
    nF = floor(nW/bs);  nT = size(Y0, 1);
    Z = reshape(Y0(:, 1:nW*nWin), nT, nW, nWin);  Z = Z(:, 1:nF*bs, :);
    Y = reshape(mean(reshape(Z, nT, bs, nF, nWin), 2), nT, nF*nWin);
end

% Author: Diellor Basha, 2026
