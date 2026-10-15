function [T, X, E] = measure_slowosc(name, S, opts)
% SCALE.MEASURE_SLOWOSC  Sleep slow oscillations as a positive control with known propagation: direction and speed per event, four estimators, two nulls (MS1 G15, G16).
%
%   [T, X, E] = rheome.scale.measure_slowosc(name)              % S = rheome.scale.sensors(name, Modality="EEG")
%   [T, X, E] = rheome.scale.measure_slowosc(name, S, MaxEvents=200, HornSchunck=[0.01 0.1], PosteriorAxis=[-1 0 0])
%   [T, X, P] = rheome.scale.measure_slowosc(name, S, Plant=[1 2 5])   % gate V: planted movers, section 4
%
% Spec: the group-analysis plan, G15 (primary: AnphySleep hd-EEG) and G16.
% Known propagation: slow oscillations travel predominantly front to back at about 1-7 m/s (Massimini et al.
% 2004). The subject is a sleep EEG cache (rheome.scale.importeeg): scored N2/N3 epochs, template forward model.
%
% 1. DETECTION (Massimini 2004 criteria, per channel, average reference, band-pass Band = 0.1-4 Hz, zero phase):
%    a negative-going zero crossing followed by a positive-going one HalfWaveS (0.3-1.0 s) later, a negative
%    peak below NegUV (-80 uV) between them, and a negative-to-positive peak-to-peak above P2PUV (140 uV).
%    Detections on different channels within ClusterS (0.2 s) of the first are one EVENT if at least MinChannels
%    channels detect it; its time is their median negative peak. An event is kept only if its whole analysis
%    window (+-(WindowS + PadS)) lies in one recorded segment and in epochs where no analysed channel carries
%    the release's artefact mark. Up to MaxEvents are drawn with a fixed Seed, BEFORE any direction is
%    computed, and every estimator reads the same events (G16).
% 2. PER EVENT, four CONDITIONS on the same window: original; reversed (the window played backwards: every
%    direction must flip exactly); surrogate (each channel's window phase-randomised independently, which keeps
%    each channel's spectrum and destroys the inter-channel lags: no propagation should be found); lagged (the
%    lagged-source null, rheome.forward.laggedpair: the frontal and the posterior generator that best explain this
%    event -- each one GenDepth tile of normal current per hemisphere, the anterior-third and posterior-third tile
%    carrying the most reconstructed current in the window, through this participant's own gain; its position
%    genFrontAPmm/genBackAPmm the energy-weighted mean of its tiles' -- refitted with no travel, the
%    event's amplitude ratio and lag (|lag| <= MaxLagS), and its own residual phase-randomised as noise). The
%    lagged null is not expected to read as "no direction": two lagged generators can look A->P. A direction
%    claim must beat it (rheome.scale.groupslowosc: real minus lagged A->P fraction, paired over participants).
%    E carries each event's fit (genFrontAPmm, genBackAPmm, fitR2, ratio, lagS) and, with no model, the
%    posterior/frontal electrode thirds' amplitude ratio and lag in the event (elecRatio, elecLagS) and in its null
%    (elecRatioLagged, elecLagSLagged): the null must reproduce them.
% 3. ESTIMATORS (column estimator; one row per event x hemisphere x condition x estimator x param):
%    framework      the reconstructed normal current (whitened MNE, average reference, read along normals
%                   smoothed over NormalSigmaMM) -> its graph-wavelet band map at the member nearest BandMM
%                   (130 mm: readout rule 1, the floor through the instrument, G1) -> tile means at TileDepth ->
%                   the best full-span Viterbi path (rheome.detect.tilepath) at TrackRate, following the band
%                   map's polarity at the event's negative peak. Direction: the displacement between the tile
%                   centres where the path LEAVES its first tile and ARRIVES at its last, projected on the tangent
%                   plane midway; speed: that geodesic distance (tile-centre ruler) over the time between them.
%                   ⚠ A least-squares slope over the whole window reads the frames where the extremum sits
%                   still and under-reads the speed ~6x (synthetic, 0.5 m/s planted); hence leave-to-arrive.
%                   A path that never leaves its tile has no direction (NaN angle, netMM 0, R 0).
%                   ⚠ THE EXTREMUM'S SPEED IS NOT THE PHASE SPEED of a wave much longer than the hemisphere: a slow
%                   oscillation (~1 Hz at 1-7 m/s: 1-7 m) is a ramp across one hemisphere, and its band map's
%                   extremum crosses it at its own pace. Synthetic (tests/tScaleSlowOsc, 0.5 m/s planted, 80 mm
%                   spheres): framework 0.17 m/s, phasereg 1.06, sensor latency 0.70, bst_of 0.01 (direction
%                   within 6 deg for all four). Compare speeds with Massimini through sensorlatency, and read
%                   the framework's speed as the speed of its readout, not of the wave.
%                   ALSO (original condition only): the Helmholtz-band SOURCE -- the extremum of the scalar
%                   potential's band map at BandMM (rheome.differential.helmholtzbands) of the 3-D current at the
%                   event's peak -- and the path's ORIGIN, as signed positions along the posterior axis (sourceAPmm,
%                   originAPmm: negative = anterior of the hemisphere's centroid).
%    bst_of         Brainstorm's bst_opticalflow (rheome.flow.bstopticalflow) on |normal current| over the window,
%                   at each HornSchunck (0.01 = Brainstorm's default); the time-mean velocity field.
%    phasereg       rheome.flow.phaseregression on the analytic normal current (Hilbert of the band-passed sensors,
%                   padded by PadS, then cropped), PatchMM patches at PhaseCentres vertices, wavelengths up to 20 m.
%    sensorlatency  Massimini's own reading, no inverse: the negative-peak latency of each channel that CARRIES the
%                   wave -- a negative half-wave (zero crossings HalfWaveS apart, his criterion 1, no amplitude
%                   threshold) whose negative peak lies within ClusterS of the event; at least MinLatencyChannels,
%                   else no direction -- regressed on the electrode's (posterior, lateral) position; the latency
%                   gradient's direction and 1/|gradient| as speed, R = the latency plane's R^2. One row per event,
%                   hemi "both", param 0 (not NaN: findgroups drops NaN keys).
%                   ⚠ Before 2026-10-11 every channel entered, and about half of them (the ones that see the wave
%                   inverted under the average reference) had their minimum pinned to the window's edge, at -0.2 or
%                   +0.2 s. On AnphySleep the posterior ones were pinned to the EARLY edge more often, which read
%                   P->A (24/25 participants).
%    Field estimators (bst_of, phasereg) give one event direction as the weighted circular mean of the per-vertex
%    directions (weight: amplitude x |velocity|, x the fit's coherence R for phasereg) and their resultant R.
% 4. GATE V (Plant = speeds in m/s; analysis "slowoscplant"): direction accuracy for movers as fast as slow
%    oscillations, in this head. The detected events are REPLACED by planted ones: a plane wave of normal current
%    along PosteriorAxis (A->P and P->A) at each speed, tSleep's SO waveform at PlantHz (0.8 Hz: wavelength
%    1.25-6.25 m at 1-5 m/s, far above 10 cm), through the same gain, scaled so its largest channel peak-to-peak is
%    the median p2p of this participant's detected events ("median") or the detection threshold ("threshold"),
%    added to a window of the participant's own clean N2/N3 EEG at least 2 windows from any detected event. The
%    same conditions and estimators read it. E is then the plant table (speedMS, dirSign, level, p2pUV, sample,
%    snrDB, and the lagged columns); T is per estimator x condition x hemi x speed x direction x level the
%    fraction with the correct A->P sign (no direction = wrong), its one-sided 95% Clopper-Pearson lower bound,
%    the circular mean and the median absolute direction error. On "lagged" it is how often the no-travel null
%    still reads the plant's direction: gate N's validation on planted travel.
%    ⚠ MEASURED LIMITS on the synthetic head (two 80 mm spheres, ~40 dB, 3 plants per cell):
%      framework      right sign at 0.25-0.5 m/s; from 1 m/s on, the same reading whichever way the wave travels
%                     (the wavelength, 1.25 m+, is ~10x the hemisphere: the band map's extremum moves with geometry)
%      phasereg       right sign to 2 m/s (error <= 27 deg); at 5 m/s REVERSED (175 deg) in both hemispheres
%      sensorlatency  right sign at every speed (error ~1 deg), and also on the lagged null: it reads lag, not travel
%    ⚠ MEASURED in the AnphySleep head (one participant, template ICBM152 + OpenMEEG, 3 plants per cell,
%      level "median"): framework at chance (0.33-0.67 correct, median error 55-114 deg) at 1, 2.7 and 5 m/s
%      for plane waves and for plants confined to the medial wall and cingulate (PlantRegion); phasereg 0.67-1.00;
%      sensorlatency (carrying channels) 1.00 within 9 deg on plane waves.
%    ⚠ A plane wave is the simplest mover: it says what the head and inverse do to a known direction, not how
%    real slow oscillations (origin, curvature, bilateral asymmetry) are read.
%    PlantRegion confines the plant to atlas parcels (PlantAtlas, labels matched case-insensitively), e.g. the
%    medial wall and cingulate where hd-EEG source models put slow waves and their anterior-posterior travel
%    (Murphy et al. 2009, PNAS 106:1608): ["superiorfrontal" "anteriorcingulate" "posteriorcingulate"
%    "isthmuscingulate" "precuneus"]. The wave and its amplitude rule are the same; only its extent changes.
% ⭐ DIRECTION IS AN ANGLE IN THE GAUGE, 0 = ANTERIOR -> POSTERIOR. At each vertex the velocity and the
%   PosteriorAxis are projected on the tangent plane and expressed in the hemisphere's trivial gauge
%   (rheome.operators.gauge, +1 at the two faces of largest separation: the MS1 gauge); angleDeg is the velocity's
%   angle minus the posterior axis's, wrapped to (-180, 180]. apCos = cos(angle): +1 A->P, -1 P->A.
% ⚠ PosteriorAxis defaults to [-1 0 0]: Brainstorm's SCS and the AnphySleep electrodes (ALS) both have +x
%   anterior. A protocol in another frame must pass its own.
% ⚠ The detection thresholds are absolute (uV) and depend on the reference; they are Massimini's, on the average
%   reference here. Events are counted per participant, so a reference that shifts amplitudes changes n, not
%   the direction rule.
%
% OUTPUT
%   X  per row: event, hemi, condition, estimator, param, angleDeg, apCos, speedMS, R, netMM, sourceAPmm, originAPmm
%   E  per detected event: sample, nightS, stage, nChannels, negUV, p2pUV, seedChannel, used, and the lagged
%      null's columns (NaN for unused events)
%   T  rheome.scale.rows("slowosc"): detection counts, and per estimator/param x condition x hemi: n, circular
%      mean angle, resultant, fraction A->P, within-participant V-test (towards 0) u and p, median speed
%
% See also: rheome.scale.groupslowosc, rheome.scale.importeeg, rheome.flow.bstopticalflow, rheome.flow.phaseregression,
%           rheome.detect.tilepath, rheome.differential.helmholtzbands, rheome.operators.gauge
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Band (1,2) double = [0.1 4]
        opts.NegUV (1,1) double = -80
        opts.P2PUV (1,1) double = 140
        opts.HalfWaveS (1,2) double = [0.3 1.0]
        opts.MinChannels (1,1) double = 5
        opts.ClusterS (1,1) double = 0.2
        opts.WindowS (1,1) double = 0.4
        opts.PadS (1,1) double = 2
        opts.MaxEvents (1,1) double = 200
        opts.Seed (1,1) double = 71
        opts.BandMM (1,1) double = 130
        opts.TileDepth (1,1) double = 6
        opts.TrackRate (1,1) double = 50
        opts.NormalSigmaMM (1,1) double = 20
        opts.PosteriorAxis (1,3) double = [-1 0 0]
        opts.UpAxis (1,3) double = [0 0 1]
        opts.HornSchunck double = [0.01 0.1]
        opts.PatchMM (1,1) double = 40
        opts.PhaseCentres (1,1) double = 300
        opts.Hemis string = ["L" "R"]
        opts.Conditions string = ["original" "reversed" "surrogate" "lagged"]
        opts.GenDepth (1,1) double = 4
        opts.MaxLagS (1,1) double = 0.3
        opts.Estimators string = ["framework" "bst_of" "phasereg" "sensorlatency"]
        opts.MinLatencyChannels (1,1) double = 8  % sensorlatency: fewer carrying channels -> no direction
        opts.Plant double = []                  % gate V: planted speeds (m/s); replaces the detected events
        opts.PlantsPerCell (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.PlantLevels string {mustBeMember(opts.PlantLevels, ["median" "threshold"])} = ["median" "threshold"]
        opts.PlantHz (1,1) double = 0.8
        opts.PlantRegion string = string.empty  % gate V: atlas labels the plant is confined to (empty: the whole cortex)
        opts.PlantAtlas (1,1) string = "Desikan-Killiany"
    end
    if isempty(S), S = rheome.scale.sensors(name, Modality="EEG"); end
    st = rheome.load.study(name);  rec = st.rec;  fs = rec.sfreq;
    names = string(st.chan.Name(S.iSel));  Loc = reshape(S.Loc, 3, []);
    Y = S.Ref * double(rec.F(S.iSel, :));  clear st
    seg = [1 size(Y, 2)];  if isfield(rec, 'seg'), seg = rec.seg; end
    [b, a] = butter(2, opts.Band / (fs/2));
    for k = 1:size(seg, 1), c = seg(k,1):seg(k,2); Y(:, c) = filtfilt(b, a, Y(:, c).').'; end
    ok = i_cleanmask(rec, names, size(Y, 2));

    E = i_detect(Y, fs, opts);
    L = round((opts.WindowS + opts.PadS) * fs);
    E.used = false(height(E), 1);
    fit = arrayfun(@(t) any(seg(:,1) <= t - L & seg(:,2) >= t + L) && all(ok(max(t-L,1):min(t+L,end))), E.sample);
    cand = find(fit);  rng(opts.Seed);
    if numel(cand) > opts.MaxEvents, cand = sort(cand(randperm(numel(cand), opts.MaxEvents))); end
    E.used(cand) = true;
    if ~isempty(opts.Plant), E = i_plants(E, seg, ok, L, opts);  cand = (1:height(E))'; end
    if isfield(rec, 'epochs') && ~isempty(rec.epochs)
        ep = rec.epochs;  [~, ie] = max((E.sample >= ep.first') & (E.sample < (ep.first + ep.n)'), [], 2);
        E.nightS = ep.onsetS(ie) + (E.sample - ep.first(ie)) / fs;  E.stage = ep.stage(ie);
    else
        E.nightS = (E.sample - 1) / fs;  E.stage = repmat("", height(E), 1);
    end

    post = opts.PosteriorAxis / norm(opts.PosteriorAxis);  lat = cross(post, opts.UpAxis);  lat = lat / norm(lat);
    C = struct([]);
    for h = opts.Hemis, C = [C, i_hemi(S, h, post, opts)]; end %#ok<AGROW>
    if ~isempty(opts.Plant) && ~isempty(opts.PlantRegion), C = i_plantregion(C, name, opts); end
    ap = Loc' * post';  q3 = quantile(ap, [1 2] / 3);  eF = ap < q3(1);  eB = ap > q3(2);   % frontal, posterior thirds
    % (not halves: under the average reference the two halves' means are exact negatives of each other)
    for v = ["genFrontAPmm" "genBackAPmm" "fitR2" "ratio" "lagS" "elecRatio" "elecLagS" "elecRatioLagged" "elecLagSLagged"]
        E.(v) = NaN(height(E), 1);
    end
    rng(opts.Seed + 1);
    w = round(opts.WindowS * fs);  mid = L + 1;  crop = mid - w : mid + w;
    rows = cell(0, 1);
    for ii = 1:numel(cand)
        t = E.sample(cand(ii));  yl0 = Y(:, t-L : t+L);
        if ~isempty(opts.Plant), [yl0, E.snrDB(cand(ii))] = i_plantwin(yl0, E(cand(ii), :), C, fs, L, opts); end
        for cnd = opts.Conditions
            switch cnd
                case "original",  yl = yl0;
                case "reversed",  yl = fliplr(yl0);
                case "surrogate", yl = i_phaserand(yl0);
                case "lagged"
                    [Gp, gFront, gAP, gH] = i_generators(C, yl0);
                    [yl, f] = rheome.forward.laggedpair(yl0, Gp, gFront, fs, MaxLagS=opts.MaxLagS, Group=gH);  e = cand(ii);
                    E{e, ["genFrontAPmm" "genBackAPmm" "fitR2" "ratio" "lagS"]} = ...
                        [gAP(f.front) * f.wFront.^2, gAP(f.back) * f.wBack.^2, f.r2, f.ratio, f.lagS];
                    E{e, ["elecRatio" "elecLagS" "elecRatioLagged" "elecLagSLagged"]} = ...
                        [i_frontback(yl0(:, crop), eF, eB, fs, opts.MaxLagS), i_frontback(yl(:, crop), eF, eB, fs, opts.MaxLagS)];
            end
            za = hilbert(yl.').';  ys = yl(:, crop);  zs = za(:, crop);
            if any(opts.Estimators == "sensorlatency")
                rows{end+1} = i_row(ii, "both", cnd, "sensorlatency", 0, i_latency(yl, fs, Loc, post, lat, opts)); %#ok<AGROW>
            end
            for c = C
                Xn = c.Kn * ys;
                if any(opts.Estimators == "framework")
                    q = i_framework(Xn, c, fs, opts);
                    if cnd == "original"
                        J = c.K3 * ys(:, w + 1);
                        Bh = rheome.differential.helmholtzbands(J, c.Sh, c.lbo, Maps=true);
                        [~, m] = min(abs(Bh.wavelengthMM - opts.BandMM));  [~, v] = max(abs(Bh.PhiBand(:, m, 1)));
                        q.sourceAPmm = 1e3 * (c.V(v,:) - c.cen) * post';
                    end
                    rows{end+1} = i_row(ii, c.h, cnd, "framework", c.bandMM, q); %#ok<AGROW>
                end
                amp = mean(abs(Xn), 2);
                if any(opts.Estimators == "bst_of")
                    for hs = opts.HornSchunck
                        of = rheome.flow.bstopticalflow(Xn, c.Sg, fs, HornSchunck=hs);
                        Vm = mean(of.velocity, 3);  sp = vecnorm(Vm, 2, 2);
                        rows{end+1} = i_row(ii, c.h, cnd, "bst_of", hs, i_fieldq(Vm, amp .* sp, sp, (1:c.nV)', c)); %#ok<AGROW>
                    end
                end
                if any(opts.Estimators == "phasereg")
                    pr = rheome.flow.phaseregression(c.Kn * zs, c.Sg, Rate=fs, Centres=c.pc, RadiusMM=opts.PatchMM, ...
                        Neighbours=c.nbr, LambdaMM=logspace(log10(40), log10(20000), 40));
                    sp = pr.speed;  sp(~isfinite(sp)) = NaN;
                    rows{end+1} = i_row(ii, c.h, cnd, "phasereg", opts.PatchMM, ...
                        i_fieldq(pr.velocity, amp(c.pc) .* pr.R .* isfinite(sp), sp, c.pc, c)); %#ok<AGROW>
                end
            end
        end
        if mod(ii, 25) == 0, fprintf('[slowosc %s] %d/%d events\n', name, ii, numel(cand)); end
    end
    X = vertcat(rows{:});
    if isempty(X), X = i_row(0, "", "", "", NaN, struct()); X(1, :) = []; end
    if isempty(opts.Plant), T = i_metrics(X, E, fs, size(Y, 2)); else, T = i_plantmetrics(X, E); end
end

%% ---------- gate V: planted movers
function P = i_plants(E, seg, ok, L, o)
% speed x direction (+1 A->P, -1 P->A) x amplitude level x PlantsPerCell plants, each on a background window of the
% participant's own clean N2/N3 EEG at least 2L from every detected event (no real slow oscillation under a plant).
% Level "median": the plant's largest channel peak-to-peak is the median p2p of this participant's detected events;
% "threshold": Massimini's P2PUV, the weakest event the detector keeps.
    v = false(size(ok));
    for k = 1:size(seg, 1), v(seg(k,1)+L : seg(k,2)-L) = true; end
    v = v & movmin(double(ok), 2*L + 1) == 1;
    for s = E.sample', v(max(1, s-2*L) : min(end, s+2*L)) = false; end
    bg = find(v);
    assert(~isempty(bg), 'rheome:slowosc:nobackground', 'no clean N2/N3 window of %d samples free of slow oscillations', 2*L + 1);
    amp = dictionary(["median" "threshold"], [median(E.p2pUV) o.P2PUV]);
    assert(isfinite(amp("median")) || ~any(o.PlantLevels == "median"), 'rheome:slowosc:noevents', ...
        'level "median" needs detected events');
    [sp, dr, lv, ~] = ndgrid(o.Plant, [1 -1], 1:numel(o.PlantLevels), 1:o.PlantsPerCell);
    n = numel(sp);  lvl = reshape(o.PlantLevels(lv(:)), [], 1);  t0 = reshape(bg(randi(numel(bg), n, 1)), [], 1);
    P = table((1:n)', sp(:), dr(:), lvl, amp(lvl), t0, NaN(n, 1), true(n, 1), ...
        'VariableNames', {'plant','speedMS','dirSign','level','p2pUV','sample','snrDB','used'});
end

function [y, snr] = i_plantwin(bg, p, C, fs, L, o)
% one planted slow oscillation (tSleep's waveform at PlantHz): a plane wave of normal current along the posterior
% axis at p.speedMS in direction p.dirSign, its trough crossing the cortex's centre at the window's centre, through
% the same gain, scaled so its largest channel peak-to-peak is p.p2pUV, added to the background window bg.
% snr: plant power over background power in the read window (+-WindowS), all channels (dB)
    tt = (-L:L) / fs;  w = @(tau) -cos(2*pi*o.PlantHz*tau) .* exp(-tau.^2 / (2 * 0.5^2));
    cen = mean(vertcat(C.cen), 1);  s = 0;
    for c = C, s = s + c.Gn * (c.plantMask .* w(tt - p.dirSign * ((c.V - cen) * c.post') / p.speedMS)); end
    s = s * p.p2pUV / max(max(s, [], 2) - min(s, [], 2));
    k = L + 1 + (-round(o.WindowS * fs) : round(o.WindowS * fs));
    y = bg + s;  snr = 10 * log10(sum(s(:, k).^2, 'all') / sum(bg(:, k).^2, 'all'));
end

function C = i_plantregion(C, name, o)
% confine the plant to the vertices of the atlas parcels whose label contains any of o.PlantRegion (case-insensitive)
    A = rheome.load.atlas(name, o.PlantAtlas);
    sel = contains(lower(string(A.Label)), lower(o.PlantRegion));
    assert(any(sel), 'rheome:slowosc:plantregion', 'no %s parcel matches %s', o.PlantAtlas, strjoin(o.PlantRegion, ', '));
    for k = 1:numel(C)
        C(k).plantMask = full(any(A.Membership(sel, C(k).gv), 1))' .* C(k).plantMask;   % and on cortex
        assert(any(C(k).plantMask), 'rheome:slowosc:plantregion', 'the plant region is empty in hemisphere %s', C(k).h);
    end
end

function T = i_plantmetrics(X, P)
% per estimator/param x condition x hemi x speed x direction x level. The expected sign is the planted direction
% (reversed: its opposite; lagged: the planted one, so frac_correct_sign there is how often the no-travel null
% still reads the plant's direction). An event with no direction (NaN) counts as WRONG in frac_correct_sign.
%   frac_correct_sign, correct_lo (one-sided 95% Clopper-Pearson lower bound), mean_error_deg (circular mean of
%   the signed error), median_abs_error_deg (finite events), frac_no_direction, median_snr_db
    X = [X, renamevars(P(X.event, {'speedMS','dirSign','level','snrDB'}), 'speedMS', 'plantMS')];   % X.speedMS is the read one
    T = rheome.scale.rows("slowoscplant", ["n_plants" "p2p_median_uV" "p2p_threshold_uV"], ...
        [height(P) max([P.p2pUV(P.level == "median"); NaN]) max([P.p2pUV(P.level == "threshold"); NaN])], ["events" "uV" "uV"]);
    key = unique(X(:, {'estimator','param','condition','hemi','plantMS','dirSign','level'}), 'rows');
    for i = 1:height(key)
        q = key(i, :);
        m = X.estimator == q.estimator & X.condition == q.condition & X.hemi == q.hemi & X.plantMS == q.plantMS & ...
            X.dirSign == q.dirSign & X.level == q.level & (X.param == q.param | (isnan(X.param) & isnan(q.param)));
        ex = q.dirSign * (1 - 2 * (q.condition == "reversed"));
        er = mod(X.angleDeg(m) - 90 * (1 - ex) + 180, 360) - 180;  n = nnz(m);  k = nnz(ex * X.apCos(m) > 0);
        f = isfinite(er);  d = "PA";  if q.dirSign > 0, d = "AP"; end
        b = q.estimator + "/" + string(q.param) + "/" + q.condition + "/" + q.hemi + "/" + q.plantMS + "ms/" + d + "/" + q.level;
        T = [T; rheome.scale.rows("slowoscplant", ["n" "frac_correct_sign" "correct_lo" "mean_error_deg" ...
            "median_abs_error_deg" "frac_no_direction" "median_snr_db"], ...
            [n k/n i_cplo(k, n) rad2deg(angle(mean(exp(1i * deg2rad(er(f)))))) median(abs(er(f))) mean(~f) ...
             median(X.snrDB(m))], ["events" "fraction" "fraction" "deg" "deg" "fraction" "dB"], b)]; %#ok<AGROW>
    end
end

function p = i_cplo(k, n)
% one-sided 95% Clopper-Pearson lower bound of k successes in n
    if k == 0, p = 0; else, p = betainv(0.05, k, n - k + 1); end
end

%% ---------- detection
function E = i_detect(Y, fs, o)
    det = zeros(0, 4);                                           % channel, negative-peak sample, negUV, p2pUV
    for k = 1:size(Y, 1)
        y = Y(k, :);  s = sign(y);  s(s == 0) = 1;
        dn = find(s(1:end-1) > 0 & s(2:end) < 0);  up = find(s(1:end-1) < 0 & s(2:end) > 0);
        for i = 1:numel(dn)
            j = find(up > dn(i), 1);  if isempty(j), break, end
            hw = (up(j) - dn(i)) / fs;  if hw < o.HalfWaveS(1) || hw > o.HalfWaveS(2), continue, end
            [ng, pn] = min(y(dn(i):up(j)));  if ng > o.NegUV, continue, end
            nx = find(dn > up(j), 1);  if isempty(nx), e = min(numel(y), up(j) + round(fs)); else, e = dn(nx); end
            ps = max(y(up(j):e));  if ps - ng < o.P2PUV, continue, end
            det(end+1, :) = [k, dn(i) + pn - 1, ng, ps - ng]; %#ok<AGROW>
        end
    end
    det = sortrows(det, 2);  E = zeros(0, 6);  i = 1;  cw = round(o.ClusterS * fs);
    while i <= size(det, 1)
        g = i:find(det(:, 2) <= det(i, 2) + cw, 1, 'last');
        if numel(unique(det(g, 1))) >= o.MinChannels
            [~, e0] = min(det(g, 2));
            E(end+1, :) = [round(median(det(g, 2))), numel(unique(det(g, 1))), min(det(g, 3)), max(det(g, 4)), det(g(e0), 1), 0]; %#ok<AGROW>
        end
        i = g(end) + 1;
    end
    E = array2table(E(:, 1:5), 'VariableNames', {'sample','nChannels','negUV','p2pUV','seedChannel'});
end

function ok = i_cleanmask(rec, names, nS)
% samples in epochs where no analysed channel is marked (no epoch table: everything is clean)
    ok = true(1, nS);
    if ~isfield(rec, 'epochs') || isempty(rec.epochs), return, end
    ep = rec.epochs;  ok(:) = false;
    for k = 1:height(ep)
        a = ep.artifact(k);
        bad = ~(a == "none" || a == "n/a") && any(ismember(strtrim(split(a, ",")), names));
        if ~bad, ok(ep.first(k) : min(ep.first(k) + ep.n(k) - 1, nS)) = true; end
    end
end

function [Gp, front, ap, grp] = i_generators(C, y)
% the event's candidate generators: per hemisphere, the anterior-third and the posterior-third GenDepth tile carrying
% the most minimum-norm current over the window. Chosen on the readout's own current, not by the best least-squares
% fit: frontal and posterior topographies are near-collinear in EEG, and the best fit of free time courses is a
% cancelling pair that one waveform per generator cannot keep (synthetic, 0.15 s planted: 0.25-0.28 s read; by
% current, 0.10 s).
    Gp = [];  front = false(1, 0);  ap = [];  grp = [];
    for k = 1:numel(C)
        c = C(k);  en = (c.av .* sum((c.Kn * y).^2, 2))' * c.Pg ./ (c.av' * c.Pg);
        f = find(c.thirdF);  [~, i] = max(en(f));  b = find(c.thirdB);  [~, j] = max(en(b));
        Gp = [Gp, c.Gp(:, [f(i) b(j)])];  front = [front, true false];  ap = [ap, c.genAP([f(i) b(j)])];  grp = [grp, k k]; %#ok<AGROW>
    end
end

function q = i_frontback(y, eF, eB, fs, maxLagS)
% the data's own amplitude ratio and lag, no model: posterior over frontal third's mean waveform (rms ratio) and
% the lag of the posterior one behind the frontal one (peak of their cross-correlation within +-maxLagS)
    a = mean(y(eF, :), 1);  b = mean(y(eB, :), 1);
    [xc, lg] = xcorr(b, a, round(maxLagS * fs));  [~, k] = max(abs(xc));
    q = [rms(b) / rms(a), lg(k) / fs];
end

function y = i_phaserand(y)
% each channel's window with its own random Fourier phases (amplitude spectrum kept, cross-channel lags destroyed)
    n = size(y, 2);  h = floor((n - 1) / 2);  Fy = fft(y, [], 2);
    ph = exp(2i * pi * rand(size(y, 1), h));
    Fy(:, 2:h+1) = Fy(:, 2:h+1) .* ph;  Fy(:, n:-1:n-h+1) = Fy(:, n:-1:n-h+1) .* conj(ph);
    y = real(ifft(Fy, [], 2));
end

%% ---------- estimators
function q = i_framework(Xn, c, fs, o)
    B = c.Phi * (c.g .* (c.Phi' * (c.M * Xn)));
    fr = 1 : max(1, round(fs / o.TrackRate)) : size(B, 2);  rate = fs / (fr(2) - fr(1));
    [~, v] = max(abs(B(:, ceil(end/2))));  sg = sign(B(v, ceil(end/2)));  if sg == 0, sg = 1; end
    Yt = rheome.geom.tilemean(c.G, sg * B(:, fr), c.av);
    p = rheome.detect.tilepath(Yt, c.G, Span="full", SampleRate=rate);
    q = i_nanq();
    if ~p.nTracks, return, end
    tl = p.tracks(1).tiles(:);  cv = c.G.centre(tl);
    q.originAPmm = 1e3 * (c.V(cv(1), :) - c.cen) * c.post';  q.netMM = 0;  q.R = 0;
    mv = find(diff(tl) ~= 0);
    if isempty(mv), return, end                                    % the extremum never left its tile: no direction
    f1 = mv(1);  f2 = mv(end) + 1;                                 % last frame at the start tile, first at the end tile
    q.netMM = 1e3 * c.G.D(tl(f1), tl(f2));  q.R = 1;
    q.angleDeg = i_angle(c.V(cv(f2), :) - c.V(cv(f1), :), cv(round((f1 + f2) / 2)), c);
    q.apCos = cosd(q.angleDeg);  q.speedMS = 1e-3 * q.netMM / ((f2 - f1) / rate);
end

function q = i_fieldq(Vv, wt, sp, vtx, c)
% one direction per event from a velocity field at vertices vtx: weighted circular mean of the gauge angles
    q = i_nanq();  wt(~isfinite(wt)) = 0;
    if ~any(wt > 0), return, end
    th = zeros(numel(vtx), 1);
    for i = 1:numel(vtx), th(i) = i_angle(Vv(i, :), vtx(i), c); end
    z = sum(wt .* exp(1i * deg2rad(th))) / sum(wt);
    q.angleDeg = rad2deg(angle(z));  q.apCos = cosd(q.angleDeg);  q.R = abs(z);
    q.speedMS = median(sp(wt > 0), 'omitnan');
end

function [ang, vt] = i_angle(v, k, c)
% the angle of v at vertex k in the gauge, relative to the posterior axis there, in (-180, 180]
    n = c.N(k, :);  vt = v - (v * n') * n;
    a = atan2(vt * c.e2(k, :)', vt * c.e1(k, :)') - c.pAng(k);
    ang = rad2deg(angle(exp(1i * a)));
end

function q = i_latency(yl, fs, Loc, post, lat, o)
% Massimini's delay map on the channels that carry the wave: each channel's negative half-wave (a negative zero
% crossing and the next positive one HalfWaveS apart, Massimini 2004 criterion 1) whose negative peak lies within
% ClusterS of the event centre; that peak's time is the channel's delay. A channel with no such half-wave -- one
% that sees the wave inverted under the average reference, or none -- has no delay. yl: the padded window.
    q = i_nanq();  mid = (size(yl, 2) + 1) / 2;  cw = round(o.ClusterS * fs);  tl = NaN(size(yl, 1), 1);
    for k = 1:size(yl, 1)
        y = yl(k, :);  s = sign(y);  s(s == 0) = 1;
        dn = find(s(1:end-1) > 0 & s(2:end) < 0);  up = find(s(1:end-1) < 0 & s(2:end) > 0);
        for i = 1:numel(dn)
            j = find(up > dn(i), 1);  if isempty(j), break, end
            hw = (up(j) - dn(i)) / fs;  if hw < o.HalfWaveS(1) || hw > o.HalfWaveS(2), continue, end
            [~, pk] = min(y(dn(i):up(j)));  pk = dn(i) + pk - 1;
            if abs(pk - mid) <= cw && ~(abs(pk - mid) >= abs(tl(k) * fs)), tl(k) = (pk - mid) / fs; end
        end
    end
    c = isfinite(tl);
    if nnz(c) < o.MinLatencyChannels, return, end
    A = [ones(nnz(c), 1), Loc(:, c)' * post', Loc(:, c)' * lat'];  tl = tl(c);
    g = A \ tl;  gr = g(2:3);
    if norm(gr) == 0, return, end
    q.angleDeg = rad2deg(atan2(gr(2), gr(1)));  q.apCos = cosd(q.angleDeg);  q.speedMS = 1 / norm(gr);
    r = tl - A * g;  q.R = 1 - sum(r.^2) / max(sum((tl - mean(tl)).^2), eps);       % R^2 of the latency plane
end

function q = i_nanq()
    q = struct('angleDeg', NaN, 'apCos', NaN, 'speedMS', NaN, 'R', NaN, 'netMM', NaN, 'sourceAPmm', NaN, 'originAPmm', NaN);
end

function r = i_row(ev, h, cnd, est, par, q)
    q0 = i_nanq();  f = fieldnames(q0);
    for i = 1:numel(f), if ~isfield(q, f{i}), q.(f{i}) = q0.(f{i}); end, end
    r = table(ev, string(h), string(cnd), string(est), par, q.angleDeg, q.apCos, q.speedMS, q.R, q.netMM, q.sourceAPmm, q.originAPmm, ...
        'VariableNames', {'event','hemi','condition','estimator','param','angleDeg','apCos','speedMS','R','netMM','sourceAPmm','originAPmm'});
end

%% ---------- per hemisphere, once
function c = i_hemi(S, h, post, o)
    H = S.B.(char(h));  Sh = H.S;  c.h = h;  c.Sh = Sh;  c.lbo = H.lbo;
    c.V = double(Sh.Vertices);  F = double(Sh.Faces);  c.nV = size(c.V, 1);  c.cen = mean(c.V, 1);  c.post = post;
    unit = @(A) A ./ max(vecnorm(A, 2, 2), eps);
    lam = H.lbo.Lambda(:);  c.Phi = H.lbo.Phi;  c.M = H.lbo.Mass;  c.av = full(sum(c.M, 2));
    hk = exp(-lam * (o.NormalSigmaMM * 1e-3)^2 / 2);  Nl = unit(Sh.VertNormals);
    Ns = unit(c.Phi * (hk .* (c.Phi' * (c.M * Nl))));
    gv = double(H.gv(:));  K = S.Res.ImagingKernel;  c.gv = gv;  c.plantMask = double(rheome.geom.cortexmask(H.S));   % plants on cortex only
    c.Kn = Ns(:,1) .* K(3*gv-2, :) + Ns(:,2) .* K(3*gv-1, :) + Ns(:,3) .* K(3*gv, :);
    c.K3 = K(reshape((3*gv' - 3) + (1:3)', [], 1), :);
    gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2);
    wl = 1e3 * wavelengths(gfb);  [~, m] = min(abs(wl - o.BandMM));  c.bandMM = wl(m);  c.g = feval(gain(gfb, m), lam);
    ge = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=H.lbo.L, M=c.M, MaxDepth=o.TileDepth);
    c.G = rheome.geom.tiles(Tr, Sh, o.TileDepth, Ruler=ge);
    gg = rheome.operators.gauge(c.V, F, Method="trivial");
    c.N = unit(gg.normal);  c.e1 = gg.e1;  c.e2 = gg.e2;
    pt = post - (c.N * post') .* c.N;  c.pAng = atan2(sum(pt .* c.e2, 2), sum(pt .* c.e1, 2));
    c.Sg = struct('Vertices', c.V, 'Faces', F, 'VertNormals', Nl, 'nV', c.nV, 'Hemi', {{(1:c.nV)'}}, ...
                  'HemiLabel', {{char("Cortex " + h)}}, 'Comment', '', 'SurfaceFile', '');
    c.pc = unique(round(linspace(1, c.nV, min(o.PhaseCentres, c.nV))))';
    Gn = Nl(:,1)' .* S.G(:, 3*gv-2) + Nl(:,2)' .* S.G(:, 3*gv-1) + Nl(:,3)' .* S.G(:, 3*gv);  c.Gn = Gn;   % normal current -> sensors
    if any(o.Conditions == "lagged")                               % candidate generators: GenDepth tiles of normal current
        Tg = rheome.geom.tiles(Tr, Sh, o.GenDepth);
        c.Pg = double(Tg.P);  c.Gp = Gn * (c.Pg .* c.av);  c.genAP = 1e3 * ((Tg.centroid - c.cen) * post')';
        q = quantile(c.genAP, [1 2] / 3);  c.thirdF = c.genAP < q(1);  c.thirdB = c.genAP > q(2);
    else, c.Pg = [];  c.Gp = [];  c.genAP = [];  c.thirdF = logical([]);  c.thirdB = logical([]);
    end
    if any(o.Estimators == "phasereg")
        pr = rheome.flow.phaseregression(zeros(c.nV, 2), c.Sg, Rate=1, Centres=c.pc, RadiusMM=o.PatchMM);  c.nbr = pr.Neighbours;
    else, c.nbr = [];
    end
end

%% ---------- per-participant summary
function T = i_metrics(X, E, fs, nS)
    T = rheome.scale.rows("slowosc", ["n_detected" "n_used" "density_per_min"], ...
        [height(E) nnz(E.used) height(E) / (nS / fs / 60)], ["events" "events" "1/min"]);
    if any(isfinite(E.ratio))                                      % the lagged-source null: its fit, and the data it matches
        v = ["fitR2" "ratio" "lagS" "elecRatio" "elecLagS" "elecRatioLagged" "elecLagSLagged"];
        T = [T; rheome.scale.rows("slowosc", "lagged_" + v + "_median", median(E{:, v}, 1, 'omitnan'), ...
            ["1" "1" "s" "1" "s" "1" "s"])];
    end
    if isempty(X), return, end
    key = unique(X(:, {'estimator','param','condition','hemi'}), 'rows');
    for i = 1:height(key)
        m = X.estimator == key.estimator(i) & X.condition == key.condition(i) & X.hemi == key.hemi(i) & ...
            (X.param == key.param(i) | (isnan(X.param) & isnan(key.param(i))));
        th = deg2rad(X.angleDeg(m));  th = th(isfinite(th));  n = numel(th);
        [u, p] = rheome.scale.vtest(th);
        z = mean(exp(1i * th));
        b = key.estimator(i) + "/" + string(key.param(i)) + "/" + key.condition(i) + "/" + key.hemi(i);
        T = [T; rheome.scale.rows("slowosc", ["n" "mean_angle_deg" "resultant" "frac_ap" "vtest_u" "vtest_p" "median_speed"], ...
            [n rad2deg(angle(z)) abs(z) mean(cos(th) > 0) u p median(X.speedMS(m), 'omitnan')], ...
            ["events" "deg" "1" "fraction" "1" "p" "m/s"], b)]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
