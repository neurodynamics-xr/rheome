classdef tScaleSlowOsc < matlab.unittest.TestCase
% The MS1 positive control (nsp cf-slowosc: G15, G16): slow oscillations planted travelling anterior -> posterior
% on the synthetic cached subject, read back by measure_slowosc. The framework must find A->P on the original
% events and P->A on the reversed ones; the sensor latency (Massimini's reading) likewise; the group statistics
% must pass the framework and the sensor latency; and the EEG-BIDS reader keeps only the scored N2/N3 epochs.
% The lagged-source null (gate N, rheome.forward.laggedpair): on the travelling plant the framework separates real
% from null; on a second plant of two lagged generators that do not travel, the fit reads the planted lag; on both,
% the null keeps the electrodes' own amplitude ratio and lag.
% bst_of runs only when a Brainstorm checkout is found (RHEOME_BRAINSTORM, or ../brainstorm3 next to this checkout).
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_slowosc'
        HS = []
        T
        X
        E
        XL
        EL
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
            tSleep(fullfile(d, tc.Name));
            est = ["framework" "phasereg" "sensorlatency"];
            b = tc.bst();
            if strlength(b) > 0, setenv('RHEOME_BRAINSTORM', b); est(end+1) = "bst_of"; end
            S = rheome.scale.sensors(tc.Name, Modality="EEG");
            [tc.T, tc.X, tc.E] = rheome.scale.measure_slowosc(tc.Name, S, PosteriorAxis=[0 -1 0], Estimators=est, ...
                HornSchunck=0.01, PatchMM=20, PhaseCentres=60, MaxEvents=12, NormalSigmaMM=5);
            tSleep(fullfile(d, tc.Name), LagS=0.15);                              % the same night from two lagged generators
            S = rheome.scale.sensors(tc.Name, Modality="EEG");
            [~, tc.XL, tc.EL] = rheome.scale.measure_slowosc(tc.Name, S, PosteriorAxis=[0 -1 0], Estimators=["framework" "sensorlatency"], ...
                Conditions=["original" "lagged"], MaxEvents=12, NormalSigmaMM=5);
        end
    end

    methods (Test)
        function detectsThePlantedEvents(tc)
            tc.verifyGreaterThanOrEqual(height(tc.E), 20);                        % 26 planted, some may merge or miss
            tc.verifyEqual(nnz(tc.E.used), 12);
            tc.verifyTrue(all(ismember(tc.E.stage(tc.E.used), ["N2" "N3"])));
        end

        function frameworkFindsAnteriorToPosteriorAndItsReversal(tc)
            [o, r] = tc.mean("framework");
            tc.verifyLessThan(abs(o), 45);
            tc.verifyGreaterThan(abs(r), 135);
        end

        function sensorLatencyFindsAnteriorToPosterior(tc)
            [o, r] = tc.mean("sensorlatency");
            tc.verifyLessThan(abs(o), 30);
            tc.verifyGreaterThan(abs(r), 150);
        end

        function everyEstimatorRunsInEveryCondition(tc)
            k = unique(tc.X(:, {'estimator','condition'}), 'rows');
            tc.verifyEqual(height(k), 4 * numel(unique(tc.X.estimator)));
            tc.verifyTrue(all(isfinite(tc.X.angleDeg(tc.X.estimator ~= "phasereg"))));
            tc.verifyTrue(all(isfinite(tc.X.sourceAPmm(tc.X.estimator == "framework" & tc.X.condition == "original"))));
            tc.verifyTrue(any(tc.T.metric == "vtest_p"));
        end

        function laggedNullSeparatesFromTravel(tc)
            ap = @(c) mean(tc.X.apCos(tc.X.estimator == "framework" & tc.X.condition == c) > 0);
            tc.verifyEqual(ap("original"), 1);
            tc.verifyLessThan(ap("lagged"), 0.5);                              % 0.25 measured
        end

        function laggedFitReadsThePlantedLag(tc)
            e = tc.EL(tc.EL.used, :);
            tc.verifyGreaterThan(median(e.lagS), 0.07);                         % 0.15 planted; 0.10 read (EEG blur
            tc.verifyLessThan(median(e.lagS), 0.2);                             % pulls the back tile forward)
            tc.verifyGreaterThan(median(e.fitR2), 0.85);
            ap = @(c) mean(tc.XL.apCos(tc.XL.estimator == "sensorlatency" & tc.XL.condition == c) > 0);
            tc.verifyGreaterThanOrEqual(ap("lagged"), 0.8);                     % the electrodes see the same front-first lag
            tc.verifyGreaterThanOrEqual(ap("original"), 0.8);
        end

        function laggedNullKeepsTheElectrodesRatioAndLag(tc)
            for e = {tc.E(tc.E.used, :), tc.EL(tc.EL.used, :)}
                q = e{1};
                tc.verifyLessThan(abs(median(q.elecRatioLagged ./ q.elecRatio) - 1), 0.15);
                tc.verifyLessThanOrEqual(median(abs(q.elecLagSLagged - q.elecLagS)), 0.05);
            end
        end

        function groupStatisticsPassTheFramework(tc)
            root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            for s = 1:6                                                         % six "participants": the same events, resampled
                d = fullfile(root, sprintf('sub%02d', s));  mkdir(d);
                ev = unique(tc.X.event);  pick = ev(randi(numel(ev), numel(ev), 1));
                Xs = tc.X(ismember(tc.X.event, pick), :);
                writetable(Xs, fullfile(d, 'slowosc.csv'));
            end
            [G, P] = rheome.scale.groupslowosc(root, fullfile(root, 'group'), MinEvents=5);
            tc.verifyEqual(numel(unique(P.subject)), 6);
            g = G(G.estimator == "framework" & G.condition == "original", :);
            tc.verifyTrue(g.passes);
            tc.verifyEqual(g.frac_ap, 1);
            tc.verifyLessThan(g.ap_lo, 1);                                      % Clopper-Pearson, n = 6
            tc.verifyEqual(g.lagged_n, 6);                                      % gate N: real minus lagged, paired
            tc.verifyGreaterThan(g.lagged_excess, 0);
            tc.verifyTrue(g.beats_lagged);
            tc.verifyTrue(isfile(fullfile(root, 'group', 'group.csv')));
        end

        function plantedMoversReadInBothDirections(tc)
            % gate V on the synthetic head: planted A->P and P->A plane waves on the night's own background between
            % its events. Validated here: the plant table, the levels, and the readings that do hold on this head --
            % the sensor latency at 0.5 and 5 m/s (and its exact flip on the reversed window), the framework at
            % 0.5 m/s. NOT asserted: the framework above ~1 m/s, which follows geometry rather than travel (measured,
            % desk evidence 1a32c6c1), so the plant mode cannot pass a readout that does not hold.
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));  tSleep(fullfile(d, tc.Name));
            S = rheome.scale.sensors(tc.Name, Modality="EEG");
            [T, X, P] = rheome.scale.measure_slowosc(tc.Name, S, Plant=[0.5 5], PlantsPerCell=2, PadS=0.5, ...
                PosteriorAxis=[0 -1 0], Estimators=["framework" "sensorlatency"], Conditions=["original" "reversed"], NormalSigmaMM=5);
            tc.verifyEqual(height(P), 2 * 2 * 2 * 2);                           % speed x direction x level x plants
            tc.verifyEqual(P.p2pUV(P.level == "threshold"), repmat(140, 8, 1));
            tc.verifyGreaterThan(P.p2pUV(P.level == "median"), 140);            % the night's events are above threshold
            tc.verifyTrue(all(isfinite(P.snrDB)));
            tc.verifyEqual(height(X), height(P) * 2 * 3);                       % framework L, R + sensorlatency
            v = @(m, pat) T.value(T.metric == m & contains(T.band, pat));
            tc.verifyEqual(v("frac_correct_sign", "sensorlatency/"), ones(16, 1));            % 2 cond x 2 speed x 2 dir x 2 level
            tc.verifyLessThan(max(v("median_abs_error_deg", "sensorlatency/")), 5);
            f = v("frac_correct_sign", "framework/" + wildcardPattern + "/original/" + wildcardPattern + "/0.5ms/");
            tc.verifyEqual(numel(f), 8);                                        % hemi x direction x level, 2 plants each
            tc.verifyGreaterThanOrEqual(mean(f), 0.75);                         % 15/16 measured: fragile even at 0.5 m/s
        end

        function vtestMatchesItsDefinition(tc)
            [u, p] = rheome.scale.vtest(zeros(8, 1));
            tc.verifyEqual(u, 4, 'AbsTol', 1e-12);                              % 8 * sqrt(2/8)
            tc.verifyEqual(p, 0.5 * erfc(4 / sqrt(2)), 'AbsTol', 1e-15);
        end

        function readerKeepsOnlyScoredNremEpochs(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            fs = 200;  n = 150 * fs;  t = (0:n-1)' / fs;
            sig = [sin(2*pi*1*t), t];                                            % channel B is the time itself
            hdr = edfheader("EDF+");  hdr.NumDataRecords = 150;  hdr.DataRecordDuration = seconds(1);  hdr.NumSignals = 2;
            hdr.SignalLabels = ["A" "B"];  hdr.PhysicalDimensions = ["uV" "uV"];
            hdr.PhysicalMin = [-200 -200];  hdr.PhysicalMax = [200 200];  hdr.DigitalMin = [-32768 -32768];  hdr.DigitalMax = [32767 32767];
            edf = fullfile(d, 'x.edf');  edfwrite(edf, hdr, sig, 'InputSampleType', 'physical');
            ev = fullfile(d, 'events.tsv');  fid = fopen(ev, 'w');
            fprintf(fid, 'onset\tduration\ttrial_type\tvalue\tsample\tartifact_channels\n');
            st = ["W" "N2" "N2" "R" "N3"];  ar = ["n/a" "none" "B" "none" "none"];
            for k = 1:5, fprintf(fid, '%g\t30\tstage/%s\t1\t%d\t%s\n', 30*(k-1), st(k), 30*(k-1)*fs, ar(k)); end
            fclose(fid);
            rec = rheome.io.read.sleepeeg(edf, ev, Rate=100, ChunkS=40);
            tc.verifyEqual(rec.epochs.stage, ["N2"; "N2"; "N3"]);
            tc.verifyEqual(size(rec.F), [2 9000]);
            tc.verifyEqual(rec.seg, [1 4000; 4001 6000; 6001 9000]);           % 60 s run in 40 + 20 s chunks, then N3
            mid = rec.seg(1,1) + 2000;                                          % 20 s into the first N2 epoch
            tc.verifyEqual(double(rec.F(2, mid)), 50, 'AbsTol', 0.05);          % channel B reads the night's time
            tc.verifyEqual(double(rec.F(2, rec.epochs.first(3) + 1000)), 130, 'AbsTol', 0.05);
            tc.verifyEqual(rec.epochs.artifact(2), "B");
        end
    end

    methods
        function [o, r] = mean(tc, est)
            a = @(c) rad2deg(angle(mean(exp(1i * deg2rad(tc.X.angleDeg(tc.X.estimator == est & tc.X.condition == c))), 'omitnan')));
            o = a("original");  r = a("reversed");
        end

        function d = bst(~)
            d = string(getenv('RHEOME_BRAINSTORM'));
            if strlength(d) == 0 || ~isfolder(d)
                d = string(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'brainstorm3'));
            end
            if ~isfile(fullfile(d, 'toolbox', 'math', 'bst_opticalflow.m')), d = ""; end
        end
    end
end

function tSleep(d, opts)
% overwrite the synthetic study with 120 s of sleep EEG: 26 slow oscillations travelling +y -> -y at 0.5 m/s
% in the normal current of both hemispheres, all epochs N2 (the last 30 s N3), units uV.
% LagS: instead, two fixed generators per hemisphere (Gaussian patches, sigma 15 mm, at the +y and -y poles), the
% -y one 0.7 x the +y one's waveform LagS later: nothing travels.
    arguments
        d
        opts.LagS = []
    end
    s = load(fullfile(d, 'study.mat'));  b = load(fullfile(d, 'bases.mat'));  sf = load(fullfile(d, 'surface.mat'));
    fs = 100;  N = 120 * fs;  t = (0:N-1) / fs;  V = sf.S.Vertices;  nV = size(V, 1);
    te = 3 + 4.4 * (0:25) + 0.3 * sin(1:26);  c = 0.5;  y0 = 0.04;
    src = zeros(nV, N);
    w = @(tau) -cos(2*pi*0.8*tau) .* exp(-tau.^2 / (2 * 0.5^2));
    if isempty(opts.LagS)
        for e = te, src = src + w(t - e - (y0 - V(:, 2)) / c); end
    else
        pf = zeros(nV, 1);  pb = zeros(nV, 1);
        for hh = [false true]
            m = find((V(:, 1) > 0) == hh);  [~, a] = max(V(m, 2));  [~, b] = min(V(m, 2));
            pf(m) = exp(-vecnorm(V(m,:) - V(m(a),:), 2, 2).^2 / (2 * 0.015^2));
            pb(m) = exp(-vecnorm(V(m,:) - V(m(b),:), 2, 2).^2 / (2 * 0.015^2));
        end
        for e = te, src = src + pf * w(t - e) + 0.7 * pb * w(t - e - opts.LagS); end
    end
    Nr = sf.S.VertNormals ./ vecnorm(sf.S.VertNormals, 2, 2);  G = s.hm.Gain;
    Gn = G(:, 1:3:end) .* Nr(:, 1)' + G(:, 2:3:end) .* Nr(:, 2)' + G(:, 3:3:end) .* Nr(:, 3)';
    F = Gn * src;  nC = size(F, 1);  F = F - mean(F, 1);
    F = F / max(abs(F(:))) * 250;  F = F + 2 * randn(size(F));
    chan = s.chan;  chan.Type = repmat({'EEG'}, 1, nC);  [chan.Channel.Type] = deal('EEG');
    rec = struct('F', single(F), 'sfreq', fs, 'ChannelFlag', ones(nC, 1), 'ChannelName', {chan.Name(:)}, ...
                 'seg', [1 N], 'epochs', table((0:3)' * 30, ["N2"; "N2"; "N2"; "N3"], repmat("none", 4, 1), (0:3)' * 30 * fs + 1, ...
                 repmat(30 * fs, 4, 1), 'VariableNames', {'onsetS','stage','artifact','first','n'}));
    hm = s.hm;  ncov = struct('NoiseCov', eye(nC), 'FourthMoment', [], 'nSamples', []); %#ok<NASGU>
    save(fullfile(d, 'study.mat'), 'rec', 'chan', 'hm', 'ncov');
    if isfile(fullfile(d, 'noise.mat')), delete(fullfile(d, 'noise.mat')); end
    clear b
end

% Author: Diellor Basha, 2026
