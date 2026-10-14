classdef tScaleFreqTrack < matlab.unittest.TestCase
% MS1 2.5 through the instrument on the synthetic cached subject (scaleSynthSubject): planted chirps tracked
% across frequency (rheome.scale.measure_freqtrack), the table's shape, and the driver's file.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_freqtrack'
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
        end
    end

    methods (Test)
        function chirpsAreTrackedAndTheNullIsHeldOut(tc)
            S = rheome.scale.sensors(tc.Name);
            [T, X] = rheome.scale.measure_freqtrack(tc.Name, S, Depth=3, JLevel=4, NPer=2, SNRs=20, Hemis="L", ...
                                                    Cases=["still" "up1"], Rates=[0 1]);
            tc.verifyEqual(sum(X.case ~= "real"), 2*2*3);                             % cases x plants x arms
            tc.verifyEqual(sort(unique(X.arm(X.case ~= "real")))', ["control" "noiseless" "snr20"]);
            tc.verifyTrue(all(ismember(["false_rate" "cross_false_rate" "real_frac_pass" "hit_rate" "sweep_ratio" ...
                "control_found" "cross_oct" "false_frac_win" "sweep_err_med" "jump_oct_med" "net_oct_med" ...
                "hit_minus_control" "null_band_occ" "real_score_ratio_med" "null_score_ratio_med" "guard_s"], T.metric)));
            g = T.value(T.metric == "guard_s");
            tc.verifyGreaterThan(g, 0.5, 'the guard covers more than cleanspan''s one-block dilation at 4 Hz');
            tc.verifyLessThan(g, 10, 'the guard stays a few seconds');
            tc.verifyGreaterThan(T.value(T.metric == "null_score_ratio_med"), 0);
            tc.verifyEqual(sum(T.metric == "null_band_occ"), 7);
            rr = T.value(T.metric == "real_band_r");  nr = T.value(T.metric == "null_band_r");
            tc.verifyGreaterThan(rr, 0.2, 'overlapping one-octave bands share envelope');
            tc.verifyEqual(nr, rr, 'AbsTol', 0.05, 'the null keeps the band coupling: one shift per LB mode, shared by every band');
            fr = T.value(T.metric == "false_rate");
            tc.verifyGreaterThanOrEqual(fr, 0);  tc.verifyLessThanOrEqual(T.value(T.metric == "cross_false_rate"), fr);
            q = X(X.arm == "noiseless", :);
            tc.verifyEqual(mean(q.pass), 1, 'a noiseless plant at 20 dB must pass the null');
            tc.verifyLessThan(median(q.octErr), 0.5, 'noiseless: the frequency centroid within half an octave');
            u = q(q.case == "up1", :);
            tc.verifyGreaterThan(median(u.sweepOctS), 0.5, 'noiseless 1 oct/s chirp reads a rising sweep');
            s = q(q.case == "still", :);
            tc.verifyLessThan(median(abs(s.sweepOctS)), 0.3, 'a still tone reads no sweep');
        end

        function anOnsetTransientDoesNotSetTheNull(tc)
            % a red common-mode transient in the first 0.5 s (cleanspan's onset finding; on cf-chirps it was the one passing
            % real window in 78/80 hemispheres and, scattered by the shift, set the null). The guard is the frame count:
            % on this small head the transient alone does not lift the threshold (mutation: only the count fails)
            d = fullfile(getenv('RHEOME_DATA'), 'synth_freqtrack_onset');  scaleSynthSubject(d);
            L = load(fullfile(d, 'study.mat'));  rec = L.rec;  n = round(0.5*rec.sfreq);
            rec.F(:, 1:n) = rec.F(:, 1:n) + 200*std(rec.F(:))*randn(size(rec.F, 1), 1).*sin(2*pi*3*(0:n-1)/rec.sfreq);
            L.rec = rec;  save(fullfile(d, 'study.mat'), '-struct', 'L');
            S = rheome.scale.sensors('synth_freqtrack_onset');
            [~, X] = rheome.scale.measure_freqtrack('synth_freqtrack_onset', S, Depth=3, JLevel=4, NPer=2, SNRs=20, ...
                                                    Hemis="L", Cases="still", Rates=0);
            tc.verifyEqual(mean(X.pass(X.arm == "noiseless")), 1, 'the onset transient must not set the threshold');
            tc.verifyLessThan(sum(X.case == "real"), 20, 'real windows come from the clean span only');
        end

        function theDriverWritesTheTable(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="freqtrack", OutDir=od, Dataset="synth");
            tc.verifyEqual(R.timing.status, "ok", strjoin(R.timing.message));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'freqtrack.csv')));
        end
    end
end
