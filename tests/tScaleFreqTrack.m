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
                "control_found"], T.metric)));
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

        function theDriverWritesTheTable(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="freqtrack", OutDir=od, Dataset="synth");
            tc.verifyEqual(R.timing.status, "ok", strjoin(R.timing.message));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'freqtrack.csv')));
        end
    end
end
