classdef tScaleEvents < matlab.unittest.TestCase
% MS1 SPEC1 §8 end to end on a synthetic cached subject (scaleSynthSubject: two icosphere hemispheres, a random
% leadfield, 40 s of noise with a 10 Hz source): the event table, its TF and timing, the planted boundaries,
% the driver's files and the group reduction.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_events'
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
        end
    end

    methods (Test)
        function eventsHaveBoundariesSpectraAndTiming(tc)
            S = rheome.scale.sensors(tc.Name);
            [T, X, Acf, TF] = rheome.scale.measure_events(tc.Name, S, Depth=3, JLevel=2, Hemis="L", ...
                                                          NRandom=2, NPoisson=5, MaxLagS=4);
            n = T.value(T.metric == "n_events");
            tc.assumeGreaterThan(n, 0);
            tc.verifyEqual(height(X), n);
            tc.verifyLessThanOrEqual(n, T.value(T.metric == "n_detections"));          % merging never adds events
            tc.verifyTrue(all(ismember(["tOn" "tOff" "tauR" "r2" "fit" "peakFreqHz" "cycles" "tiles"], X.Properties.VariableNames)));
            tc.verifyTrue(all(X.peakFreqHz(isfinite(X.peakFreqHz)) >= 8 & X.peakFreqHz(isfinite(X.peakFreqHz)) <= 16));
            tc.verifySize(TF.on_window, [45 201]);  tc.verifyEqual(TF.n_window, n);
            tc.verifyEqual(unique(Acf.series)', sort(["env_real" "env_surrogate" "train_fit" "train_poisson_median" ...
                "train_poisson_p95" "train_surrogate" "train_window"]));
            tc.verifyEqual(Acf.acf(Acf.series == "env_real" & Acf.lag_s == 0), 1, 'AbsTol', 1e-12);
            tc.verifyTrue(any(T.band == "d3/75fps_h1"));
        end

        function plantsAreReadBackPerSnrBoundAndModel(tc)
            S = rheome.scale.sensors(tc.Name);
            [T, X] = rheome.scale.measure_eventplants(tc.Name, S, Depth=3, NPer=2, SNRs=[10 30], Hemis="L");
            tc.verifyEqual(height(X), 2 * 2 * 2 * 3);                                   % plants x SNRs x bounds x models
            tc.verifyTrue(any(T.band == "30dB/tau1/quadrature" & T.metric == "participant_pass"));
            q = X(X.snr_db == 30 & X.model == "quadrature" & X.fit, :);
            tc.verifyLessThan(median(abs(q.on_err_s)), 0.165);                          % clean plants: within a tile
            q = X(X.snr_db == 30 & X.model == "posterior" & X.fit, :);
            tc.verifyLessThan(median(abs(q.on_err_s)), 0.165);
        end

        function theDriverAndTheGroupReduce(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="events", OutDir=od, Dataset="synth");
            tc.verifyEqual(R.timing.status, "ok");
            for f = ["events.csv" "eventacf.csv" "eventtf.mat" "metrics.csv"]
                tc.verifyTrue(isfile(fullfile(od, tc.Name, f)), f);
            end
            [G, ~, TF] = rheome.scale.groupevents(od, fullfile(od, 'group'), Boot=20, BootTF=5);
            tc.verifyTrue(all(G.cohortKey == "synth"));  tc.verifyEqual(max(G.n), 1);
            tc.verifyEqual(TF.synth.n, 1);
            tc.verifyTrue(isfile(fullfile(od, 'group', 'group.csv')));
        end
    end
end
