classdef tScaleFieldEvents < matlab.unittest.TestCase
% The two MS1 measures end to end on a synthetic cached subject (two icosphere hemispheres, a random
% leadfield, 40 s of noise with a 10 Hz source): fieldsmooth finds the band-limited field smoother than
% the per-vertex one, and eventsensors maps the best Viterbi path onto sensor samples and channels.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_fieldevents'
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
        end
    end

    methods (Test)

        function bandLimitedFieldIsSmoother(tc)
            S = rheome.scale.sensors(tc.Name);
            [T, X, maps] = rheome.scale.measure_fieldsmooth(tc.Name, S, CutoffsMM=[60 120], NumFrames=4);
            v = @(m, b) T.value(T.metric == m & T.band == b);
            tc.verifyGreaterThan(v("J_wavelength_mm", "bl120"), v("Jt_wavelength_mm", "raw"));
            tc.verifyGreaterThan(v("div_wavelength_mm", "bl120"), v("div_wavelength_mm", "raw"));
            tc.verifyGreaterThan(v("J_coherence", "bl120"), v("Jt_coherence", "raw"));
            tc.verifyGreaterThanOrEqual(v("J_wavelength_mm", "bl120"), v("J_wavelength_mm", "bl60"));
            tc.verifyEqual(height(X), 4 * (7 + 2*5));
            tc.verifySize(maps.bl, [1 2]);  tc.verifySize(maps.raw.Jt, [3*size(maps.surface.Vertices, 1) 1]);
        end

        function eventsLandOnTheirFramesSamples(tc)
            S = rheome.scale.sensors(tc.Name);
            [~, ~, Best] = rheome.scale.measure_grouptrack(tc.Name, S, Depths=3, JLevels=2, Hemis="L");
            tc.assertNotEmpty(Best);
            [T, ev] = rheome.scale.measure_eventsensors(tc.Name, S, Best, Depth=3, JLevel=2);
            tc.verifyEqual(T.value(T.metric == "n_steps"), numel(Best(1).track.frames));
            E = ev.E;  fs = 600;  dec = 2;  bs = 4;
            f1 = (Best(1).window - 1) * 600 + (Best(1).track.frames(1) - 1) * bs;      % 300-fps frames
            s1 = f1 * dec + 1 - round(ev.t(1) * fs);                                      % sample in the crop
            tc.verifyEqual(E.samples(1, :), [s1, s1 + bs*dec - 1]);
            tc.verifyTrue(any(E.mask(:, s1)));  tc.verifyEqual(size(E.mask), size(ev.Xraw));
            tc.verifySize(ev.pos2, [size(ev.Xraw, 1) 2]);
            h = rheome.show.eventsensors(ev.Xband, ev.t, E, ev.pos2, 'Labels', ev.labels, 'Visible', 'off');
            close(h);
        end

        function theDriverWritesTablesAndFigureData(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses=["fieldsmooth" "eventsensors"], OutDir=od);
            tc.verifyEqual(R.timing.status', ["ok" "error"]);           % eventsensors without grouptrack
            tc.verifySubstring(R.timing.message(2), "needs grouptrack first");
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'fieldsmooth.csv')));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'fieldsmooth_maps.mat')));
            tc.verifyTrue(any(R.metrics.metric == "div_wavelength_mm" & R.metrics.band == "bl92"));
        end

        function theDriverWritesPrognomeCoefficients(tc)
            % rheome.scale.run(..., Analyses="coefficients") -> rheome_coeffs.mat, Prognome's contract
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="coefficients", OutDir=od);
            tc.assertEqual(R.timing.status, "ok", R.timing.message);
            C = load(fullfile(od, tc.Name, 'rheome_coeffs.mat'));
            tc.verifyEmpty(setdiff(["W" "codes" "scales" "fs"], string(fieldnames(C))));
            tc.verifyEqual(C.fs, 100);
            tc.verifyEqual(C.codes, (2^8:2^9-1)');                    % depth 7 per hemisphere: level 8
            tc.verifySize(C.W, [256 numel(C.scales) 40*100]);          % 40 s at 100 Hz
            tc.verifyTrue(all(isfinite(C.W), 'all') && all(C.W >= 0, 'all'));
            tc.verifyEqual(R.metrics.value(R.metrics.metric == "n_tiles"), 256);
        end
    end
end

% Author: Diellor Basha, 2026
