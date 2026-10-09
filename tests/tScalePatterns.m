classdef tScalePatterns < matlab.unittest.TestCase
% The MS1 pattern measures (nsp cf-patterns: G5, G13) end to end on the synthetic cached subject: the
% evaluated surrogate equals the rebuilt one, table shapes and halves, p-values in (0, 1], and the driver
% writing both tables. Small settings: the numbers are checked for sense, not value.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_patterns'
        S
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
            tc.S = rheome.scale.sensors(tc.Name);
        end
    end

    methods (Test)

        function evaluatedSurrogateIsThePhaseRandomisedRecord(tc)
            % the shortcut measure_patterns relies on: rotating the band's bins and evaluating at t equals
            % band-passing the phase-randomised record and taking its analytic signal
            rng(1);  fs = 100;  N = 1000;  Y = randn(3, N);  band = [8 12];
            k = 1:floor((N-1)/2);  f = k*fs/N;  m = f >= band(1) & f <= band(2);  r = exp(2i*pi*rand(1, nnz(m)));
            Fy = fft(Y, [], 2);  R = ones(1, N);  R(k(m) + 1) = r;  R(N - k(m) + 1) = conj(r);
            Yp = real(ifft(Fy .* R, [], 2));                          % the record, phase-randomised in band
            Fp = fft(Yp, [], 2);  H = zeros(1, N);  H(k(m) + 1) = 2;  Zp = ifft(Fp .* H, [], 2);
            t = [0.37 2.5 7.01];  n = round(t*fs);  t = n/fs;
            Ze = (2*Fy(:, k(m) + 1)/N .* r) * exp(2i*pi*f(m)'*t);
            tc.verifyEqual(Ze, Zp(:, n + 1), 'AbsTol', 1e-12);
        end

        function patternsHaveDetectorsHalvesAndPValues(tc)
            [T, X] = rheome.scale.measure_patterns(tc.Name, tc.S, "patterns", Hemis="L", NSurr=4, FrameS=2, ...
                                                   SurrStrideS=8, PacSurrWindows=2, PacWindowS=4);
            dets = ["planarity" "source" "sink" "saddle" "vortex" "standing" "speed" "singularities"];
            nb = numel(unique(X.band(X.detector ~= "pac")));
            tc.verifyEqual(height(X), (numel(dets)*nb + 1) * 3);      % (detectors x bands + pac) x halves
            tc.verifyTrue(all(ismember([dets "pac"], X.detector)));
            tc.verifyTrue(all(X.p > 0 & X.p <= 1));
            tc.verifyEqual(X.e, X.f - 0.05, 'AbsTol', 1e-12);
            tc.verifyTrue(all(isfinite(X.fER(X.half == 0 & X.detector == "planarity"))));   % the empty-room arm ran
            tc.verifyTrue(all(X.nSurr == 4));
            tc.verifyTrue(any(T.metric == "iaf_hz"));
            tc.verifyTrue(any(T.metric == "excess_e" & T.band == "L_alpha_vortex"));
        end

        function nullsHaveTheThreeStatistics(tc)
            [T, X] = rheome.scale.measure_patterns(tc.Name, tc.S, "patternnulls", Hemis="L", NumWindows=4, ...
                                                   NumTiles=3, NPerm=10);
            tc.verifyEqual(sort(unique(X.statistic))', ["dispersion" "peakedness" "rotation_r"]);
            tc.verifyEqual(height(X), 3*3);
            tc.verifyTrue(all(X.p > 0 & X.p <= 1));
            pk = X(X.statistic == "peakedness" & X.half == 0, :);
            tc.verifyGreaterThanOrEqual(pk.value, 1);                 % max over mean
            tc.verifyTrue(any(T.metric == "z" & T.band == "L_rotation_r"));
        end

        function theDriverWritesBothTables(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses="patternnulls", OutDir=od);
            tc.verifyEqual(R.timing.status, "ok", strjoin(R.timing.message, ' | '));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'patternnulls.csv')));
        end
    end
end

% Author: Diellor Basha, 2026
