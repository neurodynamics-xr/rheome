classdef tScaleAperiodic < matlab.unittest.TestCase
% The MS1 aperiodic test (nsp cf-aperiodic: G10) end to end on the synthetic cached subject: per-cell rows for
% both arms, rates and null percentiles in range, the split's reduction in [0, 1], the control at 0 dB in band,
% a threshold and a selectivity per arm, and the driver writing aperiodic.csv. Small settings: the numbers are
% checked for sense, not value.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_aperiodic'
        S
        small = {'SigmaMM', 30, 'OffsetX', [10 100], 'SpeedsMS', [0 0.5], 'FitWindowsS', [6 Inf], 'Depths', [2 3], ...
                 'BandMM', 130, 'Reps', 2, 'NullWindows', 4, 'PathMM', 60}
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

        function cellsThresholdsAndSelectivity(tc)
            [T, X] = rheome.scale.measure_aperiodic(tc.Name, tc.S, tc.small{:});
            % emptyroom: (2 offset + exponent + control) x 2 speeds x 3 signals; rest: no exponent
            tc.verifyEqual(height(X), 4*2*3 + 3*2*3);
            tc.verifyEqual(sort(unique(X.signal))', ["p6s" "pwhole" "total"]);
            tc.verifyTrue(all(X.reps == 2));
            for c = ["onRate_d2" "onRate_d3" "passRate_d3" "null95_d3"]
                tc.verifyTrue(all(X.(c) >= 0 & X.(c) <= 1/2), c);   % per second, a 2 s window
            end
            tc.verifyLessThanOrEqual(X.onRate_d3, X.passRate_d3);
            sp = X.signal ~= "total";
            tc.verifyTrue(all(X.reduction(sp) >= -1e-9 & X.reduction(sp) <= 1 + 1e-9));
            tc.verifyEqual(X.reduction(~sp), zeros(nnz(~sp), 1));
            tc.verifyEqual(X.plantDB(X.mode == "control"), zeros(nnz(X.mode == "control"), 1), 'AbsTol', 0.5);
            o = X(X.mode == "offset" & X.signal == "total" & X.arm == "rest" & X.vMS == 0, :);
            tc.verifyLessThan(o.plantDB(o.level == 10), o.plantDB(o.level == 100));   % the ladder is monotone
            tc.verifyTrue(all(isfinite(X.phiErrMM_130) & isfinite(X.psiErrNullMM_130)));
            th = T(T.metric == "threshold_db", :);
            tc.verifyEqual(height(th), 2 * 1 * 2 * 3 * 2);                 % arm x sigma x speed x signal x depth
            se = T(T.metric == "selectivity", :);
            tc.verifyEqual(sort(se.band)', ["emptyroom_p6s" "emptyroom_pwhole" "rest_p6s" "rest_pwhole"]);
        end

        function theDriverWritesTheTable(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.verifyTrue(ismember("aperiodic", rheome.scale.analyses()));
            tc.verifyFalse(ismember("aperiodic", rheome.scale.analyses("ported")));
            R = rheome.scale.run(tc.Name, Analyses="aperiodic", OutDir=od);
            % the defaults are the production grid; on the synthetic subject they must still run or fail loudly
            tc.verifyTrue(R.timing.status == "ok" || contains(R.timing.message, "scale:aperiodic"), ...
                          strjoin(R.timing.message, ' | '));
            if R.timing.status == "ok", tc.verifyTrue(isfile(fullfile(od, tc.Name, 'aperiodic.csv'))); end
        end
    end
end

% Author: Diellor Basha, 2026
