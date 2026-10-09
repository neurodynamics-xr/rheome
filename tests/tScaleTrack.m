classdef tScaleTrack < matlab.unittest.TestCase
% The MS1 tracking measures (nsp cf-track: G7 injections, G8 the 0.52 factorial) end to end on the
% synthetic cached subject, at small settings: table shapes, balanced halves, every factorial cell present
% with both head levels, and the driver writing both tables. The numbers are checked for sense, not value.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_track'
        Ref = 'synth_trackref'
        S
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));  scaleSynthSubject(fullfile(d, tc.Ref));
            tc.S = rheome.scale.sensors(tc.Name);
        end
    end

    methods (Test)

        function injectionsAreScoredWithTheirControls(tc)
            [T, X] = rheome.scale.measure_inject(tc.Name, tc.S, Depths=[2 3], JLevels=4, NPer=2, Hemis="L", ...
                                                 Cases=["still" "move05"], Speeds=[0 0.05]);
            q = X(X.case ~= "real", :);
            tc.verifyEqual(height(q), 2*2*2);                       % cases x injections x depths
            tc.verifyEqual(sort(q.half(q.depth == 2 & q.case == "still"))', [1 2]);
            tc.verifyTrue(all(ismember(q.found, [0 1])) && all(ismember(q.c_found, [0 1])));
            r = X(X.case == "real", :);
            tc.verifyEqual(height(r), 2*floor(40/2));               % every 2 s window of the 40 s record, per depth
            tc.verifyTrue(all(isfinite(X.nullFalseS)));
            tc.verifyTrue(all(ismember(["hit_rate" "speed_ratio" "false_rate" "real_frac_pass"], T.metric)));
            tc.verifyTrue(any(T.band == "move05/d3/18.75fps_h2"));
        end

        function factorialHasEveryCellAndBothHeads(tc)
            [T, X] = rheome.scale.measure_trackfactorial(tc.Name, tc.S, RefHead=tc.Ref, Speeds=[0.05 0.5], ...
                                                         Depths=[2 3], NRep=1, NNull=4, Hemis="L");
            I = X(X.arm == "instrument", :);
            tc.verifyEqual(height(unique(I(:, {'ruler','signal','background','matching','head'}))), 32);
            tc.verifyEqual(height(I), 32 * 2*1*2);                   % cells x speeds x reps x depths
            D = X(X.arm == "direct", :);
            tc.verifyEqual(height(D), 4 * 2*1*2);                    % (ruler x matching) x speeds x reps x depths
            tc.verifyTrue(all(isnan(X.speedRatio(X.hit == 0))));     % speed only for hits
            tc.verifyEqual(T.value(T.metric == "n_head_levels"), 2);
            tc.verifyEqual(sum(T.metric == "hit_rate" & ~contains(T.band, "/d") & ~startsWith(T.band, "direct")), 32);
            [~, X1] = rheome.scale.measure_trackfactorial(tc.Name, tc.S, Speeds=0.5, Depths=3, NRep=1, NNull=4, Hemis="L");
            tc.verifyEqual(unique(X1.head), "own");                  % no RefHead: the own level alone
        end

        function theReferenceHeadRoundTrips(tc)
            f = fullfile(tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder, 'ref.tar.gz');
            rheome.scale.refhead(tc.Ref, f);
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            untar(f, fullfile(d, 'unpacked'));
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  c = onCleanup(@() setenv('RHEOME_DATA', old));
            Sr = rheome.scale.sensors('unpacked');
            tc.verifyEqual(Sr.G, tc.S.G);                             % the same synthetic head
            st = rheome.load.study('unpacked');  tc.verifyEmpty(st.rec.F);
        end

        function theDriverWritesTheFactorial(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            t0 = tic;  R = rheome.scale.run(tc.Name, Analyses="trackfactorial", OutDir=od);  fprintf('driver %.0f s\n', toc(t0));
            tc.verifyEqual(R.timing.status, "ok", strjoin(R.timing.message, ' | '));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'trackfactorial.csv')));
        end
    end
end

% Author: Diellor Basha, 2026
