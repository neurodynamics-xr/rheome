classdef tScaleCatalogue < matlab.unittest.TestCase
% The MS1 catalogue plants (nsp cf-plants: G6, G3 rule 11, G16): the two comparators on an analytic wave, and
% measure_catalogue end to end on the synthetic cached subject -- arms, estimators, ratios, strips, nulls, driver.
% Brainstorm's optical flow runs only when a Brainstorm checkout is found (RHEOME_BRAINSTORM, or ../brainstorm3
% next to this checkout); otherwise the bst_of rows are skipped (HornSchunck = []) and its own test is filtered.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_catalogue'
        S
        HS = []
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
            tc.S = rheome.scale.sensors(tc.Name);
            if strlength(tc.bst()) > 0, tc.HS = 0.01; end
        end
    end

    methods (Test)

        function phaseRegressionReadsAPlaneWave(tc)
            [S, z, c, fs] = tc.sphereWave(@(V, t) exp(1i*(2*pi*10*t - 2*pi/0.06*V(:,1))));
            out = rheome.flow.phaseregression(z, S, Rate=fs, Centres=c);
            tc.verifyEqual(median(out.speed), 0.6, 'RelTol', 0.05);                         % 60 mm at 10 Hz
            tc.verifyLessThan(median(acosd(out.velocity(:,1) ./ vecnorm(out.velocity, 2, 2))), 5);   % towards +x
            tc.verifyGreaterThan(median(out.R), 0.95);
        end

        function phaseRegressionFindsNoGradientInAStandingWave(tc)
            [S, z, c, fs] = tc.sphereWave(@(V, t) cos(2*pi/0.06*V(:,1)) * exp(1i*2*pi*10*t));
            out = rheome.flow.phaseregression(z, S, Rate=fs, Centres=c);
            tc.verifyGreaterThan(mean(isinf(out.speed)), 0.5);
        end

        function brainstormOpticalFlowFollowsTheWave(tc)
            tc.assumeNotEmpty(tc.HS, 'no Brainstorm checkout');
            [S, z, c, fs] = tc.sphereWave(@(V, t) exp(1i*(2*pi*10*t - 2*pi/0.06*V(:,1))));
            of = rheome.flow.bstopticalflow(real(z(:, 1:17)), S, fs, BrainstormDir=tc.bst());
            v = mean(of.velocity(c, :, :), 3);
            tc.verifyLessThan(acosd(mean(v(:,1) ./ vecnorm(v, 2, 2))), 10);
            tc.verifyEqual(size(of.velocity), [size(z, 1) 3 16]);
        end

        function catalogueHasEveryArmAndEstimator(tc)
            [T, X, strips] = rheome.scale.measure_catalogue(tc.Name, tc.S, "catalogue", Hemis="L", Centres=1, ...
                Patterns=["planar_macro" "rotor_meso" "source_meso"], HornSchunck=tc.HS, PatchMM=15);
            arm = X.arm + "/" + X.noise;
            tc.verifyEqual(sort(unique(arm))', ["direct/none" "meg/emptyroom" "meg/inf" "meg/rest" "meglocal/inf"]);
            tc.verifyEqual(unique(arm(X.arm == "meglocal")), "meglocal/inf");
            tc.verifyEqual(unique(X.pattern(X.arm == "meglocal")), "planar_macro");          % G3 rule 11 only
            est = ["framework" "phasereg"];  if ~isempty(tc.HS), est = ["bst_of" est]; end
            tc.verifyEqual(sort(unique(X.estimator))', est);
            tc.verifyFalse(any(X.estimator ~= "framework" & X.noise == "emptyroom"));      % comparators: 3 arms
            d = X(X.arm == "direct" & X.estimator == "framework", :);
            tc.verifyEqual(d.recovered, d.truthOnMesh);                                    % the direct arm IS the truth
            sp = X(X.quantity == "phase speed" & X.pattern == "planar_macro", :);
            tc.verifyTrue(all(isfinite(sp.ratio(sp.estimator == "framework"))));
            tc.verifyTrue(all(ismember(X.recovered(X.quantity == "propagation declared"), [0 1])));
            tc.verifyTrue(all(T.analysis == "catalogue"));
            tc.verifyTrue(any(T.metric == "preserved_frac") && any(T.metric == "ratio_median"));
            tc.verifyEqual(sort(fieldnames(strips))', {'planar_macro' 'rotor_meso' 'source_meso'});
            tc.verifyEqual(sort(fieldnames(strips.rotor_meso))', {'direct_none' 'meg_inf' 'meg_rest'});
            tc.verifySize(strips.rotor_meso.meg_inf.frames, [size(tc.S.B.L.S.Vertices, 1) 6]);
        end

        function nullsCountFalsePropagation(tc)
            [T, X] = rheome.scale.measure_catalogue(tc.Name, tc.S, "catalognulls", Hemis="L", NullsPerHemi=1, ...
                HornSchunck=tc.HS, PatchMM=15, NullSepMM=[20 40]);
            nE = 2 + numel(tc.HS);
            tc.verifyEqual(height(X), 2 * 3 * nE);                                        % types x arms x estimators
            tc.verifyTrue(all(ismember(X.declared, [0 1])));
            tc.verifyTrue(all(X.coherence >= 0 & X.coherence <= 1 + 1e-12));
            tc.verifyEqual(height(T), height(X));                                         % one rate per cell at n = 1
            tc.verifyTrue(all(T.metric == "false_prop_rate"));
        end

        function theDriverWritesTheNullTable(tc)
            if isempty(tc.HS), return, end                                                 % the default sweep needs Brainstorm
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_BRAINSTORM');  setenv('RHEOME_BRAINSTORM', tc.bst());  tc.addTeardown(@setenv, 'RHEOME_BRAINSTORM', old);
            R = rheome.scale.run(tc.Name, Analyses="catalognulls", OutDir=od);
            tc.verifyEqual(R.timing.status', "ok", strjoin(R.timing.message, ' | '));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'catalognulls.csv')));
        end
    end

    methods
        function d = bst(~)
            d = string(getenv('RHEOME_BRAINSTORM'));
            if strlength(d) == 0
                d = string(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'brainstorm3'));
            end
            if exist('bst_opticalflow', 'file') ~= 2 && ~isfolder(fullfile(d, 'toolbox', 'math')), d = ""; end
        end

        function [S, z, c, fs] = sphereWave(~, f)
            [V, F] = rheome.geom.icosphere(5);  V = 0.1 * V;
            S = struct('Vertices', V, 'Faces', F, 'VertNormals', V ./ vecnorm(V, 2, 2), 'nV', size(V, 1));
            fs = 160;  t = (0:31) / fs;  z = f(V, t);
            [~, o] = sort(V(:, 3), 'descend');  c = o(1:20);                              % near the pole: x is a chart
        end
    end
end

% Author: Diellor Basha, 2026
