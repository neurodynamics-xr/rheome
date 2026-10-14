classdef tScaleAmbient < matlab.unittest.TestCase
% MS1 1.4 and 2.3 panels: rheome.scale.measure_ambient (the 2H(J.n) coupling against its closed form) and
% rheome.scale.measure_cptrack (a drifting rotor and a waxing/waning source through rheome.detect.track), on
% the sphere and on a synthetic cached subject (scaleSynthSubject), alone and through rheome.scale.run.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_ambient'
        Root
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(tc.Root, tc.Name));
        end
    end

    methods (Test)

        function couplingMatchesItsClosedFormOnTheSphere(tc)
            [T, X] = rheome.scale.measure_ambient('');
            tc.verifyEqual(unique(X.surface), "sphere");
            tc.verifyLessThan(X.H_relerr(1), 0.01);                     % cotan H = 1/R
            k = cosd(X.tilt_deg) > 1e-9;
            tc.verifyLessThan(max(X.couple_relerr(k)), 0.03);           % C = 2H(J.n): 1.7% for a 20 mm bump at ico5
            tc.verifyGreaterThan(min(X.couple_r(k)), 0.999);
            tc.verifyLessThan(max(X.curl_normal_rel(k)), 0.01);         % the normal part makes no curl (0.3% at ico5)
            tc.verifyEqual(X.normal_share(X.tilt_deg == 90), 0, 'AbsTol', 1e-12);
            tc.verifyEqual(X.normal_share(X.tilt_deg == 0), 1, 'AbsTol', 1e-9);   % purely normal: all of it
            tc.verifyTrue(all(T.analysis == "ambient") && any(T.band == "sphere/45 deg"));
        end

        function cortexRowsAndFigureData(tc)
            [T, X, Ex] = rheome.scale.measure_ambient(tc.Name, TiltsDeg=[0 45]);
            tc.verifyEqual(X.surface', ["sphere" "sphere" "cortex" "cortex"]);
            tc.verifyTrue(all(isfinite(X.couple_relerr_60mm(X.surface == "cortex"))));
            tc.verifyTrue(isfield(Ex, 'sphere') && isfield(Ex, 'cortexL'));
            tc.verifyEqual(size(Ex.cortexL.J), [size(Ex.cortexL.vertices_mm, 1) 3]);
            tc.verifyTrue(any(T.metric == "couple_relerr_60mm"));
        end

        function rotorAndSourceAreTrackedOnTheSphere(tc)
            [T, X, Ex] = rheome.scale.measure_cptrack('');
            x = X(1, :);
            tc.verifyGreaterThan(x.rotor_cover, 0.9);
            tc.verifyLessThan(x.rotor_err_mm_median, 3);                % edge ~3.9 mm at ico5, R = 100 mm
            tc.verifyEqual(x.rotor_speednet_ms, x.speed_true_euclid_ms, 'RelTol', 0.15);
            tc.verifyLessThan(x.rotor_curl_relerr, 0.15);
            tc.verifyLessThan(x.rotor_div_leak, 0.1);
            tc.verifyEqual(x.source_on, x.source_on_true, 'AbsTol', 2);  % born and dies with its amplitude
            tc.verifyEqual(x.source_off, x.source_off_true, 'AbsTol', 2);
            tc.verifyLessThan(x.source_div_relerr, 0.15);
            tc.verifyLessThan(x.source_curl_leak, 0.1);
            tc.verifyTrue(all(T.analysis == "cptrack"));
            tc.verifyEqual(numel(Ex.sphere.t_s), 41);
        end

        function bothRunThroughTheDriver(tc)
            od = fullfile(tc.Root, 'out');
            R = rheome.scale.run(tc.Name, Analyses=["ambient" "cptrack"], OutDir=od);
            tc.verifyEqual(R.timing.status', ["ok" "ok"], strjoin(R.timing.message, ' | '));
            for f = ["ambient.csv" "ambient.mat" "cptrack.csv" "cptrack.mat" "metrics.csv"]
                tc.verifyTrue(isfile(fullfile(od, tc.Name, f)), f);
            end
            tc.verifyTrue(all(ismember(["ambient" "cptrack"], rheome.scale.analyses())));
            tc.verifyFalse(any(ismember(["ambient" "cptrack"], rheome.scale.analyses("ported"))));
        end
    end
end

% Author: Diellor Basha, 2026
