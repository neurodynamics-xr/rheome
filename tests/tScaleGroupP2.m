classdef tScaleGroupP2 < matlab.unittest.TestCase
% The MS1 group P2 measures on a synthetic cached subject (scaleSynthSubject): Table 5 per band and per
% half, Helmholtz-band recovery, the own-region fraction, the atlas geometry checks and the fused-kernel
% exactness, each alone and through rheome.scale.run.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_groupp2'
        Root
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(tc.Root, tc.Name));
            % a 1/f^2 background on every channel and a 10 Hz rhythm on top: alpha is the only rhythm
            f = fullfile(tc.Root, tc.Name, 'study.mat');  st = load(f);  rng(3);
            [C, N] = size(st.rec.F);  t = (0:N-1) / st.rec.sfreq;
            bg = filter(1, [1 -0.995], randn(N, C)).';  bg = bg - mean(bg, 2);
            st.rec.F = bg + 3 * std(bg(:)) * (0.5 + rand(C, 1)) .* sin(2*pi*10*t);
            save(f, '-struct', 'st');
            nz = load(fullfile(tc.Root, tc.Name, 'noise.mat'));  nz.nrec.F = 1e-2 * std(bg(:)) * randn(C, 20 * st.rec.sfreq);
            save(fullfile(tc.Root, tc.Name, 'noise.mat'), '-struct', 'nz');
        end
    end

    methods (Test)

        function alphaIsTheOnlyRhythmInEveryHalf(tc)
            [T, X] = rheome.scale.measure_bandperiodic(tc.Name, [], NumSeeds=20);
            tc.verifyEqual(unique(X.analysis)', ["bandperiodic" "bandperiodic_h1" "bandperiodic_h2"]);
            for a = unique(X.analysis)'
                x = X(X.analysis == a, :);  ia = x.band == "8-16 Hz";
                tc.verifyGreaterThan(x.osc_snr_dB(ia), 0, a);
                tc.verifyLessThan(max(x.osc_snr_dB(~ia)), x.osc_snr_dB(ia) - 10, a);
                tc.verifyGreaterThan(x.periodic_fraction(ia), 0.5, a);
                tc.verifyEqual(T.value(T.analysis == a & T.metric == "iaf_hz"), 10, 'AbsTol', 0.25);
            end
            v = @(m, b) T.value(T.analysis == "bandperiodic" & T.metric == m & T.band == b);
            tc.verifyGreaterThan(v("fit_r2", ""), 0.8);
            tc.verifyGreaterThan(v("total_snr_dB", "8-16 Hz"), v("osc_snr_dB", "8-16 Hz"));   % the background is not rhythm
            tc.verifyTrue(all(isfinite(X.floor_total_mm)) && all(X.floor_total_mm > 0));
            tc.verifyEqual(X.clamped_osc, double(X.SnrFixed_osc < 0.2));
        end

        function plantedHelmholtzBandsAreRecovered(tc)
            [T, X] = rheome.scale.measure_helmholtzbands(tc.Name, Plants=2);
            tc.verifyEqual(unique(X.hemi)', ["L" "R"]);
            for h = ["L" "R"]
                tc.verifyGreaterThanOrEqual(T.value(T.metric == "recovered_all_" + h), 0.8);
            end
            tc.verifyLessThan(max([X.leakSource_pct; X.leakVortex_pct]), 1);
        end

        function geometryChecksHoldOnTheCortex(tc)
            [T, X] = rheome.scale.measure_geometry(tc.Name, [], Levels=2:3, RollupDepths=2:4);
            v = @(m, b) T.value(T.metric == m & T.band == b);
            for h = ["L" "R"]
                tc.verifyEqual(v("n_singular", h), 2);  tc.verifyEqual(v("charge_sum", h), 2);
                tc.verifyLessThan(v("rollup_max_rel_error", h), 1e-12);
                tc.verifyLessThan(v("residual_rad", h), 1e-9);
            end
            tc.verifyEqual(height(X), 2 * (4 + 8));
            tc.verifyGreaterThan(v("lobe_ratio", "L3"), 0.5);  tc.verifyLessThan(v("lobe_ratio", "L3"), 1.5);
            tc.verifyGreaterThan(v("energy_in_tile", "L3"), 0.5);
        end

        function parcelAndTileSizeOnOneRuler(tc)
            [T, ~, Z] = rheome.scale.measure_geometry(tc.Name, [], Levels=2, RollupDepths=2:3);
            v = @(m, b) T.value(T.metric == m & T.band == b);
            P = Z.parcels;  t = Z.tiles;
            tc.verifyEqual(P.atlas', repmat("Desikan-Killiany", 1, 4));   % the synthetic subject has no Destrieux
            tc.verifyEqual(v("n_parcels", "Destrieux"), 0);
            tc.verifyEqual(P.hemi', ["L" "L" "R" "R"]);
            tc.verifyEqual(P.d_eq_mm, 2 * sqrt(P.area_mm2 / pi), 'RelTol', 1e-12);
            tc.verifyEqual(height(t), 2 * (1 + 2 + 4 + 8));
            for h = ["L" "R"]                  % depth 0 is the hemisphere; each depth partitions its area
                A0 = t.area_mm2(t.hemi == h & t.depth == 0);
                tc.verifyEqual(v("hemisphere_area", h), A0 / 100, 'RelTol', 1e-12);
                tc.verifyEqual(sum(P.area_mm2(P.hemi == h)), A0, 'RelTol', 1e-12);
                for d = 1:3, tc.verifyEqual(sum(t.area_mm2(t.hemi == h & t.depth == d)), A0, 'RelTol', 1e-12); end
            end
            tc.verifyEqual(v("hemisphere_area", "L"), 4 * pi * 4^2, 'RelTol', 0.02);   % icosphere of radius 40 mm
            tc.verifyEqual(v("tile_deq_median", "D0"), 2 * sqrt(4 * pi * 40^2 / pi), 'RelTol', 0.02);
            tc.verifyGreaterThan(t.geodesic_diameter_mm(t.depth == 0), 0.95 * pi * 40);   % half a great circle
            tc.verifyLessThan(t.geodesic_diameter_mm(t.depth == 0), 1.15 * pi * 40);
        end

        function ownRegionAndFusion(tc)
            ctx = rheome.flow.context(tc.Name);  K = rheome.flow.build(ctx);
            [T, X] = rheome.scale.measure_ownregion(tc.Name, [], K);
            tc.verifyEqual(height(X), 4);
            tc.verifyTrue(all(X.own_curl >= 0 & X.own_curl <= 1 & X.own_div >= 0 & X.own_div <= 1));
            tc.verifyEqual(T.value(T.metric == "frac_regions_pass" & T.band == "curl"), mean(X.own_curl >= 0.25));
            [F, Y] = rheome.scale.measure_fusion(tc.Name, ctx, K, NumFrames=3);
            tc.verifyEqual(height(Y), 3);
            % div, curl and the Grams are exact; the potentials' check is NOT asserted here: on these
            % spheres Psi is ~1e-3 of Phi and the batch-vs-single Poisson solve differs at that level
            for q = ["divVertex" "curlVertex" "divCoeff" "curlCoeff" "energy" "enstrophy"]
                tc.verifyLessThan(F.value(F.metric == q + "_worst"), 1e-10, q);
            end
        end

        function theDriverRunsEveryP2Measure(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            a = ["bandperiodic" "helmholtzbands" "ownregion" "geometry" "fusion"];
            R = rheome.scale.run(tc.Name, Analyses=a, OutDir=od);
            tc.verifyEqual(R.timing.status', repmat("ok", 1, 5), strjoin(R.timing.message, ' | '));
            for f = [a "geometry_parcels" "geometry_tiles"], tc.verifyTrue(isfile(fullfile(od, tc.Name, f + ".csv")), f); end
            tc.verifyEqual(unique(R.metrics.analysis)', sort(["bandperiodic" "bandperiodic_h1" "bandperiodic_h2" ...
                "helmholtzbands" "ownregion" "geometry" "fusion"]));
            tc.verifyTrue(all(ismember(a, rheome.scale.analyses())));
        end
    end
end

% Author: Diellor Basha, 2026
