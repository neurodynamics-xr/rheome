classdef tScaleMultimodal < matlab.unittest.TestCase
% The MS1 multimodal example (rheome.scale.measure_multimodal) on the synthetic cached subject
% (scaleSynthSubject) with a PET map and homotopic fibres: every level is an exact roll-up of the finest
% tiles, the connectome counts every fibre end once, and rheome.scale.run writes the three outputs.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_multimodal'
        Root
        nh
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            d = fullfile(tc.Root, tc.Name);  scaleSynthSubject(d);
            sf = load(fullfile(d, 'surface.mat'));  S = sf.S;  tc.nh = numel(S.Hemi{1});
            k = (1:4:tc.nh)';  Points = zeros(numel(k), 5, 3);       % homotopic: vertex i of L to vertex i of R
            Points(:, 1, :) = S.Vertices(k, :);  Points(:, 5, :) = S.Vertices(k + tc.nh, :);
            Points(1, 5, :) = S.Vertices(k(1), :);                   % one fibre ending where it starts
            save(fullfile(d, 'fibers_subject.mat'), 'Points', '-v7.3');
            pet = struct('name', 'tracer', 'x', 1 + S.VertNormals(:, 3), 'units', 'SUVR', 'file', '');
            save(fullfile(d, 'pet.mat'), 'pet');
        end
    end

    methods (Test)
        function rollupIsExact(tc)
            [T, X, Cn] = rheome.scale.measure_multimodal(tc.Name, [], Depth=4);
            v = @(m) T.value(T.metric == m);
            tc.verifyLessThan(v("rollup_max_relerr"), 1e-10);
            tc.verifyLessThan(v("connectome_rollup_max_relerr"), 1e-12);
            tc.verifyEqual(v("n_tiles_finest"), 32);  tc.verifyEqual(v("n_levels"), 6);
            tc.verifyEqual(v("n_fields"), 6);              % 4 MEG bands, 1 PET, fibre degree
            % level 0 is the whole cortex; its area is the sum of the finest tiles'
            Y = X.tiles;  r0 = Y(Y.level == 0 & Y.member == 0 & Y.field == "pet_tracer", :);
            r5 = Y(Y.level == 5 & Y.member == 0 & Y.field == "pet_tracer", :);
            tc.verifyEqual(r0.area_m2, sum(r5.area_m2), 'RelTol', 1e-12);
            tc.verifyEqual(r0.mean, sum(r5.sum_ax) / sum(r5.area_m2), 'RelTol', 1e-12);
            tc.verifyEqual(r0.code, 1);  tc.verifyEqual(sort(unique(Y.code(Y.level == 1)))', [2 3]);
            % one fibre is a self-loop; the other ends are counted twice (symmetric), all interhemispheric
            nf = numel(1:4:tc.nh);
            tc.verifyEqual(v("connectome_total"), 2 * (nf - 1));
            tc.verifyEqual(Cn.C{2}, [0 nf-1; nf-1 0]);
            tc.verifyEqual(v("interhemispheric_fraction"), 1);
            tc.verifyEqual(numel(Cn.C), 6);
            % the descriptive rho table covers every pair of fields at levels 3..5
            tc.verifyEqual(unique(X.assoc.level)', 3:5);
            tc.verifyTrue(all(abs(X.assoc.rho(~isnan(X.assoc.rho))) <= 1 + 1e-12));
        end

        function viaRun(tc)
            od = fullfile(tc.Root, 'out');
            R = rheome.scale.run(tc.Name, Analyses="multimodal", OutDir=od);
            tc.verifyEqual(R.timing.status, "ok", R.timing.message);
            for f = ["multimodal.csv" "multimodal_assoc.csv" "multimodal_connectome.mat" "metrics.csv"]
                tc.verifyTrue(isfile(fullfile(od, tc.Name, f)), f);
            end
        end
    end
end

% Author: Diellor Basha, 2026
