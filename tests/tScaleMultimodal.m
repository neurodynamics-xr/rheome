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
            % relative band power: the band's share of the four bands' stored sums, summing to 1 per tile and member
            R = Y(endsWith(Y.field, "_rel"), :);
            tc.verifyEqual(height(R), 4 * height(Y(Y.field == "meg_alpha", :)));
            [g, ~] = findgroups(R(:, {'level','code','member'}));
            tc.verifyEqual(splitapply(@sum, R.mean, g), ones(max(g), 1), 'AbsTol', 1e-12);
            a = Y(Y.level == 5 & Y.member == 2 & Y.field == "meg_alpha", :);
            M = Y(Y.level == 5 & Y.member == 2 & Y.modality == "MEG" & ~endsWith(Y.field, "_rel"), :);
            ra = R(R.level == 5 & R.member == 2 & R.field == "meg_alpha_rel", :);
            tot = splitapply(@sum, M.sum_ax, findgroups(M.code));
            tc.verifyEqual(ra.mean, a.sum_ax ./ tot, 'RelTol', 1e-12);
            tc.verifyTrue(any(X.assoc.field_a == "meg_alpha_rel" | X.assoc.field_b == "meg_alpha_rel"));
        end

        function boldInfraslow(tc)
            % a planted 0.05 Hz, 1 % BOLD oscillation on every vertex, plus a 0.2 Hz one (out of band) and a slow
            % drift that a confound carries: the infra-slow power is the planted A^2/2 = 0.5 %^2, the drift is
            % regressed out, and a spatially constant signal has no power in the wavelet members (lambda = 0)
            nm = 'synth_multimodal_bold';  src = fullfile(tc.Root, tc.Name);  d = fullfile(tc.Root, nm);
            copyfile(src, d);  delete(fullfile(d, 'fibers_subject.mat'));
            sf = load(fullfile(d, 'surface.mat'));  nV = size(sf.S.Vertices, 1);
            tr = 2;  t = (0:299) * tr;  drift = sin(2*pi*0.03*t);
            conf = [drift; 0.01 * randn(7, numel(t), 'like', t)]';
            x = 1000 + 10 * sin(2*pi*0.05*t) + 10 * sin(2*pi*0.2*t) + 50 * drift;
            x = repmat(x, nV, 1);  x(1:2, :) = 0;              % two vertices outside the brain mask
            names = ["trans_x" "trans_y" "trans_z" "rot_x" "rot_y" "rot_z" "csf" "white_matter"];
            fmri = struct('name', 'task_rest_run_1', 'x', single(x), 'tr', tr, 'units', 'a.u.', ...
                          'confounds', conf, 'confoundNames', names, 'file', '');
            save(fullfile(d, 'fmri.mat'), 'fmri');
            [T, X] = rheome.scale.measure_multimodal(nm, [], Depth=4);
            v = @(m) T.value(T.metric == m);
            tc.verifyEqual(v("n_fields"), 6);                  % 4 MEG bands, fMRI, PET
            tc.verifyEqual([v("n_bold_runs") v("bold_seconds") v("bold_lowmean_vertices")], [1 600 2]);
            Y = X.tiles(X.tiles.field == "fmri_infraslow", :);
            tc.verifyEqual(unique(Y.modality), "fMRI");  tc.verifyEqual(unique(Y.unit), "%^2");
            r0 = Y(Y.level == 0 & Y.member == 0, :);
            tc.verifyEqual(r0.mean, 0.5, 'RelTol', 0.05);
            tc.verifyLessThan(max(Y.mean(Y.member > 0)), 1e-3 * r0.mean);
            tc.verifyLessThan(v("rollup_max_relerr"), 1e-10);
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
