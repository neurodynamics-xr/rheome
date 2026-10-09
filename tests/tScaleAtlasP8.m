classdef tScaleAtlasP8 < matlab.unittest.TestCase
% MS1 group P8 (G12) on the synthetic cached subject (scaleSynthSubject), with the subject itself as the
% template: correspondence through the sphere is exact, raw coordinates are not once the template is
% shifted; atlas events and their surrogates; the shared-gauge tensor; Destrieux -> DK vs the dyadic
% roll-up and the connectome-wavelet gamma sweep on homotopic fibres; and the four through rheome.scale.run.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_atlasp8'
        Root
        Tmpl
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            d = fullfile(tc.Root, tc.Name);  scaleSynthSubject(d);
            sf = load(fullfile(d, 'surface.mat'));  S = sf.S;  nh = numel(S.Hemi{1});
            U = S.VertNormals;  S.Sphere = 0.1 * U;            % each hemisphere is a sphere: its own registration
            sf.S = S;  save(fullfile(d, 'surface.mat'), '-struct', 'sf');
            at = load(fullfile(d, 'atlas.mat'));  hL = (1:S.nV)' <= nh;  y = U(:,2) > 0;  z = U(:,3) > 0.3;
            at.atlas(2) = struct('Name', 'Destrieux', 'Label', {{'p L','q L','r L','p R','q R','r R'}}, ...
                'Membership', sparse([y & z & hL, ~y & z & hL, ~z & hL, y & z & ~hL, ~y & z & ~hL, ~z & ~hL]'));
            save(fullfile(d, 'atlas.mat'), '-struct', 'at');
            mri = struct('SCS', struct('R', eye(3), 'T', [0;0;0]), 'NCS', struct('R', eye(3), 'T', [0;0;0]));
            save(fullfile(d, 'mri.mat'), '-struct', 'mri');
            % homotopic fibres: vertex i of L to vertex i of R
            k = (1:4:nh)';  Points = zeros(numel(k), 2, 3);
            Points(:, 1, :) = S.Vertices(k, :);  Points(:, 2, :) = S.Vertices(k + nh, :);
            save(fullfile(d, 'fibers_subject.mat'), 'Points');
            % the template: the same cortex as a Brainstorm tess file, its MRI shifted 15 mm in MNI
            tc.Tmpl = fullfile(tc.Root, 'template');  mkdir(tc.Tmpl);
            sc = struct('Vertices', {S.Hemi{1}', S.Hemi{2}'}, 'Label', {'Cortex L', 'Cortex R'}, 'Region', {'LU', 'RU'});
            dk = at.atlas(1);  dks = struct('Vertices', {}, 'Label', {});
            for j = 1:numel(dk.Label), dks(j) = struct('Vertices', find(dk.Membership(j, :)), 'Label', dk.Label{j}); end
            Atlas = struct('Name', {'Structures', 'Desikan-Killiany'}, 'Scouts', {sc, dks});
            Vertices = S.Vertices;  Faces = S.Faces;  VertNormals = S.VertNormals;  Comment = 'synth template';
            Reg = struct('Sphere', struct('Vertices', S.Sphere));
            save(fullfile(tc.Tmpl, 'tess_cortex_pial_low.mat'), 'Vertices', 'Faces', 'VertNormals', 'Comment', 'Reg', 'Atlas');
            SCS = mri.SCS;  NCS = struct('R', eye(3), 'T', [15; 0; 0]);
            save(fullfile(tc.Tmpl, 'subjectimage_T1.mat'), 'SCS', 'NCS');
        end
    end

    methods (Test)

        function sphereIsExactRawIsNot(tc)
            [T, X] = rheome.scale.measure_atlas(tc.Name, [], "correspondence", Template=tc.Tmpl);
            tc.verifyEqual(X.dk_agree_sphere, [1; 1]);
            tc.verifyEqual(X.dk_agree_sphere_area, [1; 1], 'AbsTol', 1e-12);
            tc.verifyLessThan(X.dk_agree_raw, [1; 1], 'a 15 mm shift must cost raw-coordinate agreement');
            tc.verifyEqual(X.raw_nn_mm_median > 0, [true; true]);
            % on a sphere the trivial connection with poles at the poles IS the meridian frame
            tc.verifyLessThan(X.frame_dev_median_deg, [2; 2]);
            tc.verifyEqual(X.gauge_singular_faces, [2; 2]);
            tc.verifyEqual(height(T), 22);
        end

        function eventsAndSurrogates(tc)
            [T, X] = rheome.scale.measure_atlas(tc.Name, [], "atlasevents", Template=tc.Tmpl, NSurr=3);
            tc.verifyEqual(unique(X.depth)', 1:3);
            tc.verifyEqual(height(X), 2 * 3 * 4);                    % hemi x depth x (data + 3)
            r = T.value(T.metric == "events_per_min");
            tc.verifyEqual(numel(r), 6);  tc.verifyGreaterThanOrEqual(r, 0);
            p = T.value(T.metric == "p_rate");  tc.verifyTrue(all(p > 0 & p <= 1));
            % depth 1 has two tiles: an event can visit at most two
            tc.verifyLessThanOrEqual(max(X.tiles_visited_median(X.depth == 1)), 2);
        end

        function tensorIsTraceNormalised(tc)
            [T, X] = rheome.scale.measure_atlas(tc.Name, [], "gaugetensor", Template=tc.Tmpl);
            tc.verifyEqual(height(X), 16);                           % 8 depth-3 tiles per hemisphere
            tc.verifyEqual(X.t11 + X.t22, ones(16, 1), 'AbsTol', 1e-12);
            tc.verifyTrue(all(X.anisotropy >= 0 & X.anisotropy <= 1 + 1e-12));
            tc.verifyTrue(all(X.normal_share >= 0 & X.normal_share <= 1));
            tc.verifyEqual(T.value(T.metric == "n_tiles"), [8; 8]);
        end

        function dyadicRollsUpExactlyAtlasesDoNot(tc)
            [T, X] = rheome.scale.measure_atlas(tc.Name, [], "connectome", Gammas=[0 1], NSeeds=10, Kc=50);
            v = @(m, b) T.value(T.metric == m & T.band == b);
            tc.verifyLessThan(v("edge_loss", "dyadic_6_5_copower"), 1e-10);
            tc.verifyLessThan(v("edge_loss", "dyadic_6_5_fibres"), 1e-10);
            tc.verifyEqual(v("area_outside", "dyadic_6_5_copower"), 0, 'AbsTol', 1e-12);
            tc.verifyGreaterThan(v("area_outside", "destrieux_dk_copower"), 0);
            tc.verifyGreaterThan(v("edge_loss", "destrieux_dk_copower"), 1e-3);
            tc.verifyEqual(v("own_tractography", "subject-specific"), 1);
            tc.verifyEqual(v("interhemispheric_weight", "subject-specific"), 1, 'AbsTol', 1e-12);
            % gamma 0 is the surface bank: its atoms stay in their hemisphere; the connectome bank's cross.
            % ⚠ Not exactly 0 here: the two synthetic hemispheres are the same sphere, every eigenspace pairs
            % across them, and eigs returns the 50 modes with ~2 % of a pair mixed at the cut. On a cortex
            % the hemispheres differ and the surface atoms put exactly 0 across (VR1 2.1d).
            tc.verifyLessThan(max(X.contra_frac(X.gamma_mult == 0)), 0.05);
            tc.verifyGreaterThan(max(X.contra_frac(X.gamma_mult == 1)), 2 * max(X.contra_frac(X.gamma_mult == 0)));
            tc.verifyTrue(all(X.rms_width_mm > 0));
        end

        function throughTheDriver(tc)
            old = getenv('RHEOME_TEMPLATE');  setenv('RHEOME_TEMPLATE', tc.Tmpl);  tc.addTeardown(@setenv, 'RHEOME_TEMPLATE', old);
            od = fullfile(tc.Root, 'out');
            R = rheome.scale.run(tc.Name, OutDir=od, Analyses=["correspondence" "gaugetensor"]);
            tc.verifyEqual(R.timing.status, ["ok"; "ok"], strjoin(R.timing.message, ' | '));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'correspondence.csv')));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'gaugetensor.csv')));
        end
    end
end

% Author: Diellor Basha, 2026
