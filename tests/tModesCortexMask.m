classdef tModesCortexMask < matlab.unittest.TestCase
% Pattern atlas Amendment 3 (design 52798c1e section 24; decision 23628c1c; request 99da1c4f): the atlas runs on
% cortex only (decision 3882c4ce). On the synthetic cached subject (scaleSynthSubject, its own template) with one
% depth-3 group tile of the left hemisphere declared medial wall ('unknown L'): rheome.modes.operator gives that
% tile W row 0, area 0 and a finite centroid, every other tile a finite cortex-area mean, and no NaN reaches the
% tile resolution R (before the fix, a wall-only tile -- 3 of 128 at the joint F's depth 6 on OMEGA sub-0002 --
% made its W row 0/0 = NaN).
%
% Author: Diellor Basha, 2026

    properties
        Root
        Name = "synth_mask1"
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            d = fullfile(tc.Root, tc.Name);  scaleSynthSubject(d);
            sf = load(fullfile(d, 'surface.mat'));  sf.S.Sphere = 0.1 * sf.S.VertNormals;  save(fullfile(d, 'surface.mat'), '-struct', 'sf');
            S = sf.S;  tm = fullfile(tc.Root, 'template');  mkdir(tm);
            sc = struct('Vertices', {S.Hemi{1}', S.Hemi{2}'}, 'Label', {'Cortex L', 'Cortex R'}, 'Region', {'LU', 'RU'});
            Atlas = struct('Name', {'Structures', 'Desikan-Killiany'}, 'Scouts', {sc, struct('Vertices', {}, 'Label', {})});
            Vertices = S.Vertices;  Faces = S.Faces;  VertNormals = S.VertNormals;  Comment = 'synth template'; %#ok<NASGU>
            Reg = struct('Sphere', struct('Vertices', S.Sphere)); %#ok<NASGU>
            save(fullfile(tm, 'tess_cortex_pial_low.mat'), 'Vertices', 'Faces', 'VertNormals', 'Comment', 'Reg', 'Atlas');
            old = getenv('RHEOME_TEMPLATE');  setenv('RHEOME_TEMPLATE', tm);  tc.addTeardown(@setenv, 'RHEOME_TEMPLATE', old);
            % the wall: every vertex of the first depth-3 group tile of the left hemisphere
            S0 = rheome.scale.sensors(char(tc.Name));  Tp = rheome.scale.template(3);
            P = rheome.scale.grouptiles(S0, Tp, "L", 3);  gv = S0.B.L.gv(:);  w = false(S.nV, 1);  w(gv(P(:, 1) > 0)) = true;
            L = load(fullfile(d, 'atlas.mat'));  a = L.atlas;
            a.Label = [{'unknown L'}, a.Label];  a.Membership = [sparse(w'); a.Membership .* ~w'];
            atlas = a; %#ok<NASGU>
            save(fullfile(d, 'atlas.mat'), 'atlas');
        end
    end

    methods (Test)
        function wallOnlyTileReadsZero(tc)
            S = rheome.scale.sensors(char(tc.Name));  Tp = rheome.scale.template(3);
            tc.assertTrue(any(~S.Sf.Cortex) && all(S.Sf.Cortex(S.Sf.Hemi{2})), 'the wall was not masked');
            Op = rheome.modes.operator(S, Tp, CutoffMM=1, Depth=3, GenCov=true);
            t = Op.tiles;  wall = t.area_mm2 == 0;
            tc.verifyEqual(nnz(wall), 1);                                             % exactly the declared tile
            tc.verifyEqual(t.hemi(wall), "L");
            tc.verifyTrue(all(isfinite(Op.W), 'all') && all(isfinite(Op.R), 'all'));
            tc.verifyEqual(Op.W(wall, :), zeros(1, size(Op.W, 2)));
            tc.verifyTrue(all(any(Op.W(~wall, :), 2)));
            tc.verifyTrue(all(isfinite([t.x_mm t.y_mm t.z_mm]), 'all'));
            tc.verifyGreaterThan(norm([t.x_mm(wall) t.y_mm(wall) t.z_mm(wall)]), 0);   % its closed centroid, not 0
            tc.verifyEqual(squeeze(sum(abs(Op.Cgen(:, :, wall)), [1 2])), 0);        % no generator on the wall
            % cortex area: the tiles' areas sum to the left hemisphere's cortex
            m = full(sum(S.B.L.lbo.Mass, 2));
            tc.verifyEqual(sum(t.area_mm2(t.hemi == "L")), 1e6 * sum(m(S.B.L.S.Cortex)), 'RelTol', 1e-9);
        end
    end
end

% Author: Diellor Basha, 2026
