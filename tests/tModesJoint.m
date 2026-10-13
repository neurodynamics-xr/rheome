classdef tModesJoint < matlab.unittest.TestCase
% H4 ATLAS-1 (protocol section 14): the joint graph-time atlas stages on three copies of the synthetic cached
% subject (scaleSynthSubject, the subject as its own template, as tModes) -- the octave member, the scale
% transfer T(s), the joint fields with null F, discover (NMF, t_bar) -> fit -> reduce on joint atoms, and the
% J5 canary path (plants on F, rheome.modes.jointcheck). The synthetic head checks the code paths, not the numbers.
%
% Author: Diellor Basha, 2026

    properties
        Root
        Names = ["synth_joint1" "synth_joint2" "synth_joint3"]
        Opt = {'Depths', 0:2, 'GenDepth', 3, 'Octaves', [4 8; 8 16]}
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            for n = tc.Names
                d = fullfile(tc.Root, n);  scaleSynthSubject(d);
                sf = load(fullfile(d, 'surface.mat'));  sf.S.Sphere = 0.1 * sf.S.VertNormals;  save(fullfile(d, 'surface.mat'), '-struct', 'sf');
            end
            S = sf.S;  tm = fullfile(tc.Root, 'template');  mkdir(tm);
            sc = struct('Vertices', {S.Hemi{1}', S.Hemi{2}'}, 'Label', {'Cortex L', 'Cortex R'}, 'Region', {'LU', 'RU'});
            Atlas = struct('Name', {'Structures', 'Desikan-Killiany'}, 'Scouts', {sc, struct('Vertices', {}, 'Label', {})});
            Vertices = S.Vertices;  Faces = S.Faces;  VertNormals = S.VertNormals;  Comment = 'synth template'; %#ok<NASGU>
            Reg = struct('Sphere', struct('Vertices', S.Sphere)); %#ok<NASGU>
            save(fullfile(tm, 'tess_cortex_pial_low.mat'), 'Vertices', 'Faces', 'VertNormals', 'Comment', 'Reg', 'Atlas');
            old = getenv('RHEOME_TEMPLATE');  setenv('RHEOME_TEMPLATE', tm);  tc.addTeardown(@setenv, 'RHEOME_TEMPLATE', old);
        end
    end

    methods (Test)
        function octaveIsTheAnalyticBand(tc)
            fs = 200;  t = (0:fs*20-1) / fs;  y = cos(2*pi*10*t) + cos(2*pi*3*t);
            [z, rate, valid] = rheome.modes.octave(y, fs, [8 16]);
            tc.verifyEqual(rate, fs / floor(fs / 48));
            tc.verifyEqual(median(abs(z(valid))), single(1), 'RelTol', 0.02);          % the 10 Hz line, analytic, unit
            ph = unwrap(angle(double(z(valid))));  tc.verifyEqual(median(diff(ph)) * rate / (2*pi), 10, 'RelTol', 0.01);
        end

        function scaleTransferIsAShareAndHasAnEdge(tc)
            S = rheome.scale.sensors(char(tc.Names(1)));
            win = table((0:2)', [150; 60; 20], [Inf; 150; 60], [200; 90; 35], 'VariableNames', {'depth','lo_mm','hi_mm','scale_mm'});
            Tr = rheome.modes.scaletransfer(S, win);
            T = Tr.table;  tc.verifyTrue(all(T.share >= 0 & T.share <= 1.05 | isnan(T.share)));
            tc.verifyTrue(all(T.T(isfinite(T.T)) <= 1 + 1e-9));                    % normalised by the bank's peak
            tc.verifyTrue(ismember(Tr.edge_flag, ["interpolated" "beyond_finest" "below_edge"]));
        end

        function jointStagesEndToEnd(tc)
            od = fullfile(tc.Root, 'out');  files = strings(0);
            for n = tc.Names
                [Z, Y] = rheome.modes.jointfields(char(n), [], tc.Opt{:});
                f = fullfile(od, n + "_joint.mat");  if ~exist(od, 'dir'), mkdir(od); end
                save(f, '-struct', 'Y', '-v7.3');  save(f, 'Z', '-append');  files(end+1) = f; %#ok<AGROW>
            end
            tc.verifyEqual(Z.kind, "joint");
            tc.verifyEqual(Z.sources, ["data" "f_fieldspread" "f_fieldspread2" "f_fieldspread3" "n4_aperiodic" "n5_phaserand"]);
            tc.verifyEqual(height(Z.tiles{1}), sum(Z.scales.n_tiles));
            tc.verifyTrue(all(isfinite(Z.transfer.table.T_abs)));
            tc.verifyTrue(isfinite(Z.fspread.tile_logpower_err_median));
            M = rheome.modes.discover(files, fullfile(od, 'discover'), NullSources=["f_fieldspread" "f_fieldspread2" "f_fieldspread3" "n4_aperiodic"], NMF=true, ...
                                      Ks=2:3, NSplits=4, NPer=200, NPerSplit=100, Replicates=2);
            tc.verifyEqual(M.kind, "joint");
            B = M.bands(1);                                                          % F union: every K, both realisations
            tc.verifyEqual(B.nullK.f_fieldspread, [2;2;3;3;3]);  tc.verifyEqual(size(B.nullC.f_fieldspread2, 2), 5);
            tc.verifyEqual(size(B.nullC.n4_aperiodic, 2), B.K);                     % other nulls stay at K*
            tc.verifyTrue(B.r_min_f <= 0.7 && abs(B.r_min_f - min(0.7, B.f_low - 0.05)) < 1e-12);   % F calibration
            tc.verifyTrue(all(isfinite(M.bands(1).r_nmf)) && all(isfinite(M.bands(1).t_bar)));
            fits = strings(0);
            for n = tc.Names
                [~, X] = rheome.modes.fit(char(n), [], M, NSpin=10, NSurr=2);  v = X.modesfit;
                fo = fullfile(od, 'fit', n);  mkdir(fo);
                writetable([table(repmat(n, height(v), 1), 'VariableNames', {'subject'}), v], fullfile(fo, 'modesfit.csv'));
                fits(end+1) = fullfile(fo, 'modesfit.csv'); %#ok<AGROW>
            end
            tc.verifyTrue(all(ismember(["r_f_fieldspread" "t_bar" "scale_mm" "dominant_depth" "p_spin" "r_n3"], unique(v.metric)')));
            tc.verifyFalse(any(v.method == "n3cal"));
            C = rheome.modes.reduce(fits, fullfile(od, 'discover', 'centroids.mat'), fullfile(od, 'reduce'), NBoot=20);
            tc.verifyTrue(all(C.label == "not_labelled"));
            Ts = readtable(fullfile(od, 'reduce', 'modes_tests.csv'), 'TextType', 'string');
            tc.verifyTrue(all(isnan(Ts.p_holm(Ts.test == "n3_lagged"))));             % reported, not gated
            tc.verifyTrue(any(Ts.test == "i_f"));
        end

        function canaryPath(tc)
            od = fullfile(tc.Root, 'canary');  mkdir(od);  files = strings(0);
            for n = tc.Names
                [Z, Y] = rheome.modes.jointfields(char(n), [], tc.Opt{:}, Canary=true);
                f = fullfile(od, n + "_joint.mat");  save(f, '-struct', 'Y', '-v7.3');  save(f, 'Z', '-append');  files(end+1) = f; %#ok<AGROW>
            end
            tc.verifyFalse(ismember("data", Z.sources));                              % no real-data series
            tc.verifyEqual(size(Y.canary_truth, 2), 4);
            R = rheome.modes.jointcheck(files, fullfile(od, 'check'), NSplits=4, Ks=2:4);
            tc.verifyTrue(islogical(R.pass) && R.f_checked);
            fr = readtable(fullfile(od, 'check', 'joint_canary_f.csv'), 'TextType', 'string');      % F3 vs calibrated r_min
            tc.verifyEqual(logical(fr.passes_level1_i), fr.r_max_to_f12 < fr.r_min_f);
            tc.verifyTrue(all(R.k_star >= R.k_rule) && (isnan(R.k_sep) || all(R.k_star >= R.k_sep)));   % J5 amendment
            M = rheome.modes.discover(files, fullfile(od, 'kmin'), Source="canary", NullSources="f_fieldspread", ...
                                      Bands=unique(Z.canary.band), NSplits=2, Ks=2:4, KMin=4);
            tc.verifyEqual(unique([M.bands.K]), 4);  tc.verifyTrue(all([M.bands.k_rule] <= 4));
            for f = ["joint_canary.csv" "joint_canary_f.csv" "joint_canary_byK.csv" "joint_canary_transfer.csv" "joint_canary.json"]
                tc.verifyTrue(isfile(fullfile(od, 'check', f)), f);
            end
        end
    end
end

% Author: Diellor Basha, 2026
