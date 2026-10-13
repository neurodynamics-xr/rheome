classdef tModes < matlab.unittest.TestCase
% The resting-state pattern atlas (rheome.modes): the building blocks on planted data -- phase-invariant
% clustering recovers planted complex modes under random global phases, Hungarian matching, the centroid
% hash, phase randomisation keeps every cross-spectrum, the analytic band signal and the periodic split --
% and the four stages end to end on three copies of the synthetic cached subject (scaleSynthSubject) with
% the subject as its own template: fields (with nulls), discover (frozen + hashed), fit, reduce.
%
% Author: Diellor Basha, 2026

    properties
        Root
        Tmpl
        Names = ["synth_modes1" "synth_modes2" "synth_modes3"]
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            for n = tc.Names
                d = fullfile(tc.Root, n);  scaleSynthSubject(d);
                sf = load(fullfile(d, 'surface.mat'));  sf.S.Sphere = 0.1 * sf.S.VertNormals;  save(fullfile(d, 'surface.mat'), '-struct', 'sf');
            end
            S = sf.S;  tc.Tmpl = fullfile(tc.Root, 'template');  mkdir(tc.Tmpl);
            sc = struct('Vertices', {S.Hemi{1}', S.Hemi{2}'}, 'Label', {'Cortex L', 'Cortex R'}, 'Region', {'LU', 'RU'});
            Atlas = struct('Name', {'Structures', 'Desikan-Killiany'}, 'Scouts', {sc, struct('Vertices', {}, 'Label', {})});
            Vertices = S.Vertices;  Faces = S.Faces;  VertNormals = S.VertNormals;  Comment = 'synth template'; %#ok<NASGU>
            Reg = struct('Sphere', struct('Vertices', S.Sphere)); %#ok<NASGU>
            save(fullfile(tc.Tmpl, 'tess_cortex_pial_low.mat'), 'Vertices', 'Faces', 'VertNormals', 'Comment', 'Reg', 'Atlas');
            old = getenv('RHEOME_TEMPLATE');  setenv('RHEOME_TEMPLATE', tc.Tmpl);  tc.addTeardown(@setenv, 'RHEOME_TEMPLATE', old);
        end
    end

    methods (Test)

        function clusterRecoversPlantedModes(tc)
            rng(3);  d = 16;  K = 4;  n = 4000;
            A = randn(d, K) + 1i*randn(d, K);  A = A ./ vecnorm(A);
            lab = randi(K, 1, n);  X = A(:, lab) .* exp(2i*pi*rand(1, n)) .* (0.5 + rand(1, n));
            X = X + 0.15 * (randn(d, n) + 1i*randn(d, n)) / sqrt(2*d);
            [C, l2, info] = rheome.modes.cluster(X, K, Seed=1);
            [r, perm] = rheome.modes.match(A, C);
            tc.verifyGreaterThan(min(r), 0.97);
            [~, inv] = sort(perm);  tc.verifyGreaterThan(mean(reshape(inv(l2), 1, []) == lab), 0.95);
            tc.verifyTrue(info.converged);  tc.verifyEqual(vecnorm(C), ones(1, K), 'AbsTol', 1e-12);
            [~, j] = max(abs(C));  tc.verifyEqual(imag(C(sub2ind(size(C), j, 1:K))), zeros(1, K), 'AbsTol', 1e-12);
        end

        function matchIsPhaseInvariantHungarian(tc)
            rng(4);  A = randn(8, 5) + 1i*randn(8, 5);  p = [3 5 1 2 4];
            B = A(:, p) .* exp(1i * (1:5));  [r, perm] = rheome.modes.match(A, B);
            tc.verifyEqual(r, ones(5, 1), 'AbsTol', 1e-12);  tc.verifyEqual(B(:, perm) ./ A, B(:, perm) ./ A);
            [~, ip] = sort(p);  tc.verifyEqual(perm(:)', ip);
        end

        function hashIsStableAndSensitive(tc)
            C = [1+2i 3; 4 5i];
            tc.verifyEqual(rheome.modes.hash(C), rheome.modes.hash(single(C)));
            tc.verifyNotEqual(rheome.modes.hash(C), rheome.modes.hash(C + 1e-9));
            tc.verifyNotEqual(rheome.modes.hash(C), rheome.modes.hash(C(:)));
            tc.verifyEqual(strlength(string(rheome.modes.hash(C))), 64);
        end

        function phaserandKeepsEveryCrossSpectrum(tc)
            rng(5);  N = 1001;  x = randn(1, N);  Y = [x; circshift(x, 7); randn(1, N)];
            Ys = rheome.modes.phaserand(Y, RandStream('mt19937ar', 'Seed', 1));
            F = fft(Y, [], 2);  Fs = fft(Ys, [], 2);
            tc.verifyEqual(abs(Fs), abs(F), 'RelTol', 1e-9, 'AbsTol', 1e-9);
            tc.verifyEqual(Fs(1, 2:end) .* conj(Fs(2, 2:end)), F(1, 2:end) .* conj(F(2, 2:end)), 'RelTol', 1e-6, 'AbsTol', 1e-6);
            tc.verifyGreaterThan(norm(Ys - Y, 'fro'), 0.5 * norm(Y, 'fro'));
        end

        function laggedSeparatesTwoPhasesFromAGradient(tc)
            % identity blur: a two-phase source fits exactly, a smooth phase ramp over 8 tiles does not
            d = 8;  R = eye(d);
            m2 = [ones(4, 1); 0.5 * exp(1i * pi/3) * ones(4, 1)];
            L2 = rheome.modes.lagged(m2, R);  tc.verifyEqual(L2.r_n3, 1, 'AbsTol', 1e-3);
            tc.verifyEqual(L2.lag_rad, pi/3, 'AbsTol', deg2rad(10));
            tc.verifyTrue(all(L2.x1 .* L2.x2 == 0));                               % disjoint generators
            mg = exp(1i * linspace(0, 1.5*pi, d)');  Lg = rheome.modes.lagged(mg, R);
            tc.verifyLessThan(Lg.r_n3, 0.95);
        end

        function analyticBandSignal(tc)
            fs = 200;  t = (0:fs*20-1)/fs;  Y = [2*cos(2*pi*10*t); cos(2*pi*10*t + 1)];
            [z, rate, valid] = rheome.modes.analytic(Y, fs, [8 13]);
            tc.verifyGreaterThanOrEqual(rate, 3*13);  tc.verifyEqual(size(z, 2), numel(valid));
            tc.verifyEqual(double(median(abs(z(1, valid)))), 2, 'RelTol', 0.02);
            tc.verifyEqual(double(median(angle(z(1, valid) .* conj(z(2, valid))))), -1, 'AbsTol', 0.02);
            tc.verifyFalse(valid(1));  tc.verifyFalse(valid(end));
        end

        function periodicSplitFindsTheRhythm(tc)
            rs = RandStream('mt19937ar', 'Seed', 6);  fs = 200;  N = fs*120;  f = (0:N-1)/N*fs;  f = min(f, fs - f);
            pink = real(ifft(fft(randn(rs, 2, N), [], 2) ./ sqrt(max(f, 0.5)), [], 2));
            t = (0:N-1)/fs;  Y = pink + [0.3; 0] .* cos(2*pi*10*t);
            sp = rheome.modes.split(Y, fs);
            pa = @(x) bandpower(x', fs, [9 11]);
            tc.verifyGreaterThan(pa(sp.per(1, :)), 4 * pa(sp.per(2, :)));         % the rhythm is periodic
            tc.verifyLessThan(pa(sp.per(2, :)), 0.5 * pa(Y(2, :)));                % 1/f alone is mostly aperiodic
        end

        function fourStagesEndToEnd(tc)
            od = fullfile(tc.Root, 'out');
            for n = tc.Names
                R = rheome.scale.run(char(n), OutDir=fullfile(od, 'fields'), Analyses="modesfields");
                tc.verifyEqual(R.timing.status, "ok", strjoin(R.timing.message, ' | '));
            end
            files = fullfile(od, 'fields', tc.Names, 'modes_fields.mat');
            m = matfile(files(1));  Z = m.Z;
            tc.verifyEqual(Z.variants.cutoff_mm', [100 92 130 100 100]);
            tc.verifyTrue(all(ismember(["data" "n1_sources" "n4_aperiodic" "n5_phaserand"], Z.sources)), strjoin(string(Z.sources), ' '));
            lst = fullfile(od, 'fields.txt');  writelines(files, lst);              % a list file, as nsp passes it
            M = rheome.modes.discover(string(lst), fullfile(od, 'discover'), Ks=2:3, NSplits=4, NPer=200, NPerSplit=100, Replicates=2);
            tc.verifyEqual(M.sha256, rheome.modes.hash(M.bands.C));
            tc.verifyEqual(strtrim(extractBefore(string(fileread(fullfile(od, 'discover', 'centroids.sha256'))), ' ')), string(M.sha256));
            [~, X] = rheome.modes.fit(char(tc.Names(1)), [], M, NSpin=20, NSurr=3);
            v = X.modesfit;
            tc.verifyTrue(all(ismember(["data" "n1_sources" "n4_aperiodic" "n5_phaserand"], unique(v.source)')), strjoin(string(unique(v.source)'), ' '));
            tc.verifyTrue(all(ismember(["occupancy" "dwell_median_ms" "gev" "p_spin" "pg_n3" "extent_mm2" "dirichlet_wavelength_mm" "r_halves"], unique(v.metric)')), strjoin(string(unique(v.metric)'), ' '));
            occ = v.value(v.source == "data" & v.half == "all" & v.method == "backfit" & v.metric == "occupancy");
            tc.verifyEqual(sum(occ), numel(M.bands), 'AbsTol', 1e-9);              % occupancies sum to 1 per band
            tc.verifyTrue(all(v.value(v.metric == "speed_resolved") == 0));        % depth 3: never a speed
            cal = v(v.method == "n3cal" & v.metric == "n3_excess", :);                % section 13 calibration, every band
            tc.verifyEqual(height(cal), 3 * numel(M.bands));  tc.verifyTrue(all(isfinite(cal.value)));
            bad = M;  bad.bands(1).C(1) = bad.bands(1).C(1) + 1e-6;
            tc.verifyError(@() rheome.modes.fit(char(tc.Names(1)), [], bad), 'modes:fit:hash');
            % a state the back-fit gives no sample (a duplicate centroid loses every tie) has no N3: NaN, not
            % r_n3 0 / r_n3_sim 1; with one sample its mode is NaN and rheome.modes.lagged crashed
            em = M;  for b = 1:numel(em.bands), em.bands(b).C(:, end+1) = em.bands(b).C(:, 1); end
            em.sha256 = rheome.modes.hash(em.bands.C);
            [~, Xe] = rheome.modes.fit(char(tc.Names(1)), [], em, NSpin=20, NSurr=3);  ve = Xe.modesfit;
            ke = ve.source == "data" & ve.half == "all" & ve.method == "backfit" & ve.mode == size(em.bands(1).C, 2) & ve.band == em.bands(1).band;
            tc.verifyEqual(ve.value(ke & ve.metric == "occupancy"), 0);
            tc.verifyTrue(all(isnan(ve.value(ke & ismember(ve.metric, ["r_n3" "r_n3_sim" "n3_lag_rad" "n3_tiles_a"])))));
            tc.verifyEqual(nnz(ke & ve.metric == "r_n3"), 1);
            tc.verifyTrue(isnan(rheome.modes.lagged(nan(4, 1), eye(4)).r_n3));      % a non-finite mode: NaN, no crash
            % the driver writes the fit tables, then the reduction over the three
            setenv('RHEOME_MODES_CENTROIDS', fullfile(od, 'discover', 'centroids.mat'));  tc.addTeardown(@setenv, 'RHEOME_MODES_CENTROIDS', '');
            fits = strings(0);  sa = strings(0);  sb = strings(0);
            for n = tc.Names
                fo = fullfile(od, 'fit', n);  mkdir(fo);
                [~, Xn] = rheome.modes.fit(char(n), [], M, NSpin=20, NSurr=3);
                t = [table(repmat(n, height(Xn.modesfit), 1), 'VariableNames', {'subject'}), Xn.modesfit];
                writetable(t, fullfile(fo, 'modesfit.csv'));  fits(end+1) = fullfile(fo, 'modesfit.csv'); %#ok<AGROW>
                sa(end+1) = fullfile(fo, n + "_modesseries.csv");  writetable(Xn.modesseries, sa(end)); %#ok<AGROW>
                sb(end+1) = fullfile(od, 'fields', n, n + "_modesseries.csv"); %#ok<AGROW>
                copyfile(fullfile(od, 'fields', n, 'modesseries.csv'), sb(end));
            end
            % the fit recomputes the discovery's series: same commit, same hashes; a changed series is caught
            [T, ~] = rheome.modes.seriescheck(sa, sb, fullfile(od, 'seriescheck'));
            tc.verifyTrue(all(T.equal), 'the fit''s series differ from stage 1''s on one commit');
            tc.verifyEqual(numel(unique(T.subject)), numel(tc.Names));
            tc.verifyTrue(all(T.variant == 1));                                     % sensitivity variants not compared
            t = readtable(sa(1), 'TextType', 'string');  t.sha256(1) = "x";  writetable(t, sa(1));
            T = rheome.modes.seriescheck(sa, sb, fullfile(od, 'seriescheck2'));
            tc.verifyEqual(nnz(~T.equal), 1);
            j = jsondecode(fileread(fullfile(od, 'seriescheck2', 'seriescheck.json')));
            tc.verifyFalse(j.all_equal);  tc.verifyEqual(string(j.subjects_differ), tc.Names(1));
            lst = fullfile(od, 'fits.txt');  writelines(fits, lst);
            C = rheome.modes.reduce(string(lst), fullfile(od, 'discover', 'centroids.mat'), fullfile(od, 'reduce'), NBoot=50);
            tc.verifyEqual(height(C), sum([M.bands.K]));
            tc.verifyTrue(all(ismember(C.label, ["standing" "lagged-or-rotating"])));   % no propagation label (section 13)
            tc.verifyTrue(ismember("n3_excess_median", string(C.Properties.VariableNames)));
            tc.verifyTrue(all(C.speed == "unresolved"));
            tc.verifyTrue(isfile(fullfile(od, 'reduce', 'modes_tests.csv')));
            % Level 2 shares for dwell AND transition entropy, with the N5 floor: p >= 1/(1 + NSurr)
            tc.verifyTrue(all(ismember(["l2_entropy_share_sig" "l2_entropy_share_p" "l2_p_min" "l2_n_p_min" "l2_entropy_p_min" ...
                "l2_entropy_n_p_min"], string(C.Properties.VariableNames))));
            tc.verifyGreaterThanOrEqual([C.l2_p_min; C.l2_entropy_p_min], 1/4 - 1e-12);
            tc.verifyGreaterThanOrEqual(C.l2_entropy_n_p_min, 1);
        end
    end
end

% Author: Diellor Basha, 2026
