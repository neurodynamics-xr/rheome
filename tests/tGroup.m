classdef tGroup < matlab.unittest.TestCase
% rheome.group on synthetic tile maps: the GLM's t against OLS, max-|t| FWER on null data and on a planted
% tile, the summary-statistics mixed model, the spin test's false rate and detection, and calibrate end to end.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function glmTMatchesOls(tc)
            rs = RandStream('mt19937ar', 'Seed', 1);  n = 30;
            X = [ones(n, 1) (1:n)' > n/2 randn(rs, n, 1)];  Y = randn(rs, n, 4) + 0.5 * X(:, 2);
            R = rheome.group.glm(Y, X, 2, NPerm=0);
            for k = 1:4
                m = fitlm(X(:, 2:3), Y(:, k));  tc.verifyEqual(R.t(k), m.Coefficients.tStat(2), 'RelTol', 1e-10);
            end
            tc.verifyEqual(R.df, n - 3);
        end

        function glmFwerOnPlantedTile(tc)
            rs = RandStream('mt19937ar', 'Seed', 2);  n = 40;  g = double((1:n)' > n/2);
            Y = randn(rs, n, 12);  Y(:, 1) = Y(:, 1) + 3 * g;
            R = rheome.group.glm(Y, [ones(n, 1) g], 2, NPerm=999);
            tc.verifyEqual(R.p_fwer(1), 1/1000, 'AbsTol', 1e-12);
            tc.verifyLessThan(R.p_fwer(1), min(R.p_fwer(2:end)));   % the null tiles' FWER: glmNullFwerCalibrated
            tc.verifyGreaterThan(median(R.p_fwer(2:end)), 0.2);
            tc.verifyTrue(all(R.p_fwer >= R.p - 1e-12));         % FWER p is never below the tile's own p
        end

        function glmNullFwerCalibrated(tc)
            % 300 null data sets, 10 correlated tiles: the family-wise false rate stays near 0.05
            rs = RandStream('mt19937ar', 'Seed', 3);  n = 30;  m = 10;  hit = 0;  N = 300;
            L = chol(0.5 * eye(m) + 0.5);
            for d = 1:N
                Y = randn(rs, n, m) * L;  z = randn(rs, n, 1);
                R = rheome.group.glm(Y, [ones(n, 1) z], 2, NPerm=199, Seed=d);  hit = hit + any(R.p_fwer < 0.05);
            end
            [~, ci] = binofit(hit, N);
            tc.verifyLessThanOrEqual(ci(1), 0.05);  tc.verifyGreaterThanOrEqual(ci(2), 0.05 - 0.02);
        end

        function glmSubjectIsSummaryStatistics(tc)
            rs = RandStream('mt19937ar', 'Seed', 4);  ns = 16;  g = double((1:ns)' > ns/2);
            sid = repelem((1:ns)', 3);  Y = randn(rs, 3*ns, 5) + repelem(g, 3);  X = [ones(3*ns, 1) repelem(g, 3)];
            A = rheome.group.glm(Y, X, 2, NPerm=99, Subject=sid);
            Ym = splitapply(@(y) mean(y, 1), Y, sid);  B = rheome.group.glm(Ym, [ones(ns, 1) g], 2, NPerm=99);
            tc.verifyEqual(A.t, B.t, 'RelTol', 1e-12);  tc.verifyEqual(A.p_fwer, B.p_fwer);  tc.verifyEqual(A.n, ns);
            Xbad = X;  Xbad(1, 2) = 1 - Xbad(1, 2);
            tc.verifyError(@() rheome.group.glm(Y, Xbad, 2, Subject=sid), 'group:glm:within');
        end

        function glmUntestableTileIsNotSignificant(tc)
            % a tile exactly 0 or constant (round-off residuals) gets t NaN, p = p_fwer = 1 and stays out of the max
            rs = RandStream('mt19937ar', 'Seed', 7);  n = 30;  g = double((1:n)' > n/2);
            Y = randn(rs, n, 5);  Y(:, 1) = Y(:, 1) + 3 * g;  Y(:, 3) = 0;  Y(:, 4) = 0.3;
            R = rheome.group.glm(Y, [ones(n, 1) g], 2, NPerm=199);
            tc.verifyEqual(R.untestable, [false false true true false]);
            tc.verifyTrue(all(isnan(R.t(3:4))));
            tc.verifyEqual([R.p(3:4) R.p_fwer(3:4)], ones(1, 4));
            tc.verifyTrue(all(isfinite(R.maxnull)));  tc.verifyEqual(R.p_fwer(1), 1/200, 'AbsTol', 1e-12);
        end

        function glmFreedmanLaneNuisanceNull(tc)
            % the tested x correlated 0.7 with a nuisance c of strong effect (2 sd): FWER holds only if the reduced
            % model's residuals are permuted (permuting Y would carry c's effect into the null and break it)
            rs = RandStream('mt19937ar', 'Seed', 8);  n = 30;  m = 8;  hit = 0;  N = 300;
            L = chol(0.5 * eye(m) + 0.5);
            for d = 1:N
                x = randn(rs, n, 1);  c = 0.7 * x + sqrt(1 - 0.49) * randn(rs, n, 1);
                Y = randn(rs, n, m) * L + 2 * c;
                R = rheome.group.glm(Y, [ones(n, 1) x c], 2, NPerm=199, Seed=d);  hit = hit + any(R.p_fwer < 0.05);
            end
            [~, ci] = binofit(hit, N);
            tc.verifyLessThanOrEqual(ci(1), 0.05);  tc.verifyGreaterThanOrEqual(ci(2), 0.05 - 0.02);
        end

        function spinDetectsAndHoldsFalseRate(tc)
            rs = RandStream('mt19937ar', 'Seed', 5);  m = 300;
            U = randn(rs, m, 3);  U = U ./ vecnorm(U, 2, 2);  h = repmat("L", m, 1);  h(U(:, 1) > 0) = "R";
            bump = @(c) exp(-acosd(min(1, U * c')).^2 / (2 * 30^2));
            ref = bump([0 0 1]);
            S = rheome.group.spin(ref + 0.05 * randn(rs, m, 1), ref, U, h, NSpin=499, Seed=1);
            tc.verifyLessThan(S.p, 0.01);
            S2 = rheome.group.spin(ref, ref, U, h, Index=S.index);
            tc.verifyEqual(size(S2.rnull), [499 1]);
            hit = 0;  N = 100;
            for d = 1:N                                      % independent smooth maps: no correspondence
                c = randn(rs, 3, 3);  c = c ./ vecnorm(c, 2, 2);
                M = bump(c(1, :)) - bump(c(2, :)) + bump(c(3, :)) + 0.1 * randn(rs, m, 1);
                Sd = rheome.group.spin(M, ref, U, h, Index=S.index);  hit = hit + (Sd.p < 0.05);
            end
            [~, ci] = binofit(hit, N);  tc.verifyLessThanOrEqual(ci(1), 0.05);
            F = rheome.group.spin(ref, [ref bump([1 0 0])], U, h, Index=S.index);
            tc.verifyTrue(all(F.p_fwer >= F.p));
            Mn = ref;  Mn(1) = NaN;                          % a NaN in the map is not evidence
            N1 = rheome.group.spin(Mn, ref, U, h, Index=S.index);  tc.verifyEqual([N1.p N1.p_fwer], [1 1]);
        end

        function calibrateEndToEnd(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            rs = RandStream('mt19937ar', 'Seed', 6);  hemi = repelem(["L"; "R"], 6);  node = repmat((8:13)', 2, 1);
            U = randn(rs, 12, 3);  U = U ./ vecnorm(U, 2, 2);  U(1:6, 1) = -abs(U(1:6, 1));  U(7:12, 1) = abs(U(7:12, 1));
            Ct = table(hemi, node, U(:, 1), U(:, 2), U(:, 3), 'VariableNames', {'hemi','node_id','u1','u2','u3'});
            for ds = ["omega" "camcan"]                      % the same subject IDs in two cohorts: two participants each
                p = fullfile(d, "cf-atlas-" + ds);  mkdir(p);
                for s = 1:24
                    pre = sprintf('sub-%02d_', s);
                    t = table(hemi, repmat(3, 12, 1), node, rand(rs, 12, 1), rand(rs, 12, 1), 1 + rand(rs, 12, 1), ...
                              'VariableNames', {'hemi','depth','node_id','anisotropy','normal_share','trace'});
                    writetable(t, fullfile(p, [pre 'gaugetensor.csv']));
                    writetable(table(string(pre), ds, 'VariableNames', {'subject','dataset'}), fullfile(p, [pre 'metrics.csv']));
                end
            end
            o = fullfile(d, 'out');
            [C, D] = rheome.group.calibrate(d, o, Centres=Ct, NDraw=4, NPerm=99, NSpin=99, Effects=[0 4], Measures="anisotropy");
            for f = ["draws.csv" "calibration.csv" "tiles.csv" "participants.csv"], tc.verifyTrue(isfile(fullfile(o, f))); end
            tc.verifyEqual(height(D), 2 * 4 * 2 * 2);  tc.verifyEqual(height(C), 2 * 2 * 2);
            tc.verifyEqual(unique(D.n), 24);  tc.verifyEqual(unique(D.n_tiles), 12);  tc.verifyEqual(unique(D.nuisance), 0.5);
            tc.verifyEqual(C.detected_rate(C.effect == 4), ones(4, 1));   % 4 sd at the seed tile, 12 vs 12
            tc.verifyEqual(unique(C.K), 2 * 2 + 2 * 2);               % fp_outside x 4 conditions + any_fwer, spin x 2
            tc.verifyTrue(all(ismember({'any_fwer_lo','fp_outside_hi','spin_rate','gate_lo_spin','calibrated'}, C.Properties.VariableNames)));
            % paired conditions: every Transform x Nuisance on the same draws; rankit removes the skew
            [C2, D2] = rheome.group.calibrate(d, fullfile(d, 'out2'), Centres=Ct, NDraw=2, NPerm=19, NSpin=19, Effects=0, ...
                                              Measures="log10_trace", Nuisance=[0 0.5], Transform=["none" "rankit"]);
            tc.verifyEqual(height(D2), 2 * 2 * 2 * 2 * 2);  tc.verifyEqual(height(C2), 2 * 2 * 2 * 2);
            k = D2(D2.dataset == "omega" & D2.draw == 1 & D2.kind == "group", :);
            tc.verifyEqual(numel(unique(k.seed_tile)), 1);           % the same draw across conditions
            Sh = readtable(fullfile(d, 'out2', 'shape.csv'), 'TextType', 'string');
            tc.verifyEqual(height(Sh), 2 * 2);
            tc.verifyLessThan(Sh.skew_absmax(Sh.transform == "rankit"), 0.05 * [1; 1]);
            delete(fullfile(d, 'cf-atlas-camcan', 'sub-03_metrics.csv'));   % a participant without a cohort
            tc.verifyError(@() rheome.group.calibrate(d, o, Centres=Ct, NDraw=1, NPerm=9, NSpin=9, Measures="anisotropy"), ...
                           'group:calibrate:cohort');
        end
    end
end

% Author: Diellor Basha, 2026
