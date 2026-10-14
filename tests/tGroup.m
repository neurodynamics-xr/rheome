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
        end

        function calibrateEndToEnd(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            rs = RandStream('mt19937ar', 'Seed', 6);  hemi = repelem(["L"; "R"], 6);  node = repmat((8:13)', 2, 1);
            U = randn(rs, 12, 3);  U = U ./ vecnorm(U, 2, 2);  U(1:6, 1) = -abs(U(1:6, 1));  U(7:12, 1) = abs(U(7:12, 1));
            Ct = table(hemi, node, U(:, 1), U(:, 2), U(:, 3), 'VariableNames', {'hemi','node_id','u1','u2','u3'});
            for s = 1:24
                p = fullfile(d, 'cf-atlas');  if ~isfolder(p), mkdir(p); end;  pre = sprintf('sub-%02d_', s);
                t = table(hemi, repmat(3, 12, 1), node, rand(rs, 12, 1), rand(rs, 12, 1), 1 + rand(rs, 12, 1), ...
                          'VariableNames', {'hemi','depth','node_id','anisotropy','normal_share','trace'});
                writetable(t, fullfile(p, [pre 'gaugetensor.csv']));
                writetable(table(string(pre), "omega", 'VariableNames', {'subject','dataset'}), fullfile(p, [pre 'metrics.csv']));
            end
            o = fullfile(d, 'out');
            [C, D] = rheome.group.calibrate(d, o, Centres=Ct, NDraw=4, NPerm=99, NSpin=99, Effects=[0 4], Measures="anisotropy");
            for f = ["draws.csv" "calibration.csv" "tiles.csv" "participants.csv"], tc.verifyTrue(isfile(fullfile(o, f))); end
            tc.verifyEqual(height(D), 4 * 2 * 2);  tc.verifyEqual(height(C), 2 * 2);
            tc.verifyEqual(unique(D.n), 24);  tc.verifyEqual(unique(D.n_tiles), 12);
            tc.verifyEqual(C.detected_rate(C.effect == 4), [1; 1]);   % 4 sd at the seed tile, 12 vs 12
            tc.verifyTrue(all(ismember({'any_fwer_lo','fp_outside_hi','spin_rate','calibrated'}, C.Properties.VariableNames)));
        end
    end
end

% Author: Diellor Basha, 2026
