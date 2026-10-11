classdef tUpstreamDriver < matlab.unittest.TestCase
% rheome.scale.measure_upstream -- the AGE-1 participant readout from sensors -- on a synthetic cached subject
% (scaleSynthSubject): a 4 Hz frontal and posterior generator per hemisphere with a lag and no travel, so the
% lagged-source surrogate (rheome.forward.laggedrest) must read what the data read; the periodic gain switched
% off (U_total) reads U when there is no background; a short span is not read; and the passes, halves and
% upstream.csv come out of rheome.scale.run.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_upstream'
        Root
        G1
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(tc.Root, tc.Name));
            % the gradient: z on each hemisphere's registration sphere (low below, high above)
            [V0, ~] = rheome.geom.icosphere(3);
            L = struct('sphere', 100 * V0, 'g', V0(:, 3));  R = L;  %#ok<NASGU>
            tc.G1 = fullfile(tc.Root, 'g1.mat');  save(tc.G1, 'L', 'R');
            old1 = getenv('RHEOME_G1');  setenv('RHEOME_G1', tc.G1);  tc.addTeardown(@setenv, 'RHEOME_G1', old1);
            % Desikan-Killiany names for the G1 sign check: the upper half of each hemisphere is "association"
            a = load(fullfile(tc.Root, tc.Name, 'atlas.mat'));
            a.atlas.Label = {'inferiorparietal L', 'pericalcarine L', 'inferiorparietal R', 'pericalcarine R'};
            save(fullfile(tc.Root, tc.Name, 'atlas.mat'), '-struct', 'a');
            tc.plant(0.02);
        end
    end

    methods
        function plant(tc, noise)
        % 60 s: per hemisphere the most frontal and the most posterior depth-3 tile carry a uniform normal current
        % at 4 Hz, the posterior one 0.8 rad behind (a lagged pair with no travel, shaped as the null's candidates)
            f = fullfile(tc.Root, tc.Name, 'study.mat');  st = load(f);  B = rheome.load.bases(tc.Name);
            sf = load(fullfile(tc.Root, tc.Name, 'surface.mat'));
            fs = st.rec.sfreq;  N = 60 * fs;  t = (0:N-1) / fs;  rng(5);
            J = zeros(size(st.hm.Gain, 2), 2);  Nr = sf.S.VertNormals;
            for h = ["L" "R"]
                Hh = B.(h);  Tr = rheome.geom.tree(Hh.S, L=Hh.lbo.L, M=Hh.lbo.M, MaxDepth=3);  G = rheome.geom.tiles(Tr, Hh.S, 3);
                xc = G.centroid(:, 1);  [~, a] = max(xc);  [~, p] = min(xc);  ar = full(sum(Hh.lbo.Mass, 2));
                for k = [a p]
                    v = Hh.gv(G.tileOf == k);  w = ar(G.tileOf == k);
                    for c = 1:3, J(3*v - 3 + c, 1 + (k == p)) = w .* Nr(v, c); end
                end
            end
            F = st.hm.Gain * J * [sin(2*pi*4*t); 0.8 * sin(2*pi*4*t - 0.8)];
            F = F + noise * std(F(:)) * randn(size(F));
            st.rec.F = F;  st.rec.Time = t;  save(f, '-struct', 'st');
        end
    end

    methods (Test)
        function totalReadsUWithoutBackground(tc)
            % no aperiodic background: the periodic gain is ~1 where the power is, so h = 1 changes nothing
            [T, X] = rheome.scale.measure_upstream(tc.Name, [], Passes=["periodic" "total"], Halves=false, MinSpanS=10, OwnMin=0);
            v = @(m) T.value(T.metric == m);
            tc.verifyEqual(v("read"), 1);
            tc.verifyEqual(v("U_total"), v("U"), 'AbsTol', 0.01);
            tc.verifyEqual(unique(X.pass)', ["periodic" "total"]);
            tc.verifyEqual(v("aperiodic_ok"), 1);
        end

        function surrogateOfALaggedPairReadsTheData(tc)
            % the data ARE a lagged pair that does not travel: the null must read the same share, U* ~ 0
            [T, X] = rheome.scale.measure_upstream(tc.Name, [], Passes=["periodic" "surrogate"], MinSpanS=10, OwnMin=0);
            v = @(m) T.value(T.metric == m);
            tc.verifyNotEqual(v("U"), 0.5);                                     % the lag alone gives a direction
            tc.verifyEqual(v("U_star"), 0, 'AbsTol', 0.1);
            tc.verifyEqual(v("U_star"), v("U") - v("U_surrogate"), 'AbsTol', 1e-12);
            tc.verifyTrue(all(ismember(["U_h1" "U_h2" "U_surrogate_h1" "U_star_h2"], T.metric)));
            tc.verifyEqual(sort(unique(X.span))', ["h1" "h2" "whole"]);
        end

        function shortSpanIsNotRead(tc)
            [T, X] = rheome.scale.measure_upstream(tc.Name, [], MinSpanS=1000, OwnMin=0);
            tc.verifyEqual(T.value(T.metric == "read"), 0);
            tc.verifyFalse(any(startsWith(T.metric, "U")));
            tc.verifySubstring(char(X.reason), 'clean span');
        end

        function halfWithNoCleanSegmentIsNotRead(tc)
            % OMEGA sub-CONP0136 (job 64400011_387): artefact blocks every 3 s through the second half leave h2 no
            % clean Welch segment; that half is recorded as not read and the whole span and h1 are still read
            f = fullfile(tc.Root, tc.Name, 'study.mat');  copyfile(f, [f '.bak']);  tc.addTeardown(@movefile, [f '.bak'], f);
            st = load(f);  fs = st.rec.sfreq;  s = std(st.rec.F(:));
            for t0 = 31:3:48, b = round(t0*fs) + (1:round(0.25*fs));  st.rec.F(:, b) = st.rec.F(:, b) + 50 * s; end
            save(f, '-struct', 'st');
            [T, X] = rheome.scale.measure_upstream(tc.Name, [], Passes="periodic", MinSpanS=10, OwnMin=0);
            tc.verifyEqual(T.value(T.metric == "read"), 1);
            tc.verifyTrue(all(ismember(["U" "U_h1"], T.metric)));
            tc.verifyFalse(ismember("U_h2", T.metric));
            tc.verifySubstring(char(X.reason(X.span == "h2")), 'Welch segments');
        end

        function runWritesUpstreamCsv(tc)
            % through the per-subject driver with the frozen defaults: 60 s is under the 180 s minimum, so the
            % recording is recorded as not read, and the analysis itself succeeds
            od = fullfile(tc.Root, 'out');
            R = rheome.scale.run(tc.Name, Analyses="upstream", OutDir=od);
            tc.verifyEqual(R.timing.status, "ok");
            U = readtable(fullfile(od, tc.Name, 'upstream.csv'), 'TextType', 'string');
            tc.verifySubstring(char(U.reason), 'clean span 60');
            tc.verifyEqual(R.metrics.value(R.metrics.metric == "g1_tilecheck_L"), 1);
        end

        function theGradientSignIsFixedByAnatomy(tc)
            % a G1 map with the opposite sign reads the same U: the sign is set by the regions, not the file
            [T1] = rheome.scale.measure_upstream(tc.Name, [], Passes="periodic", Halves=false, MinSpanS=10, OwnMin=0);
            L = load(tc.G1);  L.L.g = -L.L.g;  L.R.g = -L.R.g;  f = fullfile(tc.Root, 'g1neg.mat');  save(f, '-struct', 'L');
            [T2] = rheome.scale.measure_upstream(tc.Name, [], Passes="periodic", Halves=false, MinSpanS=10, OwnMin=0, Gradient=f);
            tc.verifyEqual(T1.value(T1.metric == "read"), 1);
            tc.verifyEqual(T2.value(T2.metric == "U"), T1.value(T1.metric == "U"), 'AbsTol', 1e-12);
        end

        function everyPassThroughTheDriver(tc)
            [T, X] = rheome.scale.measure_upstream(tc.Name, [], MinSpanS=10, OwnMin=0, OfChunk=200);
            for m = ["U" "U_surrogate" "U_total" "U_of" "U_star" "U_star_h1" "U_star_h2"]
                u = T.value(T.metric == m);
                tc.verifyTrue(isscalar(u) && u >= -1 && u <= 1, m);
            end
            tc.verifyEqual(sort(unique(X.pass))', ["of" "periodic" "surrogate" "total"]);
        end

        function crosstalkReadsAllThreeComponents(tc)
            rng(1);  nV = 30;  K = randn(3*nV, 12);  G = randn(12, 3*nV);  M = sparse(1, 1:10, true, 1, nV);
            own = rheome.flow.crosstalk(K, G, M);
            p = sum((K(1:30, :) * G).^2, 1);  pv = p(1:3:end) + p(2:3:end) + p(3:3:end);
            tc.verifyEqual(own, sum(pv(1:10)) / sum(pv), 'RelTol', 1e-12);
        end

        function gateGReadsReliabilityAndSpreadOnly(tc)
            % three read recordings (+1 not read) with known halves: ICC(2,1) by hand, and no mean in the output
            d = fullfile(tc.Root, 'gateg');  mkdir(d);
            h = [0.10 0.12; -0.05 -0.02; 0.30 0.25];
            for i = 1:3
                M = rheome.scale.rows("upstream", ["read" "U_star_h1" "U_star_h2" "U_star" "aperiodic_ok"], [1 h(i,:) mean(h(i,:)) 1], "share");
                writetable(M, fullfile(d, sprintf('sub-%d_metrics.csv', i)));
            end
            writetable(rheome.scale.rows("upstream", "read", 0, "flag"), fullfile(d, 'sub-9_metrics.csv'));
            G = rheome.scale.groupupstream(d);
            n = 3;  k = 2;  m = mean(h, 'all');  msr = k * sum((mean(h, 2) - m).^2) / (n-1);
            msc = n * sum((mean(h, 1) - m).^2) / (k-1);  E = h - mean(h, 2) - mean(h, 1) + m;  mse = sum(E.^2, 'all') / ((n-1)*(k-1));
            tc.verifyEqual(G.icc, (msr - mse) / (msr + (k-1)*mse + k*(msc - mse)/n), 'AbsTol', 1e-12);
            tc.verifyEqual(G.sd, std(mean(h, 2)), 'AbsTol', 1e-12);
            tc.verifyEqual([G.n_read G.n_not_read G.n_aperiodic_failed], [3 1 0]);
            tc.verifyFalse(any(contains(lower(fieldnames(G)), "mean")));
        end

        function laggedRestTilesTheSpan(tc)
            rng(2);  fs = 50;  t = (0:fs*30-1) / fs;  Gp = randn(20, 4);  front = logical([1 0 1 0]);
            Y = Gp(:, 1) * sin(2*pi*3*t) + 0.6 * Gp(:, 2) * sin(2*pi*3*t - 1) + 0.01 * randn(20, numel(t));
            [Ys, info] = rheome.forward.laggedrest(Y, Gp, front, fs);
            tc.verifyTrue(all(info.read));
            W = info.windows;  tc.verifyEqual(W.first(2:end), W.last(1:end-1) + 1);
            tc.verifyEqual([W.first(1) W.last(end)], [1 numel(t)]);
            tc.verifyEqual(median(W.lagS), 1 / (2*pi*3), 'AbsTol', 0.021);     % the planted lag, 1 rad at 3 Hz
            tc.verifyEqual(size(Ys), size(Y));
            mk = false(1, numel(t));  mk(700) = true;
            [Ym, im] = rheome.forward.laggedrest(Y, Gp, front, fs, Mask=mk);
            tc.verifyFalse(im.read(700));  tc.verifyEqual(Ym(:, 700), zeros(20, 1));
        end
    end
end
