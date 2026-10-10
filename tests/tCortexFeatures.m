classdef tCortexFeatures < matlab.unittest.TestCase
% The cortical features roll up exactly: a node is the sum of its two children and a tile of its two
% halves, to 1e-12 relative, and rheome.select.rows reads the sidecar in the schema's columns.
% (Exit criterion 2 of docs/2026-09-29-alpha-labelling-design.md, synthetic half.)
%
% Author: Diellor Basha, 2026

    methods (Test)

        function aNodeIsTheSumOfItsChildren(tc)
            [R, ~] = i_synth();
            for L = 1:numel(R.lev)
                X = R.lev{L};
                for id = R.ids(R.ids < 8)'
                    [~, p] = ismember(id, R.ids);  [~, c] = ismember([2*id; 2*id+1], R.ids);
                    tc.verifyLessThan(max(abs(X(p,:,:) - sum(X(c,:,:), 1)), [], 'all') / max(abs(X(p,:,:)), [], 'all'), 1e-12);
                end
            end
        end

        function aTileIsTheSumOfItsHalves(tc)
            [R, K] = i_synth();
            for L = 2:numel(K)
                A = R.lev{L};  B = R.lev{L-1};  kk = 1:K(L);
                lo = 2*kk - 1;  hi = min(2*kk, K(L-1));  odd = hi == lo;
                S = B(:, lo, :) + B(:, hi, :) .* reshape(~odd, 1, [], 1);
                tc.verifyLessThan(max(abs(A - S), [], 'all') / max(abs(A), [], 'all'), 1e-12, sprintf('level %d', L-1));
            end
        end

        function theRootHoldsEverything(tc)
            [R, K, X] = i_synth();
            top = R.lev{end};
            tc.verifyEqual(R.ids(1), 1);  tc.verifyEqual(size(top, 2), K(end));
            tc.verifyEqual(squeeze(sum(top(1,:,:), 2)), squeeze(sum(X, [1 2])), 'RelTol', 1e-12);
        end

        function mixedDepthLeavesAreRefused(tc)
            tc.verifyError(@() rheome.select.rollupcortex(ones(2, 3), [8; 4], [3 2]), 'select:rollupcortex:depth');
        end

        function rowsReadTheSidecarInSchemaColumns(tc)
            fx = selFixture();  db = fx.db;  K = db.grid.K(:)';
            leaves = (8:15)';  rng(3);
            X = rand(numel(leaves), K(1), 5);
            R = rheome.select.rollupcortex(X, leaves, K);
            d = tempname;  mkdir(d);  cl = onCleanup(@() rmdir(d, 's'));
            db.file = fullfile(d, 'store.mat');
            F = struct('ids', R.ids, 'lev', {R.lev}, 'band_id', 6, 'nodes', table());
            save(rheome.select.cortexfile(db), '-struct', 'F', '-v7.3');
            sc = rheome.select.schema();
            for name = ["feature_cortex" "feature_cortex_space"]
                T = rheome.select.rows(db, name, Level=1);
                tc.verifyEqual(T.Properties.VariableNames, sc(strcmp({sc.name}, name)).columns, name);
            end
            T = rheome.select.rows(db, "feature_cortex", Level=1);
            tc.verifyEqual(height(T), numel(R.ids) * K(2));
            r = T(T.node_id == 1 & T.k == 1, :);
            tc.verifyEqual(r.energy, R.lev{2}(1, 1, 3), 'RelTol', 1e-14);
            T = rheome.select.rows(db, "feature_cortex_space", Level=0);
            tc.verifyEqual(unique(T.cband_id)', 1:2);
        end

        function theRealSidecarRollsUpExactly(tc)
            % exit criterion 2 on the real store: depth-3 nodes = their children at level 3, a level-4
            % tile = its two level-3 tiles; needs rheome.flow.cortexfeatures('sub01') to have run
            Ct = rheome.select.catalog();
            tc.assumeFalse(isempty(Ct), 'no tile store in the cache');
            rheomeTestSubject(tc, 'subject');   % skips, with the reason, when the cache or the subject is absent
            hit = find(Ct.recording_id == string(rheomeTestSubject()) & Ct.bank == "frame" & Ct.store == "default", 1);
            tc.assumeNotEmpty(hit, 'test subject store absent');
            db = rheome.select.open(char(Ct.file(hit)));
            tc.assumeTrue(isfile(rheome.select.cortexfile(db)), 'cortical sidecar not built');
            F = load(rheome.select.cortexfile(db), 'ids', 'lev', 'K');
            A = F.lev{4};  B = F.lev{5};
            for id = F.ids(floor(log2(F.ids)) == 4)'
                [~, p] = ismember(id, F.ids);  [~, c] = ismember([2*id; 2*id+1], F.ids);
                tc.verifyLessThan(max(abs(A(p,:,:) - sum(A(c,:,:), 1)), [], 'all') / max(abs(A(p,:,:)), [], 'all'), 1e-12);
            end
            K3 = size(A, 2);  S = A(:, 1:2:end, :);  n = floor(K3/2);
            S(:, 1:n, :) = S(:, 1:n, :) + A(:, 2:2:end, :);
            tc.verifyLessThan(max(abs(B - S), [], 'all') / max(abs(B), [], 'all'), 1e-12);
        end

        function aFragmentDoesNotStopTheTree(tc)
            % another subject's left hemisphere: a one-vertex fragment peeled off at depth 4 used to end the tree there
            rheomeTestSubject(tc, 'second');   % skips, with the reason, when the cache or the subject is absent
            try, B = rheome.load.bases(rheomeTestSubject('second')); catch, tc.assumeFail('second test subject bases absent'); end
            for h = ["L" "R"]
                T = rheome.geom.tree(B.(h).S, L=B.(h).lbo.L, M=B.(h).lbo.Mass, MaxDepth=7);
                tc.verifyEqual(height(T), 255, h);
                tc.verifyGreaterThan(min(T.n_vertices), 1, h);
            end
        end

        function classFollowsAncestryAndWritesIdempotently(tc)
            % leaf 8 on throughout, sibling 9 never: sustained, not widespread -> class 2 everywhere.
            % leaves 10, 11 on for the first second only: brief, their parent (5) full -> class 3.
            fx = selFixture();  db = fx.db;  K = db.grid.K(:)';  n0 = K(1);
            X = zeros(8, n0, 3);  X(:,:,1) = 1;
            X(1,:,2) = 1;  X(3:4,1:4,2) = 1;  X(:,:,3) = X(:,:,2);
            R = rheome.select.rollupcortex(X, (8:15)', K);
            F = struct('ids', R.ids, 'lev', {R.lev}, 'K', K, 'tExtent', db.grid.tExtent, 'band_id', 6, ...
                'leafIds', (8:15)', 'leafBimodal', true(8, 1), 'nodes', i_emptynodes());
            f = rheome.select.cortexfile(db);  save(f, '-struct', 'F', '-v7.3');  cl = onCleanup(@() delete(f));
            L = rheome.flow.labelalpha("x", Families="class", Store=string(db.file), WideDepth=1, SustainS=4, Write=false);
            c = L.rows.alpha_class;  o = L.rows.alpha_onFrac;
            tc.verifyTrue(all(c.value(c.unit_id == 8) == 2));
            tc.verifyEqual(c.value(c.unit_id == 10 & c.k == 1), 3);
            tc.verifyEmpty(c.value(c.unit_id == 9));
            tc.verifyEqual(o.value(o.unit_id == 10 & o.k == 1), 0.5, 'AbsTol', 1e-12);
            m0 = height(rheome.select.measures(db));
            L1 = rheome.flow.labelalpha("x", Families="class", Store=string(db.file), WideDepth=1, SustainS=4);
            L2 = rheome.flow.labelalpha("x", Families="class", Store=string(db.file), WideDepth=1, SustainS=4);
            tc.verifyFalse(L1.txn.alpha_class.existed);  tc.verifyTrue(L2.txn.alpha_class.existed);
            tc.verifyEqual(L2.txn.alpha_class.txn_id, L1.txn.alpha_class.txn_id);
            % the query: leaf 8 (left) is class 2 in every 2 s tile -> one state, holding one episode
            ids = [];
            for nm = ["tOn" "tOff" "tauR" "tauF" "plateauS" "r2"]
                t = rheome.select.measure(db, table(2, 0, 3, 0.6, 'VariableNames', {'unit_id','level','k','value'}), ...
                    Kind="alpha_episode_" + nm, Scope="cortex", Source="rheome.flow.labelalpha");
                ids(end+1) = t.txn_id; %#ok<AGROW>
            end
            Q = rheome.flow.alphastates("x", Class=2, Hemi="L", Store=string(db.file), Export=false);
            K3 = db.grid.K(4);
            tc.verifyEqual(height(Q.states), 1);  tc.verifyEqual(Q.states.node_id, 8);
            tc.verifyEqual([Q.states.k_start Q.states.k_end], [1 K3]);
            tc.verifyEqual(Q.states.n_episodes, 1);  tc.verifyEqual(Q.episodes.t_on, 0.6);
            tc.verifyEqual(Q.union.dur_s, K3 * db.grid.tExtent(4));
            for i = ids, rheome.select.rollback(db, i); end
            for k = ["alpha_class" "alpha_onFrac"], rheome.select.rollback(db, L1.txn.(k).txn_id); end
            tc.verifyEqual(height(rheome.select.measures(db)), m0);
        end

        function rowsWithoutASidecarAreEmptyAndTyped(tc)
            fx = selFixture();  db = fx.db;  db.file = fullfile(tempname, 'none.mat');
            T = rheome.select.rows(db, "feature_cortex", Level=0);
            tc.verifyEqual(height(T), 0);  tc.verifyClass(T.recording_id, 'string');
        end
    end
end

function T = i_emptynodes()
    sc = rheome.select.schema();  f = sc(strcmp({sc.name}, 'cortex_node'));
    T = table('Size', [0 numel(f.columns)], 'VariableTypes', repmat({'double'}, 1, numel(f.columns)), 'VariableNames', f.columns);
end

function [R, K, X] = i_synth()
% depth-3 whole-cortex leaves (ids 8..15), an odd level-0 count so the last coarse tile is partial
    K = [75 38 19 10 5 3 2 1];  rng(7);
    X = rand(8, K(1), 4) .* 10.^(3 * rand(8, 1, 4));
    R = rheome.select.rollupcortex(X, (8:15)', K);
end

% Author: Diellor Basha, 2026
