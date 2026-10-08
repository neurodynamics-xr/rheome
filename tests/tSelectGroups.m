classdef tSelectGroups < matlab.unittest.TestCase
% The position pyramid in the query side: group relations materialise, group scope agrees
% between executors, the tree descent equals the scan at its depth, and the SQL names
% the group relations.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function groupRelationsMaterialiseWithExactSums(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            S = rheome.select.rows(db, "sensor_group");
            tc.verifyEqual(height(S), 2 * fx.C - 1);
            tc.verifyEqual(S.parent_id(1), 0);
            L = g.Lmax - 1;
            G = rheome.select.rows(db, "feature_group", Level=L);  F = rheome.select.rows(db, "feature", Level=L);
            tc.verifyEqual(height(G), g.K(L+1) * numel(db.groupNodes));
            root = G(G.group_id == 1, :);
            tc.verifyEqual(root.sum_x2, accumarray(F.k, F.sum_x2), 'RelTol', 1e-12);
            GB = rheome.select.rows(db, "feature_group_band", Level=L);
            tc.verifyEqual(GB.Properties.VariableNames, {'recording_id','group_id','level','k','band_id','energy','env_max'});
        end

        function groupScopeAgreesBetweenExecutors(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            shapes = {{"envMax", ">", "p80", Lb, alpha, []}, {"energy", ">=", "p50", Lb+1, alpha, 1}, ...
                      {"absMax", ">", "p90", 0, [], []}, {"rms", ">", "p90", 1, [], 0}, {"sumX2", "<", "p30", 1, [], []}};
            for i = 1:numel(shapes)
                s = shapes{i};
                q = rheome.select.query(Stat=s{1}, Op=s{2}, Threshold=s{3}, Level=s{4}, Band=s{5}, Scope="group", Depth=s{6});
                [Rf, cf] = rheome.select.frames(db, q);  Rs = rheome.select.scan(db, q);
                tc.verifyEqual(height(Rf), height(Rs), sprintf('%s', s{1}));
                if height(Rf) > 0
                    tc.verifyEqual(Rf(:, {'group_id','band_id','k'}), Rs(:, {'group_id','band_id','k'}));
                    tc.verifyEqual(Rf.value, Rs.value, 'RelTol', 1e-12);
                end
                tc.verifyTrue(ismember('group_id', Rf.Properties.VariableNames));
                tc.verifyGreaterThan(height(Rs), 0, 'vacuous shape');
                tc.verifyEqual(cf.executor, 'frames');
            end
        end

        function theTreeDescentEqualsTheScanAtItsDepth(tc)
            fx = selFixture();  db = fx.db;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            for d = [0 1]
                q = rheome.select.query(Stat="envMax", Op=">", Threshold="p60", Level=Lb, Band=alpha, Scope="group", Depth=d);
                [Rt, ct] = rheome.select.tree(db, q);  Rs = rheome.select.scan(db, q);
                tc.verifyEqual(height(Rt), height(Rs), sprintf('depth %d', d));
                if height(Rt) > 0, tc.verifyEqual(Rt(:, {'group_id','band_id','k'}), Rs(:, {'group_id','band_id','k'})); end
                tc.verifyEqual(ct.executor, 'tree');
                tc.verifyEqual(ct.plan.depth', 0:d);
            end
            q = rheome.select.query(Stat="energy", Op=">", Threshold=1, Level=Lb, Band=alpha, Scope="group", Depth=1);   % impossible
            [Rt, ct] = rheome.select.tree(db, q);
            tc.verifyEqual(height(Rt), 0);  tc.verifyEqual(height(ct.plan), 1);
        end

        function theRootIsTheWholeArray(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            q = rheome.select.query(Stat="sumX2", Op=">", Threshold=0, Level=g.Lmax, Scope="group", Groups=1);
            R = rheome.select.scan(db, q);
            tc.verifyEqual(height(R), 1);
            tc.verifyEqual(R.value, sum(fx.X.^2, 'all'), 'RelTol', 1e-10);
        end

        function sqlNamesTheGroupRelations(tc)
            q = rheome.select.query(Stat="envMax", Op=">", Threshold="p95", Level=3, Band=5, Scope="group", Depth=2);
            txt = rheome.select.sql(q);
            tc.verifyTrue(contains(txt, 'FROM feature_group_band'));
            tc.verifyTrue(contains(txt, 'group_id IN (SELECT group_id FROM sensor_group WHERE recording_id = :rec AND depth = 2)'));
            tc.verifyTrue(contains(txt, 'p.group_id = f.group_id'));
            tc.verifyTrue(endsWith(strtrim(txt), 'ORDER BY group_id, band_id, k;'));
            txt = rheome.select.ddl();
            tc.verifyTrue(contains(txt, 'CREATE TABLE sensor_group ('));
            tc.verifyTrue(contains(txt, 'REFERENCES sensor_group (recording_id, group_id)'));
        end

        function channelsAndGroupScopeAreExclusive(tc)
            tc.verifyError(@() rheome.select.query(Stat="energy", Op=">", Threshold=0, Level=3, Band=1, Scope="group", Channels=1), 'select:query:scope');
        end

    end
end
