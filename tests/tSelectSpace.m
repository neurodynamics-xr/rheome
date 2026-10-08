classdef tSelectSpace < matlab.unittest.TestCase
% The wavelength axis in the query side: the space relations materialise, spatial stats
% agree between executors in both scopes, and the SQL names the space relations.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function spaceRelationsMaterialise(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            S = rheome.select.rows(db, "space_band");
            tc.verifyGreaterThan(height(S), 2);
            tc.verifyEqual(S.kind(1), "longer");  tc.verifyEqual(S.kind(end), "shorter");
            L = g.Lmax - 1;
            R = rheome.select.rows(db, "feature_space", Level=L);
            tc.verifyEqual(height(R), g.K(L+1) * fx.C * height(S));
            Rg = rheome.select.rows(db, "feature_group_space", Level=L);
            tc.verifyEqual(height(Rg), g.K(L+1) * numel(db.groupNodes) * height(S));
            % root partition through the relations
            F = rheome.select.rows(db, "feature", Level=L);
            root = Rg(Rg.group_id == 1, :);
            tc.verifyEqual(accumarray(root.k, root.energy), accumarray(F.k, F.sum_x2), 'RelTol', 1e-10);
        end

        function spatialStatsAgreeBetweenExecutors(tc)
            fx = selFixture();  db = fx.db;
            nS = height(db.sbands);
            shapes = {{"spaceEnergy", ">", "p70", 2, [], "channel", []}, {"spaceEnvMax", ">=", "p50", 1, 1, "channel", []}, ...
                      {"spaceEnergy", ">", "p60", 2, nS, "group", 1}, {"spaceEnergy", "<", "p40", 3, [], "group", []}};
            for i = 1:numel(shapes)
                s = shapes{i};
                q = rheome.select.query(Stat=s{1}, Op=s{2}, Threshold=s{3}, Level=s{4}, Band=s{5}, Scope=s{6}, Depth=s{7});
                [Rf, cf] = rheome.select.frames(db, q);  Rs = rheome.select.scan(db, q);
                tc.verifyEqual(height(Rf), height(Rs), s{1});
                if height(Rf) > 0
                    tc.verifyEqual(Rf.value, Rs.value, 'RelTol', 1e-12);
                    tc.verifyEqual(Rf.band_id, Rs.band_id);
                end
                tc.verifyGreaterThan(height(Rs), 0, 'vacuous');
            end
        end

        function theSqlNamesTheSpaceRelations(tc)
            q = rheome.select.query(Stat="spaceEnergy", Op=">", Threshold="p90", Level=2, Band=2);
            txt = rheome.select.sql(q);
            tc.verifyTrue(contains(txt, 'FROM feature_space'));
            tc.verifyTrue(contains(txt, 'sband_id IN (2)'));
            tc.verifyTrue(contains(txt, 'channel_id, sband_id, k'));
            qg = rheome.select.query(Stat="spaceEnvMax", Op=">", Threshold=1e-13, Level=2, Scope="group", Depth=1, MinDuration=1);
            txt = rheome.select.sql(qg);
            tc.verifyTrue(contains(txt, 'FROM feature_group_space'));
            tc.verifyTrue(contains(txt, 'PARTITION BY group_id, sband_id ORDER BY k'));
            tc.verifyTrue(endsWith(strtrim(txt), 'ORDER BY group_id, sband_id, k;'));
            tc.verifyTrue(contains(rheome.select.ddl(), 'CREATE TABLE space_band ('));
        end

        function aStoreWithoutTheAxisSaysSo(tc)
            fx = selFixture();  db = fx.db;
            db.sbands = table();
            q = rheome.select.query(Stat="spaceEnergy", Op=">", Threshold=0, Level=2);
            tc.verifyError(@() rheome.select.scan(db, q), 'select:space');
        end

    end
end
