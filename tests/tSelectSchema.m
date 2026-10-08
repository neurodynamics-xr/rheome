classdef tSelectSchema < matlab.unittest.TestCase
% The schema is consistent, the DDL says it, and every relation materialises from a
% store with the right columns and the diagonal's row counts.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function keysAndForeignKeysAreConsistent(tc)
            S = rheome.select.schema();
            names = {S.name};
            for r = S
                tc.verifyTrue(all(ismember(r.key, r.columns)), r.name);
                fk = fieldnames(r.foreign);
                for i = 1:numel(fk)
                    tc.verifyTrue(ismember(r.foreign.(fk{i}), names), sprintf('%s -> %s', r.name, r.foreign.(fk{i})));
                end
                tc.verifyEqual(numel(r.columns), numel(r.types));
            end
        end

        function ddlNamesEveryRelationKeyAndIndex(tc)
            txt = rheome.select.ddl();
            for r = rheome.select.schema()
                tc.verifyTrue(contains(txt, sprintf('CREATE TABLE %s (', r.name)));
                tc.verifyTrue(contains(txt, sprintf('PRIMARY KEY (%s)', strjoin(r.key, ', '))));
            end
            tc.verifyTrue(contains(txt, 'REFERENCES tile (recording_id, level, k)'));
            tc.verifyTrue(contains(txt, 'CREATE INDEX feature_band_env'));
            tc.verifyTrue(contains(txt, 'CREATE UNIQUE INDEX txn_content'));
        end

        function everyRelationMaterialisesWithItsColumns(tc)
            fx = selFixture();  db = fx.db;  S = rheome.select.schema();
            for r = S
                % the 'note' relations live in the sidecar file and have their own readers
                % (rheome.select.labels, rheome.select.measures); `rows` materialises the store's own.
                if strcmp(r.kind, 'note'), continue; end
                if any(strcmp(r.name, {'tile','tile_cone','feature','feature_band','feature_group','feature_group_band','feature_space','feature_group_space','feature_cortex','feature_cortex_space'}))
                    R = rheome.select.rows(db, r.name, Level=db.grid.Lmax);
                else
                    R = rheome.select.rows(db, r.name);
                end
                tc.verifyEqual(R.Properties.VariableNames, r.columns, r.name);
            end
        end

        function rowCountsFollowTheDiagonal(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            for L = [0 2 g.Lmax]
                K = g.K(L+1);  nb = numel(g.bandsAt{L+1});
                tc.verifyEqual(height(rheome.select.rows(db, "tile", Level=L)), K);
                tc.verifyEqual(height(rheome.select.rows(db, "feature", Level=L)), K * fx.C);
                tc.verifyEqual(height(rheome.select.rows(db, "feature_band", Level=L)), K * fx.C * nb);
                tc.verifyEqual(height(rheome.select.rows(db, "tile_cone", Level=L)), K * nb);
            end
            tc.verifyEqual(height(rheome.select.rows(db, "channel")), fx.C);
            tc.verifyEqual(height(rheome.select.rows(db, "band")), height(db.bands));
        end

        function featureBandRowsEqualTheStoreArrays(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;  L = g.Lmax - 2;
            R = rheome.select.rows(db, "feature_band", Level=L);
            present = g.bandsAt{L+1};
            E = reshape(db.m.(sprintf('L%02d_energy', L)), g.K(L+1), fx.C, numel(present));
            for i = 1:height(R)
                tc.verifyEqual(R.energy(i), E(R.k(i), R.channel_id(i), find(present == R.band_id(i))));
            end
        end

        function readsAreCounted(tc)
            fx = selFixture();  db = fx.db;
            b0 = db.cost('bytes');  r0 = db.cost('reads');
            rheome.select.level(db, 0, 'energy');  rheome.select.level(db, 0, 'energy');       % second hit is cached
            tc.verifyGreaterThan(db.cost('bytes'), b0);
            tc.verifyEqual(db.cost('reads') - r0, 1);
        end

    end
end
