classdef tSelectSql < matlab.unittest.TestCase
% The SQL rendering names the right relations, keys and semi-joins, and is pinned as
% text for one fixed query so a change is a deliberate change.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function aBandQueryRendersTheSemiJoinChain(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            q = rheome.select.query(Stat="envMax", Op=">", Threshold="p95", Level=2, Band=5);
            txt = rheome.select.sql(q, db);
            tc.verifyTrue(contains(txt, 'percentile_cont(0.95) WITHIN GROUP (ORDER BY env_max)'));
            tc.verifyTrue(contains(txt, sprintf('s%d AS (SELECT channel_id, band_id, k FROM feature_band WHERE recording_id = ''%s''', g.Lmax, db.recording_id)));
            tc.verifyTrue(contains(txt, 'p.k = (f.k - 1) / 2 + 1'));
            tc.verifyTrue(contains(txt, sprintf('s%d p ON', 3)));                                % the last semi-join is to level 3
            tc.verifyTrue(contains(txt, 'JOIN tile t ON t.recording_id = f.recording_id AND t.level = f.level AND t.k = f.k'));
            tc.verifyTrue(contains(txt, 'f.env_max > (SELECT theta FROM thr)'));
            tc.verifyTrue(endsWith(strtrim(txt), 'ORDER BY channel_id, band_id, k;'));
            tc.verifyEqual(count(txt, ' AS (SELECT'), 1 + (g.Lmax - 2) + 1);                 % thr + survivors + hits
        end

        function aDerivedStatJoinsTileAndUsesTheProxyBound(tc)
            q = rheome.select.query(Stat="rms", Op=">", Threshold=2e-13, Level=1);
            txt = rheome.select.sql(q);
            tc.verifyTrue(contains(txt, 'sqrt(f.sum_x2 / NULLIF(t.n, 0)) > (SELECT theta FROM thr)'));
            tc.verifyTrue(contains(txt, 'abs_max >= (SELECT theta FROM thr)'));            % the bound
            tc.verifyTrue(contains(txt, 'thr AS (SELECT 2e-13 AS theta)'));
        end

        function durationsRenderGapsAndIslands(tc)
            fx = selFixture();  db = fx.db;
            q = rheome.select.query(Stat="energy", Op=">", Threshold="p80", Level=3, Band=5, MinDuration=1.5);
            txt = rheome.select.sql(q, db);
            tc.verifyTrue(contains(txt, 'ROW_NUMBER() OVER (PARTITION BY channel_id, band_id ORDER BY k)'));
            tc.verifyTrue(contains(txt, 'WHERE run_length >= CEIL(1.5 / 2)'));
        end

        function aNoPruningQueryHasNoChain(tc)
            q = rheome.select.query(Stat="energy", Op="<", Threshold="p20", Level=3, Band=5);
            txt = rheome.select.sql(q);
            tc.verifyFalse(contains(txt, 'p.k = (f.k - 1) / 2 + 1'));
            tc.verifyTrue(contains(txt, 'f.energy < (SELECT theta FROM thr)'));
        end

        function thePinnedQueryIsStable(tc)
            q = rheome.select.query(Stat="envMax", Op=">", Threshold=3e-13, Level=1, Band=4, Channels=[1 2]);
            txt = rheome.select.sql(q);
            expected = sprintf(['WITH thr AS (SELECT 3e-13 AS theta),\n' ...
                's_top AS (SELECT channel_id, band_id, k FROM feature_band WHERE recording_id = :rec AND band_id IN (4) AND channel_id IN (1, 2) AND level = :top AND env_max >= (SELECT theta FROM thr)),\n' ...
                '-- repeat for each level from :top - 1 down to 2: s_L AS (SELECT f.channel_id, f.band_id, f.k FROM feature_band f JOIN s_(L+1) p ON <parent key> WHERE f.level = L AND env_max >= (SELECT theta FROM thr)),\n' ...
                'hits AS (SELECT f.recording_id, f.channel_id, f.level, f.k, t.t_center, t.t_extent, f.band_id, ''envMax'' AS stat, f.env_max AS value FROM feature_band f JOIN tile t ON t.recording_id = f.recording_id AND t.level = f.level AND t.k = f.k JOIN s_top p ON p.channel_id = f.channel_id AND p.band_id = f.band_id AND p.k = (f.k - 1) / 2 + 1 WHERE f.recording_id = :rec AND f.band_id IN (4) AND f.channel_id IN (1, 2) AND f.level = 1 AND f.env_max > (SELECT theta FROM thr))\n' ...
                'SELECT * FROM hits ORDER BY channel_id, band_id, k;']);
            tc.verifyEqual(txt, expected);
        end

    end
end
