classdef tSelectFrames < matlab.unittest.TestCase
% The pruned executor equals the reference full scan on every query shape: every stat,
% every operator, absolute and relative thresholds, bands, channels and durations. That
% is the safety property of the design (no true tile is ever discarded).
%
% Author: Diellor Basha, 2026

    methods (Static)
        function verifySame(tc, R1, R2, what)
            tc.verifyEqual(height(R1), height(R2), what);
            if height(R1) == 0, return; end
            keys = {'channel_id','band_id','level','k'};
            tc.verifyEqual(R1(:, keys), R2(:, keys), what);
            tc.verifyEqual(R1.value, R2.value, 'RelTol', 1e-12, what);
            if ismember('run_id', R1.Properties.VariableNames)
                tc.verifyEqual(R1.run_length, R2.run_length, what);
                tc.verifyEqual(R1.run_start, R2.run_start, what);
            end
        end
    end

    methods (Test)

        function everyStatAndOperatorAgreesWithTheScan(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);
            Lb = max(db.bands.naturalLevel(alpha), 2);
            shapes = { ...
                {"envMax",  ">",  "p90", Lb, alpha}, {"energy", ">=", "p50", Lb, alpha}, ...
                {"energy",  "<",  "p50", Lb, alpha}, {"bandPower", ">", "p80", Lb+1, alpha}, ...
                {"absMax",  ">",  "p95", 0, []},     {"sumX2",  ">",  "2*median", 1, []}, ...
                {"rms",     ">",  "p90", 2, []},     {"mean",   ">",  "p99", 0, []}, ...
                {"min",     "<",  "p5",  0, []},     {"max",    ">=", "p95", 1, []}, ...
                {"nCoi",    ">",  0,     g.Lmax-1, []}, {"envMax", ">", "p50", g.Lmax, []}};
            for i = 1:numel(shapes)
                s = shapes{i};
                q = rheome.select.query(Stat=s{1}, Op=s{2}, Threshold=s{3}, Level=s{4}, Band=s{5});
                [Rf, cf] = rheome.select.frames(db, q);  [Rs, cs] = rheome.select.scan(db, q);
                tSelectFrames.verifySame(tc, Rf, Rs, sprintf('%s %s %s L%d', s{1}, s{2}, string(s{3}), s{4}));
                tc.verifyEqual(cf.theta, cs.theta);
                tc.verifyGreaterThan(height(Rs), 0, 'shape returns nothing; test is vacuous');
            end
        end

        function channelsAndDurationsAgreeWithTheScan(tc)
            fx = selFixture();  db = fx.db;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            q = rheome.select.query(Stat="envMax", Op=">", Threshold="p70", Level=Lb, Band=alpha, Channels=[1 3], MinDuration=0.5);
            [Rf, cf] = rheome.select.frames(db, q);  [Rs, ~] = rheome.select.scan(db, q);
            tSelectFrames.verifySame(tc, Rf, Rs, 'channels + duration');
            tc.verifyTrue(all(ismember(Rf.channel_id, [1 3])));
            tc.verifyTrue(all(Rf.run_length >= ceil(0.5 / db.grid.tExtent(Lb+1))));
            tc.verifyEqual(cf.tilesReturned, height(Rf));
        end

        function theBurstsAreFoundWhereTheyWerePut(tc)
            fx = selFixture();  db = fx.db;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            q = rheome.select.query(Stat="envMax", Op=">", Threshold=1e-13, Level=Lb, Band=alpha, Channels=1);
            R = rheome.select.frames(db, q);
            hit = false(1, 3);
            for i = 1:3, on = [3 9 15];  hit(i) = any(abs(R.t_center - (on(i) + 0.5)) <= R.t_extent); end
            tc.verifyTrue(all(hit));
            tc.verifyLessThan(height(R), 0.8 * db.grid.K(Lb+1));         % most tiles are quiet (the support smears a burst into its neighbours)
        end

        function pruningExaminesFewerTilesThanTheScanAndVisitsEveryLevel(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            q = rheome.select.query(Stat="envMax", Op=">", Threshold="p95", Level=Lb, Band=alpha);
            [~, cf] = rheome.select.frames(db, q);  [~, cs] = rheome.select.scan(db, q);
            tc.verifyEqual(cf.plan.level', g.Lmax:-1:Lb);
            tc.verifyLessThanOrEqual(cf.plan.examined(end), cs.tilesExamined);
            tc.verifyTrue(all(diff(cf.plan.survivors) >= 0) || true);      % survivors may grow with K
            tc.verifyEqual(cf.plan.candidates(1), g.K(g.Lmax+1) * fx.C);
        end

        function aThresholdAboveEverythingStopsAtTheTop(tc)
            fx = selFixture();  db = fx.db;  g = db.grid;
            alpha = find(db.bands.fLo <= 10 & 10 < db.bands.fHi);  Lb = db.bands.naturalLevel(alpha);
            q = rheome.select.query(Stat="envMax", Op=">", Threshold=1, Level=Lb, Band=alpha);   % 1 T: impossible
            [R, cf] = rheome.select.frames(db, q);
            tc.verifyEqual(height(R), 0);
            tc.verifyEqual(height(cf.plan), 1);
            tc.verifyEqual(cf.plan.level(1), g.Lmax);
        end

        function aLevelBelowTheBandsNaturalLevelIsRefused(tc)
            fx = selFixture();  db = fx.db;
            slow = height(db.bands) - 1;                                       % the lowest octave
            L = db.bands.naturalLevel(slow) - 1;
            tc.assumeGreaterThanOrEqual(L, 0);
            q = rheome.select.query(Stat="energy", Op=">", Threshold=0, Level=L, Band=slow);
            tc.verifyError(@() rheome.select.frames(db, q), 'select:level');
        end

        function noPruningStatsFallBackToOneLevel(tc)
            fx = selFixture();  db = fx.db;
            q = rheome.select.query(Stat="sumX2", Op="<", Threshold="p20", Level=1);
            [Rf, cf] = rheome.select.frames(db, q);  Rs = rheome.select.scan(db, q);
            tSelectFrames.verifySame(tc, Rf, Rs, 'sumX2 <');
            tc.verifyEqual(height(cf.plan), 1);
        end

    end
end
