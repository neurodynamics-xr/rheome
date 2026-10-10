classdef tFlowSweep < matlab.unittest.TestCase
% rheome.flow.sweep: the melt to the measurement matrix, and the round-trip through the store.
%
% ⚠ The sweep tests run with Write=false so the suite does not append to the shared store. The
% write path is tested separately on three rows, and rolled back.
%
% Author: Diellor Basha, 2026
    properties
        M; R; T; db
    end
    methods (TestClassSetup)
        function build(t)
            C = [];
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try, C = rheome.select.catalog(); catch, end
            if isempty(C) || ~any(C.recording_id == string(rheomeTestSubject()) & C.bank == "frame")
                t.assumeFail('no default frame store for test subject');
            end
            hit = find(C.recording_id == string(rheomeTestSubject()) & C.bank == "frame" & C.store == "default", 1);
            t.db = rheome.select.open(char(C.file(hit)));
            [t.M, t.R, t.T] = rheome.flow.sweep(rheomeTestSubject(), Level=3, Depth=2, BandIds=6, ...
                MaxWindows=4, Write=false, Verbose=false);
        end
    end
    methods (Test)
        function everyRowNamesATile(t)
            % ⚠ not every row names a NODE: the global eigenmode rows are scope="none", unit_id 0,
            % because a mode's support is the whole hemisphere. Every row does name a tile.
            t.verifyTrue(all(ismember(t.M.scope, ["cortex","none"])));
            t.verifyTrue(all(t.M.level == 3));
            t.verifyTrue(all(t.M.k >= 1 & t.M.k <= t.db.grid.K(4)));
            c = t.M.scope == "cortex";
            t.verifyTrue(all(t.M.unit_id(c) >= 1 & t.M.unit_id(c) == fix(t.M.unit_id(c))));
            t.verifyTrue(all(t.M.unit_id(~c) == 0));
        end
        function theWindowsAreTheStoresOwnTiles(t)
            % ⭐ the whole point: no second definition of "window". tStart/tEnd must be tile spans.
            g = t.db.grid;
            for i = 1:height(t.T)
                k = t.T.k(i);
                t.verifyEqual(t.T.tStart(i), g.tCenter{4}(k) - g.tExtent(4)/2, 'AbsTol', 1e-9);
            end
        end
        function bandIdIsTheStoresBandId(t)
            t.verifyTrue(all(t.M.band_id == 6));
            t.verifyEqual(unique(t.T.fLo), t.db.bands.fLo(6));
            t.verifyEqual(unique(t.T.fHi), t.db.bands.fHi(6));
        end
        function theMeltCoversEveryKindAndRow(t)
            % ⚠ M holds TWO families with different row counts: the per-node cortical rows and the
            % GLOBAL eigenmode-octave rows (scope="none", unit_id 0, one per tile per octave). An
            % assertion of height(M) == height(T)*nKinds was correct until the second family was
            % added and then failed 192 vs 312, which is the test doing its job.
            c = t.M(t.M.scope == "cortex", :);
            t.verifyEqual(height(c), height(t.T)*numel(unique(c.kind)));
            g = t.M(t.M.scope == "none", :);
            t.verifyTrue(all(g.unit_id == 0), 'a global row must carry unit_id 0');
            t.verifyTrue(~isempty(g), 'the global eigenmode rows are missing');
            key = unique(t.M(:, {'k','unit_id','band_id','kind'}));
            t.verifyEqual(height(key), height(t.M), 'the measurement key is not unique');
        end
        function theGlobalOctaveSharesPartitionTheSpectrum(t)
            g = t.M(startsWith(t.M.kind, "spaceOctaveShare"), :);
            t.assumeNotEmpty(g);
            u = unstack(g(:, {'k','kind','value'}), 'value', 'kind');
            s = sum(u{:, 2:end}, 2);
            t.verifyEqual(s, ones(size(s)), 'AbsTol', 1e-9, ...
                'the spatial octaves must partition the spectrum');
        end
        function writeAndReadBackRoundTrips(t)
            rows = table([3;3;3], [1;2;3], [2;2;2], [6;6;6], [1.5;2.5;3.5], ...
                'VariableNames', {'level','k','unit_id','band_id','value'});
            tx = rheome.select.measure(t.db, rows, Kind="tFlowSweepProbe", Scope="cortex", Source="test");
            cleanup = onCleanup(@() rheome.select.rollback(t.db, tx.txn_id));
            got = rheome.select.measures(t.db, Kind="tFlowSweepProbe", Scope="cortex");
            t.verifyEqual(height(got), 3);
            t.verifyEqual(sort(got.value), [1.5;2.5;3.5], 'AbsTol', 1e-12);
            t.verifyTrue(all(got.unit_id == 2));
            % ⭐ idempotent by content hash: the same rows again must not write twice
            tx2 = rheome.select.measure(t.db, rows, Kind="tFlowSweepProbe", Scope="cortex", Source="test");
            t.verifyEqual(tx2.txn_id, tx.txn_id);
            t.verifyTrue(tx2.existed);
        end
        function aBadCortexUnitIsRejected(t)
            rows = table(3, 1, 0.5, -4, 6, 'VariableNames', {'level','k','value','unit_id','band_id'});
            t.verifyError(@() rheome.select.measure(t.db, rows, Kind="tBad", Scope="cortex"), ...
                'select:measure:unit');
        end
        function theSensorSideCarriesNoInvalidGroupStatistic(t)
            % ⚠ crest/kurtosis/skewness divide a max-over-channels by a sum-over-channels and read
            % crest 0.443-3.027 in group scope, where absMax/rms is >= 1 by definition.
            bad = ["crest","kurtosis","skewness","shapeFactor","impulseFactor","clearanceFactor"];
            t.verifyEmpty(intersect(string(t.R.Properties.VariableNames), bad));
        end
    end
end
