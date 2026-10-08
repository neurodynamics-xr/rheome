classdef tSelectLabel < matlab.unittest.TestCase
% Label transactions: atomic, idempotent, readable by kind, rolled back by id.
%
% Author: Diellor Basha, 2026

    methods (TestMethodSetup)
        function clean(tc)
            fx = selFixture();
            if exist(fx.db.labelFile, 'file') == 2, delete(fx.db.labelFile); end
        end
    end

    methods (Test)

        function aTransactionWritesBothTablesTogether(tc)
            fx = selFixture();  db = fx.db;
            rows = table([1; 1; 2], [3; 3; 3], [2; 3; 5], 'VariableNames', {'channel_id','level','k'});
            t = rheome.select.label(db, rows, Kind="bad", Source="test");
            tc.verifyFalse(t.existed);
            tc.verifyEqual(t.count, 3);
            [L, T] = rheome.select.labels(db);
            tc.verifyEqual(height(L), 3);  tc.verifyEqual(height(T), 1);
            tc.verifyTrue(all(L.txn_id == t.txn_id));
            tc.verifyEqual(L.kind, repmat("bad", 3, 1));
            tc.verifyEqual(exist([db.labelFile '.tmp'], 'file'), 0);            % no temp file left behind
        end

        function theSameRowsAreIdempotent(tc)
            fx = selFixture();  db = fx.db;
            rows = table([1; 2], [3; 3], [2; 5], 'VariableNames', {'channel_id','level','k'});
            t1 = rheome.select.label(db, rows, Kind="bad");
            t2 = rheome.select.label(db, rows([2 1], :), Kind="bad");                 % other order, same content
            tc.verifyTrue(t2.existed);
            tc.verifyEqual(t2.txn_id, t1.txn_id);
            tc.verifyEqual(height(rheome.select.labels(db)), 2);
            t3 = rheome.select.label(db, rows, Kind="artefact");                       % other kind: a new transaction
            tc.verifyFalse(t3.existed);
            tc.verifyEqual(height(rheome.select.labels(db)), 4);
        end

        function labelsReadBackByKindLevelAndChannel(tc)
            fx = selFixture();  db = fx.db;
            rheome.select.label(db, table([1; 2], [3; 4], [2; 1], 'VariableNames', {'channel_id','level','k'}), Kind="bad");
            rheome.select.label(db, table(0, 2, 7, 'VariableNames', {'channel_id','level','k'}), Kind="note");
            tc.verifyEqual(height(rheome.select.labels(db, Kind="bad")), 2);
            tc.verifyEqual(height(rheome.select.labels(db, Kind="bad", Level=4)), 1);
            tc.verifyEqual(height(rheome.select.labels(db, Channels=0)), 1);
        end

        function rollbackRestoresThePriorState(tc)
            fx = selFixture();  db = fx.db;
            t1 = rheome.select.label(db, table(1, 3, 2, 'VariableNames', {'channel_id','level','k'}), Kind="bad");
            t2 = rheome.select.label(db, table(1, 3, 4, 'VariableNames', {'channel_id','level','k'}), Kind="bad");
            n = rheome.select.rollback(db, t1.txn_id);
            tc.verifyEqual(n, 1);
            [L, T] = rheome.select.labels(db);
            tc.verifyEqual(height(L), 1);  tc.verifyEqual(L.txn_id, t2.txn_id);
            tc.verifyEqual(T.txn_id, t2.txn_id);
            tc.verifyError(@() rheome.select.rollback(db, t1.txn_id), 'select:rollback:txn');
        end

        function aTileOutsideTheGridIsRefused(tc)
            fx = selFixture();  db = fx.db;
            tc.verifyError(@() rheome.select.label(db, table(1, 0, db.grid.K0 + 1, 'VariableNames', {'channel_id','level','k'}), Kind="bad"), 'select:label:tile');
        end

    end
end
