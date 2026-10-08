classdef tSelectMeasure < matlab.unittest.TestCase
% The SECOND bookkeeping matrix: measurements hung on pyramid nodes, and the join back.
%
% The moments are dense and mergeable; these are sparse, heterogeneous and merge not at all.
% What is asserted here is the contract that keeps them apart: a measurement stays at the
% node it was made on, nothing rolls it up, and the only thing that crosses levels in a
% context join is the merged statistic, never the measurement.
%
% Author: Diellor Basha, 2026

    properties
        fx; db; rows
    end

    methods (TestClassSetup)
        function open(tc)
            tc.fx = selFixture();
            tc.db = tc.fx.db;
        end
    end

    methods (TestMethodSetup)
        function cleanSidecar(tc)
            if exist(tc.db.labelFile, 'file') == 2, delete(tc.db.labelFile); end
            b = tc.db.bands.j(tc.db.bands.fLo == 8);
            tc.rows = table(repmat(2, 5, 1), [3 4 9 10 16]', [2 3 1 4 2]', ones(5, 1), repmat(b, 5, 1), ...
                            'VariableNames', {'level','k','value','unit_id','band_id'});
        end
    end

    methods (TestMethodTeardown)
        function wipe(tc)
            if exist(tc.db.labelFile, 'file') == 2, delete(tc.db.labelFile); end
        end
    end

    methods (Test)

        function theSchemaSaysWhichMatrixEachRelationIs(tc)
            % ⭐ The split is machine-readable, so "this is a summary of its children" and
            % "this was measured here" cannot be confused by a reader or by a query planner.
            S = rheome.select.schema();
            k = {S.kind};
            tc.verifyTrue(all(ismember(k, {'dimension','mergeable','note'})));
            nm = string({S.name});
            tc.verifyEqual(string(k(nm == "feature")), "mergeable");
            tc.verifyEqual(string(k(nm == "feature_group_band")), "mergeable");
            tc.verifyEqual(string(k(nm == "measure")), "note");
            tc.verifyEqual(string(k(nm == "label")), "note");
            tc.verifyEqual(string(k(nm == "tile")), "dimension");
            m = S(nm == "measure");
            tc.verifyEqual(m.columns, {'measure_id','recording_id','scope','unit_id','level','k','band_id','kind','value','source','txn_id'});
            d = rheome.select.ddl();
            tc.verifySubstring(d, 'CREATE TABLE measure');
            tc.verifySubstring(d, 'CREATE INDEX measure_node');
        end

        function aWriteIsATransactionAndIsIdempotent(tc)
            t = rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel", Source="flow");
            tc.verifyFalse(t.existed);
            tc.verifyEqual(t.count, 5);
            t2 = rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel", Source="flow");
            tc.verifyTrue(t2.existed);
            tc.verifyEqual(t2.txn_id, t.txn_id);
            tc.verifyEqual(height(rheome.select.measures(tc.db)), 5);        % not written twice
            t3 = rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel", Source="manual");
            tc.verifyFalse(t3.existed);                               % a different source is a different write
            M = rheome.select.measures(tc.db);
            tc.verifyEqual(height(M), 10);
            tc.verifyEqual(sort(unique(M.measure_id)), (1:10)');
        end

        function measurementsAreAddressedByScopeAndUnit(tc)
            % `label` knows only channels; a flow measurement usually belongs to a patch.
            g = tc.rows;  g.unit_id = repmat(tc.db.groupNodes(1), 5, 1);
            rheome.select.measure(tc.db, g, Kind="rotation", Scope="group");
            rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            M = rheome.select.measures(tc.db, Scope="group");
            tc.verifyEqual(height(M), 5);
            tc.verifyTrue(all(M.unit_id == tc.db.groupNodes(1)));
            tc.verifyEqual(unique(M.kind), "rotation");
            tc.verifyEqual(height(rheome.select.measures(tc.db, Scope="channel")), 5);
        end

        function aBadNodeOrUnitIsRefusedAtTheWrite(tc)
            bad = tc.rows;  bad.k(1) = 1e6;
            tc.verifyError(@() rheome.select.measure(tc.db, bad, Kind="x", Scope="channel"), 'select:measure:tile');
            bad = tc.rows;  bad.level(1) = 99;
            tc.verifyError(@() rheome.select.measure(tc.db, bad, Kind="x", Scope="channel"), 'select:measure:tile');
            bad = tc.rows;  bad.unit_id(1) = 999;
            tc.verifyError(@() rheome.select.measure(tc.db, bad, Kind="x", Scope="channel"), 'select:measure:unit');
            bad = tc.rows;  bad.unit_id(1) = 999;
            tc.verifyError(@() rheome.select.measure(tc.db, bad, Kind="x", Scope="group"), 'select:measure:unit');
            tc.verifyError(@() rheome.select.measure(tc.db, tc.rows(:, {'level','k'}), Kind="x"), 'select:measure:rows');
        end

        function readingFiltersOnEveryKey(tc)
            rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            tc.verifyEqual(height(rheome.select.measures(tc.db, Kind="nothing")), 0);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Level=2)), 5);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Level=3)), 0);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Min=3)), 2);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Max=1)), 1);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Units=1)), 5);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Units=2)), 0);
            ext = tc.db.grid.tExtent(3);
            w = [(3-1)*ext, 4*ext];                                   % tiles 3 and 4
            tc.verifyEqual(height(rheome.select.measures(tc.db, Window=w)), 2);
        end

        function nothingRollsUpAndTheContractSaysSo(tc)
            % ⚠ THE WHOLE POINT OF THE SECOND MATRIX. Five measurements at level 2 leave the
            % level-3 nodes above them empty: a value here is not a summary of anything.
            rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            tc.verifyEqual(height(rheome.select.measures(tc.db, Level=3)), 0);
            tc.verifyEqual(height(rheome.select.measures(tc.db, Level=1)), 0);
            S = rheome.select.schema();
            tc.verifyEqual(string({S(string({S.name}) == "measure").kind}), "note");
        end

        function oneTransactionLogCoversLabelsAndMeasurements(tc)
            % A rollback must not leave half a write behind, and a label write must not drop
            % measurements written earlier -- they share one sidecar file.
            tm = rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            tl = rheome.select.label(tc.db, tc.rows(:, {'level','k'}), Kind="bad");
            tc.verifyNotEqual(tl.txn_id, tm.txn_id);
            tc.verifyEqual(height(rheome.select.measures(tc.db)), 5);        % survived the label write
            tc.verifyEqual(height(rheome.select.labels(tc.db)), 5);
            n = rheome.select.rollback(tc.db, tm.txn_id);
            tc.verifyEqual(n, 5);
            tc.verifyEqual(height(rheome.select.measures(tc.db)), 0);
            tc.verifyEqual(height(rheome.select.labels(tc.db)), 5);          % labels untouched
            n = rheome.select.rollback(tc.db, tl.txn_id);
            tc.verifyEqual(n, 5);
            tc.verifyEqual(height(rheome.select.labels(tc.db)), 0);
            tc.verifyError(@() rheome.select.rollback(tc.db, 99), 'select:rollback:txn');
        end

        function aSidecarWrittenBeforeMeasurementsExistedStillReads(tc)
            % Old stores carry only label and txn; loading must not error.
            L = table(1, string(tc.db.recording_id), 0, 2, 3, 0, "bad", 1, "manual", 1, ...
                      'VariableNames', {'label_id','recording_id','channel_id','level','k','band_id','kind','value','source','txn_id'});
            T = table(1, string(tc.db.recording_id), "bad", 1, "deadbeef", "", "now", ...
                      'VariableNames', {'txn_id','recording_id','kind','count','content_hash','author','created'});
            label = L;  txn = T;                                      %#ok<NASGU>
            save(tc.db.labelFile, 'label', 'txn', '-v7');
            tc.verifyEqual(height(rheome.select.measures(tc.db)), 0);
            tc.verifyEqual(height(rheome.select.labels(tc.db)), 1);
            t = rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            tc.verifyEqual(t.txn_id, 2);                              % the log continues
            tc.verifyEqual(height(rheome.select.labels(tc.db)), 1);
        end

        function theContextJoinAttachesTheStateEachMeasurementOccurredIn(tc)
            rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            M = rheome.select.measures(tc.db);
            C = rheome.select.context(tc.db, M, Stats=["rms","spectralCentroid"], Ancestors=0);
            tc.verifyEqual(height(C), height(M));
            R = rheome.select.derive(tc.db, Level=2, Channels=1, Stats=["rms","spectralCentroid"]);
            for i = 1:height(C)
                j = find(R.k == C.k(i), 1);
                tc.verifyEqual(C.rms(i), R.rms(j), 'RelTol', 1e-12);
                tc.verifyEqual(C.spectralCentroid(i), R.spectralCentroid(j), 'RelTol', 1e-12);
            end
        end

        function anAncestorsContextComesFromTheNodeAboveByKeyArithmetic(tc)
            rheome.select.measure(tc.db, tc.rows, Kind="vortex_count", Scope="channel");
            M = rheome.select.measures(tc.db);
            C = rheome.select.context(tc.db, M, Stats="rms", Ancestors=[0 2]);
            tc.verifyEqual(C.k_up2, floor((C.k - 1) / 4) + 1);
            tc.verifyEqual(C.level_up2, C.level + 2);
            R4 = rheome.select.derive(tc.db, Level=4, Channels=1, Stats="rms");
            for i = 1:height(C)
                j = find(R4.k == C.k_up2(i), 1);
                tc.verifyEqual(C.rms_up2(i), R4.rms(j), 'RelTol', 1e-12);
            end
            tc.verifyNotEqual(C.rms, C.rms_up2);                      % two scales, two numbers
            tc.verifyError(@() rheome.select.context(tc.db, M, Ancestors=-1), 'select:context:ancestor');
        end

        function aPerBandContextIsTakenInTheRowsOwnBand(tc)
            % A vortex found in one band is scored against that band, not against band 1.
            b = tc.db.bands.j(tc.db.bands.fLo == 8);
            L = tc.db.bands.naturalLevel(tc.db.bands.j == b);
            r = tc.rows;  r.level = repmat(L, 5, 1);  r.k = (1:5)';  r.band_id = repmat(b, 5, 1);
            rheome.select.measure(tc.db, r, Kind="vortex_count", Scope="channel");
            M = rheome.select.measures(tc.db);
            C = rheome.select.context(tc.db, M, Stats="share", Ancestors=0);
            R = rheome.select.derive(tc.db, Level=L, Channels=1, Stats="share");
            bAt = R.Properties.CustomProperties.bands;
            col = find(bAt == b, 1);
            for i = 1:height(C)
                j = find(R.k == C.k(i), 1);
                tc.verifyEqual(C.share(i), R.share(j, col), 'RelTol', 1e-12);
            end
        end

        function aContextBelowTheDiagonalIsNaNRatherThanAnError(tc)
            % The band is not carried at that level, so there is nothing to join to. The
            % join must not fail, and must not borrow a neighbouring band's number.
            b = tc.db.bands.j(tc.db.bands.naturalLevel == max(tc.db.bands.naturalLevel));
            b = b(1);
            r = tc.rows;  r.band_id = repmat(b, 5, 1);  r.level = zeros(5, 1);  r.k = (1:5)';
            rheome.select.measure(tc.db, r, Kind="vortex_count", Scope="channel");
            M = rheome.select.measures(tc.db);
            C = rheome.select.context(tc.db, M, Stats=["share","rms"], Ancestors=0);
            tc.verifyTrue(all(isnan(C.share)));
            tc.verifyTrue(all(~isnan(C.rms)));                        % a time stat needs no band
        end

    end
end

% Author: Diellor Basha, 2026
