classdef tScalePlantFloors < matlab.unittest.TestCase
% The MS1 group plant measures (nsp cf-plant-floors: G1, G2, G9, G11) end to end on the synthetic cached
% subject: table shapes, balanced cells and halves, the no-instrument arm beating chance, and the driver
% writing each analysis's per-plant table. Small settings: the numbers are checked for sense, not value.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_plantfloors'
        S
    end

    methods (TestClassSetup)
        function cache(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, tc.Name));
            tc.S = rheome.scale.sensors(tc.Name);
        end
    end

    methods (Test)

        function floorsHaveBalancedCellsAndHalves(tc)
            [T, X] = rheome.scale.measure_plantfloors(tc.Name, tc.S, LambdaMM=[30 60], PlantsPerCell=2, ...
                                                      SNRdB=[Inf 0], Hemis="L");
            tc.verifyEqual(height(X), 2*2 * 2 * 3);                 % fields x types x (2 SNR + direct)
            for ty = ["source" "vortex"]
                for l = [30 60]
                    k = X.type == ty & X.lambdaMM == l & X.condition == "direct";
                    tc.verifyEqual(sum(k), 2);  tc.verifyEqual(sort(X.half(k))', [1 2]);
                end
            end
            d = X(X.condition == "direct", :);
            tc.verifyLessThan(median(d.errMM), median(d.chanceMM));  % the planted field locates itself
            tc.verifyTrue(all(ismember(["err_mm" "floor_mm" "chance_mm" "direct_err_mm" "normal_frac"], T.metric)));
            tc.verifyEqual(sum(T.metric == "floor_mm"), 2 * 2 * 3);  % type x SNR x (all, h1, h2)
        end

        function movingVortexReadsEveryFrame(tc)
            [T, X] = rheome.scale.measure_movingvortex(tc.Name, tc.S, Placements=2, SigmaMM=[20 40], SNRdB=[Inf 0], Hemis="L");
            tc.verifyEqual(height(X), 2*2*2);
            tc.verifyGreaterThan(min(X.nFrames), 2);
            tc.verifyTrue(all(isfinite([X.coreErrMM; X.directErrMM; X.chanceMM; X.netRatio])));
            tc.verifyEqual(sort(unique(X.half))', [1 2]);
            tc.verifyTrue(any(T.metric == "core_err_mm" & T.band == "s20_0dB_h2"));
        end

        function noiseFloorItemsReturnTheirMetrics(tc)
            want = struct('composition', "sol_gap_vortex_emptyroom", 'sizeruler', "slope", ...
                          'vortexscale', "wavelength_mm", 'rotation', "icoh_rotating", 'diracangles', "mean_cos2");
            for it = string(fieldnames(want))'
                [T, X] = rheome.scale.measure_noisefloor(tc.Name, tc.S, it, Hemis="L", Seeds=2, NumFrames=3, ...
                                                         SNRdB=[Inf 10], MomentNAm=[Inf 10]);
                tc.verifyNotEmpty(X, it);
                tc.verifyTrue(any(T.metric == want.(it)), it);
                tc.verifyTrue(all(T.analysis == it), it);
            end
            [T, X] = rheome.scale.measure_noisefloor(tc.Name, tc.S, "composition", Hemis="L", Seeds=2, NumFrames=3);
            v = X(X.source == "vortex_planted" & X.region == "whole", :);
            tc.verifyGreaterThan(median(v.solFrac), median(v.irrFrac));   % a seeded vortex is solenoidal before the inverse
            tc.verifyLessThan(abs(median(v.normalFrac)), 1e-3);           % and tangential (gauge vs vertex normals)
            tc.verifyTrue(any(T.band == "leadfield_rows_whole"));
            [~, R] = rheome.scale.measure_noisefloor(tc.Name, tc.S, "rotation", Hemis="L", Seeds=2, MomentNAm=Inf);
            % noiseless: a rotating source has quadrature between the two patterns, a standing one none.
            % (The synthetic gain is not magnetic and nearly blind to a vortex, so the size, 0.97 on
            % one participant, is not checked here; the canary on real heads is.)
            tc.verifyGreaterThan(median(R.iCoh(R.kind == "rotating")), 0.05);
            tc.verifyLessThan(median(abs(R.iCoh(R.kind == "standing"))), 1e-6);
        end

        function theDriverWritesThePlantTables(tc)
            od = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            R = rheome.scale.run(tc.Name, Analyses=["diracangles" "rotation"], OutDir=od);
            tc.verifyEqual(R.timing.status', ["ok" "ok"], strjoin(R.timing.message, ' | '));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'diracangles.csv')));
            tc.verifyTrue(isfile(fullfile(od, tc.Name, 'rotation.csv')));
            tc.verifyTrue(any(R.metrics.analysis == "rotation"));
        end
    end
end

% Author: Diellor Basha, 2026
