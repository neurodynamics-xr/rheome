classdef tRatePlan < matlab.unittest.TestCase
% rheome.flow.rateplan -- the per-filter output rate of a constant-Q bank.
% Author: Diellor Basha, 2026

    methods (Test)
        function ratesAreSamplesPerCycleTimesCentreFrequency(tc)
            R = rheome.flow.rateplan([10 20], 2400, 30);
            tc.verifyEqual(R.rate, [300 600], 'AbsTol', 1e-9);
            tc.verifyEqual(R.capped, [false false]);
        end

        function theAcquisitionRateIsAHardCeiling(tc)
            % 30 samples/cycle needs 30*f <= fs. Above fs/30 the recording simply cannot
            % supply it, and no processing choice recovers it -- so it must be REPORTED,
            % not silently rounded down.
            R = rheome.flow.rateplan([10 60], 600, 30);
            tc.verifyEqual(R.rate, [300 600], 'AbsTol', 1e-9);
            tc.verifyEqual(R.capped, [false true]);
            tc.verifyEqual(R.ceilingHz, 20, 'AbsTol', 1e-9);
            tc.verifyEqual(R.actualPerCycle(2), 10, 'AbsTol', 1e-9);   % 600/60, not 30
        end

        function decimationFactorsNeverGoBelowOne(tc)
            R = rheome.flow.rateplan([100 200], 600, 30);      % both want more than fs
            tc.verifyTrue(all(R.rate <= 600));
            tc.verifyTrue(all(R.decim >= 1));
        end

        function reportsTheBankTotalAgainstFlatAndCritical(tc)
            R = rheome.flow.rateplan([10 20 40], 2400, 30, [2 4 8]);
            tc.verifyEqual(R.flatTotal, 3*2400, 'AbsTol', 1e-9);
            tc.verifyEqual(R.total, 300+600+1200, 'AbsTol', 1e-9);
            tc.verifyEqual(R.criticalTotal, 14, 'AbsTol', 1e-9);
        end

        function rejectsANonPositiveTarget(tc)
            tc.verifyError(@() rheome.flow.rateplan([10 20], 600, 0), 'flow:rateplan:target');
        end
    end
end

% Author: Diellor Basha, 2026
