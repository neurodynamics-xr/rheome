classdef tFbJoint < matlab.unittest.TestCase

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function runsAndSelfChecks(tc)
            close all;
            out = rheome.demos.filterbank_joint();
            tc.verifyTrue(out.ok, 'every self-check in the joint demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function opensFiveFigures(tc)
            close all;
            rheome.demos.filterbank_joint();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function recoversThePlantedSpeed(tc)
            close all;
            out = rheome.demos.filterbank_joint();
            tc.verifyEqual(out.recoveredSlope, out.plantedSlope, 'RelTol', 0.10);
            tc.verifyGreaterThan(out.slopeR2, 0.95);
        end

        function aClippedBandCannotMeasureTheSlope(tc)
            % The documented failure mode, demonstrated rather than asserted: a retention
            % that does not span the diagonal has a narrower support and a worse fit.
            close all;
            out = rheome.demos.filterbank_joint();
            tc.verifyLessThan(diff(out.clippedSupport), diff(out.fSupport));
        end

        function needsNoWaveletToolbox(tc)
            src = fileread(which('rheome.demos.filterbank_joint'));
            tc.verifyEmpty(regexp(src, 'cwtfilterbank', 'once'));
        end

        function runsWithoutWarnings(tc)
            close all;
            tc.verifyWarningFree(@() rheome.demos.filterbank_joint());
        end

    end
end

% Author: Diellor Basha, 2026
