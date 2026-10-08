classdef tFbGraph < matlab.unittest.TestCase

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function runsAndSelfChecks(tc)
            close all;
            out = rheome.demos.filterbank_graph();
            tc.verifyTrue(out.ok, 'every self-check in the graph demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function opensFiveFigures(tc)
            close all;
            rheome.demos.filterbank_graph();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function recoversEveryPlantedWidthAtTheEnergyCalibration(tc)
            % vertexSpectrum is an ENERGY marginal, so the peak member sits at
            % sqrt(2)*sigma, not sigma. Within one member spacing at 4 voices/octave.
            close all;
            out = rheome.demos.filterbank_graph();
            tc.verifyNumElements(out.recoveredSigma, numel(out.plantedSigma));
            tc.verifyEqual(out.expectedSigma, sqrt(2)*out.plantedSigma, 'RelTol', 1e-12);
            tc.verifyEqual(out.recoveredSigma, out.expectedSigma, 'RelTol', 2^(1/4)-1);
        end

        function needsNoWaveletToolbox(tc)
            src = fileread(which('rheome.demos.filterbank_graph'));
            tc.verifyEmpty(regexp(src, 'cwtfilterbank', 'once'));
        end

        function runsWithoutWarnings(tc)
            % A demo that warns on every run trains people to ignore warnings.
            close all;
            tc.verifyWarningFree(@() rheome.demos.filterbank_graph());
        end

    end
end

% Author: Diellor Basha, 2026
