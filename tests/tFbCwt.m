classdef tFbCwt < matlab.unittest.TestCase

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function runsAndSelfChecks(tc)
            if isempty(ver('wavelet')), tc.assumeFail('Wavelet Toolbox not installed'); end
            close all;
            out = rheome.demos.filterbank_cwt();
            tc.verifyFalse(out.skipped);
            tc.verifyTrue(out.ok, 'every self-check in the CWT demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function opensFiveFigures(tc)
            if isempty(ver('wavelet')), tc.assumeFail('Wavelet Toolbox not installed'); end
            close all;
            rheome.demos.filterbank_cwt();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function qfactorIsScalarBecauseTheCwtIsConstantQ(tc)
            % The defining property of a wavelet transform, and the contrast
            % graphfilterbank's per-member Q is measured against.
            if isempty(ver('wavelet')), tc.assumeFail('Wavelet Toolbox not installed'); end
            close all;
            out = rheome.demos.filterbank_cwt();
            tc.verifySize(out.qfactor, [1 1]);
        end

        function ridgeTracksTheChirpClosely(tc)
            if isempty(ver('wavelet')), tc.assumeFail('Wavelet Toolbox not installed'); end
            close all;
            out = rheome.demos.filterbank_cwt();
            tc.verifyLessThan(out.ridgeErrHz, 5);       % measured 1.86 Hz
        end

        function skipsCleanlyWithoutTheToolbox(tc)
            % Documented behaviour: absence is a skip, not a failure. Verified by reading
            % the guard rather than uninstalling a toolbox.
            src = fileread(which('rheome.demos.filterbank_cwt'));
            tc.verifySubstring(src, "ver('wavelet')");
            tc.verifySubstring(src, 'skipped');
        end

    end
end

% Author: Diellor Basha, 2026
