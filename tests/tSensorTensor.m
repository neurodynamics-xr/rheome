classdef tSensorTensor < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function tensorDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_tensor();
            tc.verifyTrue(out.ok);
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function tensorDemoOpensFiveFigures(tc)
            close all;
            rheome.demos.sensor_tensor();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function theBankIsAnApertureLadderInMillimetres(tc)
            % Without the calibration a filterbank on a graph is a list of dimensionless
            % numbers. The claim is that these are LENGTHS, ordered, and on the array's scale.
            close all;
            out = rheome.demos.sensor_tensor();
            tc.verifyTrue(all(diff(out.bank) < 0), 'apertures must descend with member index');
            tc.verifyGreaterThan(min(out.bank), 5);      % mm
            tc.verifyLessThan(max(out.bank), 400);
        end

        function sizeAndRateAreSeparatedIntoDifferentCells(tc)
            % Two components differing in BOTH scale and rate must land in distinct tensor
            % cells, and be recovered in physical units.
            close all;
            out = rheome.demos.sensor_tensor();
            tc.verifyNumElements(out.separated, 2);
            tc.verifyGreaterThan(out.separated(1).lam, 2*out.separated(2).lam);
            tc.verifyLessThan(out.separated(1).f, out.separated(2).f);
        end

        function aSpeedMemberDoesWhatNoRectangleCan(tc)
            % speedSel = [slow|diagonal, fast|diagonal, slow|rectangle, fast|rectangle].
            % The rectangle admits both trains; the diagonal admits one and rejects the
            % other. That contrast is the reason a non-separable member exists at all.
            close all;
            out = rheome.demos.sensor_tensor();
            s = out.speedSel;
            tc.verifyGreaterThan(s(3), 0.6, 'rectangle must admit the slow train');
            tc.verifyGreaterThan(s(4), 0.6, 'rectangle must admit the fast train');
            tc.verifyGreaterThan(s(1), 0.4, 'diagonal must admit its own train');
            tc.verifyLessThan(s(2), 0.05,   'diagonal must reject the other train');
            tc.verifyGreaterThan(s(1)/max(s(2),1e-6), 20);
        end

        function tensorDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_tensor'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
