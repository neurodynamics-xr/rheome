classdef tSensorPhase < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function phaseDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_phase();
            tc.verifyTrue(out.ok, 'every self-check in the phase demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function phaseDemoOpensFiveFigures(tc)
            close all;
            rheome.demos.sensor_phase();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function threeIndependentRoutesAgreeOnTheWavenumber(tc)
            % The montage phase slope, the per-face phase gradient, and the pairwise lag are
            % three different computations over the same field. If they ever disagree, one of
            % them has broken -- and which one is not obvious from any single number.
            close all;
            out = rheome.demos.sensor_phase();
            k = [out.kFromPhase, out.kFromGradient, out.lagSlope];
            tc.verifyLessThan((max(k)-min(k))/mean(k), 0.10, ...
                'phase slope, phase gradient and pairwise lag must agree on k');
        end

        function theRotorHasExactlyOneCoreOfUnitCharge(tc)
            close all;
            out = rheome.demos.sensor_phase();
            tc.verifyNumElements(out.rotor.charge, 1);
            tc.verifyEqual(abs(out.rotor.charge), 1);
        end

        function theAnalyticSignalHasTheRightEnvelopeAndPhase(tc)
            % A pure tone must come back with a flat envelope and a phase advancing at
            % exactly omega. This is the FFT route standing in for hilbert(), which is
            % Signal Processing Toolbox and unavailable here.
            fs = 200; t = (0:511)/fs; f = 12;
            z = rheome.demos.sensor_dynamics_selftest('analytic', cos(2*pi*f*t));
            mid = 100:400;                       % away from the transform's edges
            tc.verifyEqual(abs(z(mid)), ones(1,numel(mid)), 'RelTol', 0.02);
            dph = diff(unwrap(angle(z(mid))));
            tc.verifyEqual(mean(dph), 2*pi*f/fs, 'RelTol', 0.01);
        end

        function theAnalyticSignalIsInvariantToAmplitudeInItsPhase(tc)
            fs = 200; t = (0:255)/fs;
            u = cos(2*pi*9*t);
            z1 = rheome.demos.sensor_dynamics_selftest('analytic', u);
            z2 = rheome.demos.sensor_dynamics_selftest('analytic', 37*u);
            mid = 60:200;
            tc.verifyEqual(angle(z1(mid)), angle(z2(mid)), 'AbsTol', 1e-9);
        end

        function phaseDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_phase'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
