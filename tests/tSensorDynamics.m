classdef tSensorDynamics < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function dynamicsDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_dynamics();
            tc.verifyTrue(out.ok, 'every self-check in the dynamics demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function dynamicsDemoOpensFiveFigures(tc)
            close all;
            rheome.demos.sensor_dynamics();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function standingAndTravellingSeparateAtOppositeEndsOfTheIndex(tc)
            % The phase's central claim. A standing wave collapses every sensor onto one line
            % through the origin of the complex plane, so the smaller eigenvalue vanishes; a
            % travelling wave spreads them around a circle, so the two are comparable. If
            % these ever converge, the discriminator has stopped working.
            close all;
            out = rheome.demos.sensor_dynamics();
            tc.verifyLessThan(out.standingIndex(2), 0.02);          % standing
            tc.verifyGreaterThan(out.standingIndex(1), 0.35);       % travelling
            tc.verifyGreaterThan(out.standingIndex(1) - out.standingIndex(2), 0.3);
        end

        function theIndexIsInvariantToAmplitude(tc)
            % It must measure PHASE STRUCTURE, not size. Scaling the field by 100 cannot
            % change whether it travels.
            z = exp(1i*linspace(0, 6*pi, 60)).';
            s1 = rheome.demos.sensor_dynamics_selftest('travelindex', z);
            s2 = rheome.demos.sensor_dynamics_selftest('travelindex', 100*z);
            tc.verifyEqual(s1, s2, 'RelTol', 1e-10);
            tc.verifyGreaterThan(s1, 0.4);
        end

        function aRealFieldReadsAsStanding(tc)
            % A purely real spatial pattern has one common phase by construction.
            z = randn(80,1);
            tc.verifyLessThan(rheome.demos.sensor_dynamics_selftest('travelindex', z), 1e-10);
        end

        function plantedWaveSpeedComesBack(tc)
            close all;
            out = rheome.demos.sensor_dynamics();
            tc.verifyEqual(out.speedFit.slope, 0.35, 'RelTol', 0.15);
            tc.verifyGreaterThan(out.speedFit.r2, 0.9);
        end

        function dynamicsDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_dynamics'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
