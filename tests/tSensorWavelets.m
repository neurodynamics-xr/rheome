classdef tSensorWavelets < matlab.unittest.TestCase
% Scale, rate and speed as three controls, for GENERATION and for ANALYSIS.

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function demoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_wavelets();
            tc.verifyTrue(out.ok, 'every self-check in the wavelet-tensor demo must pass');
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function demoOpensFiveFigures(tc)
            close all;
            rheome.demos.sensor_wavelets();
            tc.verifyEqual(numel(findobj('Type','figure')), 5);
        end

        function everyPlantedScaleComesBack(tc)
            % The graph member IS the spatial scale control. Planting one and recovering a
            % different wavelength would mean the ladder is decorative.
            close all;
            out = rheome.demos.sensor_wavelets();
            tc.verifyGreaterThanOrEqual(numel(out.scale.err), 4);
            tc.verifyLessThan(max(abs(out.scale.err)), 0.08);
        end

        function everyPlantedRateComesBack(tc)
            close all;
            out = rheome.demos.sensor_wavelets();
            tc.verifyGreaterThanOrEqual(numel(out.rate.err), 4);
            tc.verifyLessThan(max(abs(out.rate.err)), 0.03);
        end

        function everyPlantedSpeedComesBack(tc)
            % ⚠ SPEED IS THE ONE THAT NEEDS THE FULL WAVENUMBER AXIS. dispersion says the
            % members are not involved and the band must span the diagonal; a demo that
            % fits speed inside one member or one narrow band recovers a confident wrong
            % number, R2 near 1 and all.
            close all;
            out = rheome.demos.sensor_wavelets();
            tc.verifyGreaterThanOrEqual(numel(out.speed.err), 3);
            tc.verifyLessThan(max(abs(out.speed.err)), 0.10);
        end

        function theControlsAreIndependent(tc)
            % Changing the spatial scale must not move the recovered speed.
            close all;
            out = rheome.demos.sensor_wavelets();
            tc.verifyLessThan(out.independence, 0.10);
        end

        function theDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_wavelets'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
