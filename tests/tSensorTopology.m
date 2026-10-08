classdef tSensorTopology < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function topologyDemoRunsAndEverySelfCheckPasses(tc)
            close all;
            out = rheome.demos.sensor_topology();
            tc.verifyTrue(out.ok);
            tc.verifyTrue(all(isfinite([out.checks.measured])));
        end

        function topologyDemoOpensFourFigures(tc)
            close all;
            rheome.demos.sensor_topology();
            tc.verifyEqual(numel(findobj('Type','figure')), 4);
        end

        function curlOfAGradientVanishes(tc)
            % Not a tolerance, a theorem: curl(grad(phi)) is identically zero wherever phi is
            % single-valued. If this ever rises above machine noise, the phase gradient is
            % not a gradient and every wavevector in the project is suspect.
            close all;
            out = rheome.demos.sensor_topology();
            tc.verifyLessThan(max([out.stats.curlBulk]), 1e-9);
        end

        function theTwoInvariantsSeparateAllFourFamilies(tc)
            % The claim: winding splits {rotating, spiral} from {travelling, target}, and
            % divergence splits {target, spiral} from {travelling, rotating}. Neither alone
            % isolates the spiral; together they do.
            close all;
            out = rheome.demos.sensor_topology();
            g = @(n) out.stats(strcmp({out.stats.name}, n));
            tc.verifyEqual([g('travelling').winding g('target').winding], [0 0]);
            tc.verifyEqual([g('rotating').winding   g('spiral').winding], [1 1]);
            lo = max(g('travelling').divBulk, g('rotating').divBulk);
            hi = min(g('target').divBulk,     g('spiral').divBulk);
            tc.verifyGreaterThan(hi, 3*lo, 'divergence must separate the two pairs');
        end

        function theSpiralIsTheOnlyFamilyCarryingBoth(tc)
            close all;
            out = rheome.demos.sensor_topology();
            g = @(n) out.stats(strcmp({out.stats.name}, n));
            both = arrayfun(@(s) s.winding > 0 && s.divBulk > 0.8, out.stats);
            tc.verifyEqual(nnz(both), 1);
            tc.verifyEqual(out.stats(both).name, 'spiral');
        end

        function topologyDemoMakesNoSourceClaim(tc)
            src = fileread(which('rheome.demos.sensor_topology'));
            tc.verifyEmpty(regexpi(src, 'leadfield|cortex|cortical|inverse|dipole|anatom', 'once'));
        end

    end
end

% Author: Diellor Basha, 2026
