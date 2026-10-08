classdef tSensorWindow < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (Test)

        function theWindowFloorIsTwiceThePitch(tc)
            w = rheome.sensors.window(rheome.sensors.grid('utah'));
            tc.verifyEqual(w.LambdaMin, 2*400e-6, 'RelTol', 1e-12);
        end

        function theWindowCeilingIsTheApertureWhichIsTheDiagonal(tc)
            w = rheome.sensors.window(rheome.sensors.grid('utah'));
            tc.verifyEqual(w.LambdaMax, 9*400e-6*sqrt(2), 'RelTol', 1e-12);
        end

        function utahAndEcogWindowsAreDisjoint(tc)
            % Same geometry class, pitches two orders apart. That their windows do not
            % overlap at all is the suite's cleanest statement that a window belongs to the
            % instrument and not to the method.
            wu = rheome.sensors.window(rheome.sensors.grid('utah'));
            we = rheome.sensors.window(rheome.sensors.grid('ecog'));
            tc.verifyLessThan(wu.LambdaMax, we.LambdaMin);
        end

        function aProbeSupportsNoWindingBecauseItHasNoFaces(tc)
            w = rheome.sensors.window(rheome.sensors.linear('NumSites', 32));
            tc.verifyFalse(w.HasFaces);
            tc.verifyFalse(w.Supports.winding);
            tc.verifyFalse(w.Supports.chirality);
            tc.verifyFalse(w.Supports.direction);
        end

        function aProbeStillSupportsDispersionAndSpeed(tc)
            % The probe's limits are of two KINDS. Direction and winding are unavailable;
            % omega(lambda) is not -- it is available in projection along the shank.
            w = rheome.sensors.window(rheome.sensors.linear('NumSites', 32));
            tc.verifyTrue(w.Supports.dispersion);
            tc.verifyTrue(w.Supports.speed);
        end

        function aGridSupportsEverything(tc)
            w = rheome.sensors.window(rheome.sensors.grid('ecog'));
            tc.verifyTrue(w.HasFaces);
            tc.verifyTrue(all(cell2mat(struct2cell(w.Supports))));
        end

        function octavesCountsTheUsableRange(tc)
            w = rheome.sensors.window(rheome.sensors.grid('utah'));
            tc.verifyEqual(w.Octaves, log2(w.LambdaMax/w.LambdaMin), 'RelTol', 1e-12);
            tc.verifyGreaterThan(w.Octaves, 2);
        end

    end
end

% Author: Diellor Basha, 2026
