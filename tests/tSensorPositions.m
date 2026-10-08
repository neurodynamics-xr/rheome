classdef tSensorPositions < matlab.unittest.TestCase
% An array from bare coordinates. The general case the four named constructors specialise.

    methods (Test)

        function itBuildsTheSameContractTheOtherConstructorsDo(tc)
            ref = rheome.sensors.grid('ecog');
            arr = rheome.sensors.positions('copy', ref.Pos);
            tc.verifyEqual(arr.Pos, ref.Pos);
            tc.verifyEqual(arr.nCh, ref.nCh);
            tc.verifyEqual(arr.Dim, ref.Dim);
            % Pitch and aperture are measured from the coordinates, so they must agree
            % with the preset that produced them.
            tc.verifyEqual(arr.Pitch,    ref.Pitch,    'RelTol', 1e-12);
            tc.verifyEqual(arr.Aperture, ref.Aperture, 'RelTol', 1e-12);
        end

        function twoColumnCoordinatesAreAcceptedAsPlanar(tc)
            % The report carries each array as a PROJECTED plane, two columns. Requiring
            % three would mean every caller pads by hand and one of them gets it wrong.
            ref = rheome.sensors.grid('ecog');
            arr = rheome.sensors.positions('flat', ref.Pos(:,1:2));
            tc.verifySize(arr.Pos, [ref.nCh, 3]);
            tc.verifyEqual(arr.Pos(:,3), zeros(ref.nCh,1));
            tc.verifyEqual(arr.Dim, 2);
        end

        function dimensionIsInferredFromTheCoordinates(tc)
            % A line of contacts is Dim 1 however it is oriented; a lattice is Dim 2.
            n = 32;
            line = [linspace(0,1,n).', zeros(n,1), zeros(n,1)] * 1e-3;
            tc.verifyEqual(rheome.sensors.positions('line', line).Dim, 1);
            tc.verifyEqual(rheome.sensors.positions('lattice', rheome.sensors.grid('utah').Pos).Dim, 2);
        end

        function anObliqueLineIsStillOneDimensional(tc)
            % ⚠ RANK, NOT AXIS ALIGNMENT. A probe inserted at an angle has spread on all
            % three coordinates and is still a chain.
            n = 40;  t = linspace(0, 1, n).';
            obliq = [t, 2*t, -0.5*t] * 1e-3;
            tc.verifyEqual(rheome.sensors.positions('oblique', obliq).Dim, 1);
        end

        function theArrayItReturnsDrivesTheOperatorUnchanged(tc)
            % The point of matching the contract: everything downstream accepts it.
            arr = rheome.sensors.positions('copy', rheome.sensors.grid('ecog').Pos);
            G = rheome.sensors.graph(arr);
            tc.verifyEqual(G.nV, arr.nCh);
            M = rheome.sensors.modes(G);
            tc.verifyEqual(numel(M.Lambda), arr.nCh);
            tc.verifyLessThan(abs(M.Lambda(1)), 1e-8);
        end

        function labelsDefaultButCanBeGiven(tc)
            P = rheome.sensors.grid('ecog').Pos;
            tc.verifyEqual(numel(rheome.sensors.positions('a', P).Labels), size(P,1));
            lab = compose('S%d', 1:size(P,1));
            tc.verifyEqual(rheome.sensors.positions('a', P, 'Labels', lab).Labels{7}, 'S7');
        end

        function badCoordinatesAreRefused(tc)
            tc.verifyError(@() rheome.sensors.positions('a', zeros(4,5)), 'sensors:positions:pos');
            tc.verifyError(@() rheome.sensors.positions('a', zeros(1,3)), 'sensors:array:count');
        end

    end
end

% Author: Diellor Basha, 2026
