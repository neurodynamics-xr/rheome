classdef tGeomLadder < matlab.unittest.TestCase
% rheome.geom.ladder: the spatial ladder's structure, which is derivable and must therefore be asserted.
%
% Author: Diellor Basha, 2026
    properties
        L; U
    end
    methods (TestClassSetup)
        function build(t)
            try
                t.L = rheome.geom.ladder(rheomeTestSubject(), MaxDepth=8);
            catch
                t.assumeFail('cached bases for test subject are not present');
            end
            t.U = t.L.Properties.UserData;
        end
    end
    methods (Test)
        function theTreeIsDyadicInAreaSoDiameterGoesBySqrtTwo(t)
            % ⚠ the factor that makes one spatial octave TWO tree depths
            r = t.U.diameterMM(1:end-1) ./ t.U.diameterMM(2:end);
            t.verifyEqual(mean(r), sqrt(2), 'RelTol', 0.05, ...
                'rheome.geom.tree bisects by area, so diameter must fall by sqrt(2) per depth');
        end
        function nodeCountDoublesPerDepth(t)
            t.verifyEqual(t.U.nNodes(:)', 2.^(1:numel(t.U.nNodes)));
        end
        function weylPredictsTheModeCount(t)
            % counted from the cached spectrum vs A k^2 / 4pi. The finest octave is truncated by the
            % basis (1000 modes), so only the resolved ones are checked.
            ok = t.L.nModes > 5 & t.L.nModes < 300;
            r = t.L.nModes(ok) ./ t.L.nModesWeyl(ok);
            t.verifyEqual(r, ones(size(r)), 'AbsTol', 0.12, 'Weyl should hold to ~10%');
        end
        function theDiagonalStepsTwoDepthsPerOctave(t)
            d = diff(t.L.naturalDepth);
            d = d(d > 0);                       % the first edge saturates at the coarsest tile
            t.verifyTrue(all(d == 2), 'a factor of 2 in wavelength must cost 2 tree depths');
        end
        function modesPerCellIsConstantAlongTheDiagonal(t)
            % ⭐ the twin of constant cycles-per-tile: 2^-D k^2 is constant when D steps 2 per octave
            m = t.L.modesPerCell(t.L.naturalDepth > min(t.L.naturalDepth));
            t.verifyLessThan(std(m)/mean(m), 0.12, ...
                'modes per (tile, octave) cell must not drift along the diagonal');
        end
        function theCoarsestOctaveIsForbidden(t)
            % ⚠ correct, not a bug: it needs a tile larger than the hemisphere
            t.verifyTrue(t.L.forbidden(1));
        end
        function theTwoAdmissibilityTestsAreDifferent(t)
            % ⚠ node size and wavelength are separate questions and must not be conflated
            t.verifyEqual(t.L.nodeAdmissible, t.L.nodeDiameterMM > 2*t.U.R50mm);
            t.verifyEqual(t.L.waveAdmissible, t.L.wavelengthCenterMM > 2*t.U.R50mm);
            t.verifyFalse(isequal(t.L.nodeAdmissible, t.L.waveAdmissible));
        end
    end
end
