classdef tJfbBand < matlab.unittest.TestCase
% The band member is a FILTER, not a truncation. Its edges matter.

    methods (Test)

        function reachesExactlyZeroAtBothEdges(tc)
            % This is what makes truncating to the retained bins lossless rather than
            % approximately so: the discarded bins lie beyond a gain already zero.
            h = rheome.jointfilterbank.band(8, 13);
            tc.verifyEqual(h(2*pi*8),  0, 'AbsTol', 0);
            tc.verifyEqual(h(2*pi*13), 0, 'AbsTol', 0);
        end

        function isUnityInTheInterior(tc)
            h = rheome.jointfilterbank.band(8, 13);
            tc.verifyEqual(h(2*pi*10.5), 1, 'RelTol', 1e-12);
        end

        function isZeroOutside(tc)
            h = rheome.jointfilterbank.band(8, 13);
            tc.verifyEqual(h(2*pi*[1 5 20 40]), zeros(1,4), 'AbsTol', 0);
        end

        function widthZeroIsABrickWall(tc)
            % Required for any pass whose spectrum will be FITTED: a taper removes the
            % lowest bins, which carry most of the leverage on a 1/f slope.
            h = rheome.jointfilterbank.band(8, 13, 0);
            tc.verifyEqual(h(2*pi*8.0001), 1, 'RelTol', 1e-9);
        end

        function defaultWidthIsCappedOnAWideBand(tc)
            % A pure fraction-of-bandwidth rule is catastrophic on a wide band: 10% of
            % [1 45] Hz would be 4.4 Hz of taper eating the low-frequency end.
            h = rheome.jointfilterbank.band(1, 45);
            tc.verifyGreaterThan(h(2*pi*2.5), 0.99);   % 1.5 Hz in, past a 1 Hz taper
        end

        function widthIsClampedToHalfTheBandwidth(tc)
            % The two edges must not overlap.
            h = rheome.jointfilterbank.band(10, 12, 50);
            v = h(2*pi*11);
            tc.verifyGreaterThan(v, 0);
            tc.verifyLessThanOrEqual(v, 1);
        end

        function gainsStayInUnitInterval(tc)
            h = rheome.jointfilterbank.band(8, 13);
            v = h(2*pi*linspace(0, 40, 501));
            tc.verifyTrue(all(v >= 0 & v <= 1));
        end

        function decreasingBandIsRefused(tc)
            tc.verifyError(@() rheome.jointfilterbank.band(13, 8), 'jointfilterbank:band');
        end

        function shapeFollowsTheInputShape(tc)
            h = rheome.jointfilterbank.band(8, 13);
            tc.verifySize(h(2*pi*(1:10)),  [1 10]);
            tc.verifySize(h(2*pi*(1:10)'), [10 1]);
        end

    end
end

% Author: Diellor Basha, 2026
