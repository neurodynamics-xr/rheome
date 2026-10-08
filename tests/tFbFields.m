classdef tFbFields < matlab.unittest.TestCase
% The synthetic fields carry the demos' ground truth, so they are tested on their own.

    methods (Test)

        function bumpIsExactlyGaussianInGEODESICDistance(tc)
            % Assert the functional form at EVERY vertex rather than sampling the nearest
            % vertex to d = sigma: ico3 spacing is ~14 mm against sigma = 20 mm, so
            % "nearest" is several mm off and would test the mesh, not the bump.
            S = rheome.demos.fbd_selftest('sphere', 3, 60);
            sigma = 0.020;
            u = rheome.demos.fbd_selftest('bump', S, 1, sigma);
            d = S.R * acos(min(1, max(-1, S.V * S.V(1,:).' / S.R^2)));
            tc.verifyEqual(u, exp(-d.^2/(2*sigma^2)), 'AbsTol', 1e-12);
        end

        function bumpUsesGeodesicNotChordalDistance(tc)
            % The distinction matters on a sphere: at 90 degrees the chord is R*sqrt(2)
            % while the geodesic is R*pi/2, a 10% difference that would bias every
            % recovered width.
            S = rheome.demos.fbd_selftest('sphere', 3, 60);
            sigma = 0.050;
            u = rheome.demos.fbd_selftest('bump', S, 1, sigma);
            chord = vecnorm(S.V - S.V(1,:), 2, 2);
            uChord = exp(-chord.^2/(2*sigma^2));
            tc.verifyGreaterThan(max(abs(u - uChord)), 0.01);
        end

        function bumpPeaksAtItsSeed(tc)
            S = rheome.demos.fbd_selftest('sphere', 3, 60);
            u = rheome.demos.fbd_selftest('bump', S, 77, 0.02);
            [~, ipk] = max(u);
            tc.verifyEqual(ipk, 77);
        end

        function sectoralHarmonicLandsOnItsOwnEigenvalue(tc)
            % A degree-l harmonic on a sphere of radius R has lambda = l(l+1)/R^2. Its LBO
            % spectrum must concentrate there -- this is what makes the joint demo's ridge
            % analytic rather than fitted.
            S = rheome.demos.fbd_selftest('sphere', 4, 200);
            l = 5;
            u = rheome.demos.fbd_selftest('sectoral', S, l, 0);
            c = S.T.forward(u);
            [~, ipk] = max(abs(c));
            tc.verifyEqual(S.Lambda(ipk), l*(l+1)/S.R^2, 'RelTol', 0.08);
        end

        function rotatingFieldCarriesItsAnalyticRidge(tc)
            S = rheome.demos.fbd_selftest('sphere', 3, 60);
            Om = 2*pi*2;  tv = (0:255)/256;
            [U, gt] = rheome.demos.fbd_selftest('rotating', S, 3:6, Om, tv);
            tc.verifySize(U, [S.nV, 256]);
            tc.verifyEqual(gt.omega, (3:6)*Om, 'RelTol', 1e-12);
            tc.verifyEqual(gt.k, sqrt((3:6).*(4:7))/S.R, 'RelTol', 1e-12);
            tc.verifyTrue(isfinite(gt.slope) && gt.slope > 0);
        end

        function rotatingFieldActuallyRotates(tc)
            % A rigid rotation must preserve the field's norm over time.
            S = rheome.demos.fbd_selftest('sphere', 3, 60);
            tv = (0:63)/64;
            U = rheome.demos.fbd_selftest('rotating', S, 4, 2*pi, tv);
            n = vecnorm(U, 2, 1);
            tc.verifyEqual(max(n)/min(n), 1, 'RelTol', 0.05);
        end

    end
end

% Author: Diellor Basha, 2026
