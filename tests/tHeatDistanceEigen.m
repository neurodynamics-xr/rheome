classdef tHeatDistanceEigen < matlab.unittest.TestCase
% rheome.geom.measure.heatDistanceEigen without opts.t: t is estimated from the mean edge length, read from fg's
% .Vertices when it has them, else from face_gradient's operators (the same h to rounding); with neither, a clear
% error (rheome:heatDistanceEigen:noT) asks for opts.t.
%
% Author: Diellor Basha, 2026

    properties
        V; F; basis; fg
    end

    methods (TestClassSetup)
        function sphere(tc)
            [tc.V, tc.F] = rheome.geom.icosphere(3);
            [L, M] = rheome.operators.laplace_beltrami(tc.V, tc.F);
            [Phi, Lam] = eigs(-L, M, 200, 'smallestabs');  lam = abs(diag(Lam));  [lam, o] = sort(lam);
            tc.basis = struct('Phi', Phi(:, o), 'Lambda', lam, 'Mass', M);
            tc.fg = rheome.operators.face_gradient(tc.V, tc.F);
        end
    end

    methods (Test)
        function withoutTFromTheFaceGradients(tc)
            h = mean(vecnorm(tc.V(tc.F(:,1),:) - tc.V(tc.F(:,2),:), 2, 2));
            want = rheome.geom.measure.heatDistanceEigen(tc.basis, tc.fg, [1 7], struct('t', max(h^2, 50 / max(tc.basis.Lambda))));
            got = rheome.geom.measure.heatDistanceEigen(tc.basis, tc.fg, [1 7]);           % fg has no .Vertices
            tc.verifyEqual(got, want, 'AbsTol', 1e-12);
            fgV = tc.fg;  fgV.Vertices = tc.V;                                             % with .Vertices: the same
            tc.verifyEqual(rheome.geom.measure.heatDistanceEigen(tc.basis, fgV, [1 7]), want, 'AbsTol', 1e-12);
        end
        function hFromGradientsIsTheEdgeLength(tc)
            % eigenvalues scaled up so that 50/max(Lambda) < h^2: t = h^2, so h itself is what is compared
            b = tc.basis;  b.Lambda = b.Lambda * 1e5;
            h = mean(vecnorm(tc.V(tc.F(:,1),:) - tc.V(tc.F(:,2),:), 2, 2));
            tc.assertGreaterThan(h^2, 50 / max(b.Lambda));
            t = h^2;
            fgV = tc.fg;  fgV.Vertices = tc.V;
            tc.verifyEqual(rheome.geom.measure.heatDistanceEigen(b, tc.fg, 3), ...
                rheome.geom.measure.heatDistanceEigen(b, fgV, 3, struct('t', t)), 'AbsTol', 1e-12);
        end
        function noWayToEstimateTIsAClearError(tc)
            bare = struct('Gx', tc.fg.Gx, 'Gy', tc.fg.Gy, 'Gz', tc.fg.Gz, 'FaceArea', tc.fg.FaceArea);   % no .Faces
            tc.verifyError(@() rheome.geom.measure.heatDistanceEigen(tc.basis, bare, 1), 'rheome:heatDistanceEigen:noT');
            noGrad = struct('Faces', tc.F, 'Gx', tc.fg.Gx, 'Gy', tc.fg.Gy, 'Gz', tc.fg.Gz);              % no .FaceArea
            tc.verifyError(@() rheome.geom.measure.heatDistanceEigen(tc.basis, noGrad, 1), 'rheome:heatDistanceEigen:noT');
            tc.verifyWarningFree(@() rheome.geom.measure.heatDistanceEigen(tc.basis, bare, 1, struct('t', 1e-3)));
        end
    end
end
