classdef tGeodesic < matlab.unittest.TestCase
% rheome.geom.geodesic (heat method) against the great-circle distance on a sphere: accurate at long range
% (the old version saturated beyond ~50 edges), the same in metres and millimetres, defined on meshes with
% zero-area triangles, and Inf on another connected component. The fsaverage comparison against
% geometry-central is in the evidence of request 07dbab95 (compare_geodesic.m).
%
% Author: Diellor Basha, 2026
    methods (Test)
        function matchesGreatCircleOutToTheAntipode(t)
            [V, F] = rheome.geom.icosphere(6);  R = 0.07;  V = R*V;        % 70 mm sphere, ~1.2 mm edges
            d = rheome.geom.geodesic(struct('Vertices', V, 'Faces', F), 1);
            truth = R*acos(min(1, max(-1, V*V(1,:)'/R^2)));
            far = truth > 0.02;                                            % beyond ~16 edges
            t.verifyLessThan(max(abs(d(far) - truth(far))./truth(far)), 0.03);
            t.verifyGreaterThan(max(d), 0.97*pi*R, 'no saturation: the antipode is reached');
        end
        function unitsDoNotMatter(t)
            [V, F] = rheome.geom.icosphere(4);  V = 0.07*V;
            dm  = rheome.geom.geodesic(struct('Vertices', V,      'Faces', F), 5);
            dmm = rheome.geom.geodesic(struct('Vertices', 1000*V, 'Faces', F), 5);
            t.verifyLessThan(max(abs(1000*dm - dmm)), 1e-6*max(dmm));
        end
        function zeroAreaTrianglesStillFactor(t)
            % split one face at its edge midpoint: two triangles, one with zero area
            [V, F] = rheome.geom.icosphere(4);
            f = F(1,:);  m = (V(f(1),:) + V(f(2),:))/2;  V(end+1,:) = m;  k = size(V,1);
            F(1,:) = [f(1) k f(3)];  F(end+1,:) = [k f(2) f(3)];
            F(end+1,:) = [f(1) f(2) k];                                    % degenerate (collinear)
            d = rheome.geom.geodesic(struct('Vertices', V, 'Faces', F), 10);
            t.verifyTrue(all(isfinite(d)));
            t.verifyGreaterThan(max(d), 0.9*pi);
        end
        function anotherComponentIsInfinitelyFar(t)
            [V, F] = rheome.geom.icosphere(3);  n = size(V,1);
            S = struct('Vertices', [V; V + [3 0 0]], 'Faces', [F; F + n]);
            [d, solver] = rheome.geom.geodesic(S, 1);
            t.verifyTrue(all(isinf(d(n+1:end))));
            t.verifyEqual(d(1), 0, 'AbsTol', 1e-12);
            d2 = rheome.geom.geodesic(S, n+1, solver);                      % the reused solver
            t.verifyTrue(all(isinf(d2(1:n))) && all(isfinite(d2(n+1:end))));
        end
    end
end
