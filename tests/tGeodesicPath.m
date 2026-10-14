classdef tGeodesicPath < matlab.unittest.TestCase
% rheome.geom.geodesicpath: a great circle on a sphere (its length and its points on the surface), and a path
% around a fold is not a chord. The fsaverage comparison against the exact geodesic is in request 07dbab95's evidence.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function aGreatCircleOnTheSphere(t)
            [V, F] = rheome.geom.icosphere(5);  S = struct('Vertices', V, 'Faces', F);
            a = 1;  [~, b] = min(V * V(a,:)' + 0.3*abs(V(:,3) - V(a,3)));     % far from a, not antipodal
            [P, info] = rheome.geom.geodesicpath(S, a, b);
            truth = acos(V(a,:)*V(b,:)');
            t.verifyEqual(info.length, truth, 'RelTol', 0.02);
            t.verifyGreaterThan(info.length, truth*(1 - 0.005), 'a path on the mesh is no shorter than the arc');
            t.verifyEqual(P([1 end],:), V([a b],:), 'AbsTol', 1e-12);
            n = cross(V(F(info.faces,2),:) - V(F(info.faces,1),:), V(F(info.faces,3),:) - V(F(info.faces,1),:), 2);
            off = abs(sum((P(1:end-1,:) - V(F(info.faces,1),:)).*n, 2))./vecnorm(n, 2, 2);
            t.verifyLessThan(max(off), 1e-10, 'every point lies in its face');
            [~, info2] = rheome.geom.geodesicpath(S, b, a, info.solver);   % the reused solver, reversed
            t.verifyEqual(info2.length, info.length, 'RelTol', 0.02);
        end
    end
end
