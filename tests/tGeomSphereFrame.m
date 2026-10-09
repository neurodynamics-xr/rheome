classdef tGeomSphereFrame < matlab.unittest.TestCase
% rheome.geom.sphereframe is the meridian frame with its poles at the sphere's +-Axis, pushed to any mesh
% sharing the sphere's vertices (ported from nxr-cortical-flow-matlab feat/sphere-pole-frame, 21cc8d9b).
%
% Author: Diellor Basha, 2026
    methods (Test)
        function onTheSphereItIsTheMeridian(t)
            [V, F] = rheome.geom.icosphere(3);
            fr = rheome.geom.sphereframe(V, F, V, V);
            mer = [0 0 1] - V(:,3).*V;  mer = mer ./ vecnorm(mer,2,2);
            ok = ~fr.singular;
            t.verifyLessThan(max(acosd(min(1, sum(fr.e1(ok,:).*mer(ok,:),2)))), 1e-4);
            t.verifyEqual(sort([fr.north fr.south]), sort([find(V(:,3)==max(V(:,3)),1) find(V(:,3)==min(V(:,3)),1)]));
            t.verifyLessThan(max(abs(sum(fr.e2(ok,:).*cross(V(ok,:), fr.e1(ok,:), 2), 2) - 1)), 1e-10);
        end
        function aScaledCortexGetsTheSameDirections(t)
            % an anisotropic stretch of the sphere: e1 must follow the push-forward of the meridian
            [V, F] = rheome.geom.icosphere(3);  A = diag([0.07 0.09 0.05]);
            Vc = V*A;  N = V/A;  N = N ./ vecnorm(N,2,2);        % normals of an ellipsoid
            fr = rheome.geom.sphereframe(Vc, F, N, V);
            mer = [0 0 1] - V(:,3).*V;  pf = mer*A;  pf = pf ./ vecnorm(pf,2,2);
            ok = ~fr.singular & abs(V(:,3)) < 0.9;
            t.verifyLessThan(median(acosd(min(1, sum(fr.e1(ok,:).*pf(ok,:),2)))), 1);
        end
    end
end

% Author: Diellor Basha, 2026
