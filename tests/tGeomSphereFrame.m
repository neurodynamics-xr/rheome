classdef tGeomSphereFrame < matlab.unittest.TestCase
% rheome.geom.sphereframe is the meridian frame with its poles at the sphere's +-Axis, pushed to any mesh
% sharing the sphere's vertices; rheome.geom.spherepatches says where pooling in it breaks; and
% rheome.scale.reducegaugeflow pools subjects in the shared frame, switching charts at the poles.
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
        function patchesBreakOnlyNearThePolesAndTheOtherChartCoversThem(t)
            [V, ~] = rheome.geom.icosphere(5);
            Pz = rheome.geom.spherepatches(V, Level=3);  Px = rheome.geom.spherepatches(V, Level=3, Axis=[1 0 0]);
            t.verifyEqual(Pz.patch, Px.patch, 'the patching does not depend on the chart');
            t.verifyEqual(sum(Pz.hasPole), 2);
            t.verifyGreaterThan(min(abs(90 - Pz.colat(Pz.excluded))), 60, 'only polar caps are excluded');
            t.verifyFalse(any(Pz.excluded & Px.excluded), 'chart x must cover every z-excluded patch');
        end
        function reduceReadsHeadingAndAgreement(t)
            d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
            % two subjects: patch 1 flows north in both (R = 1); patch 2 north vs south (R = 0, axial)
            for s = 1:2
                sg = 3 - 2*s;                                   % +1, -1
                T = table(["z";"z";"x";"x"], repmat("periodic",4,1), [1;2;1;2], [5;5;5;5], [10;10;10;10], ...
                    [0.01; sg*0.01; 0; 0], [0;0;0;0], [1e-4;1e-4;0;0], [0;0;0;0], [0;0;0;0], [1;1;1;1], ...
                    [60;60;90;90], [0;10;0;10], [1;1;1;1], [0;0;0;0], [0;0;0;0], 'VariableNames', ...
                    {'chart','band','patch','n_vertices','n_samples','vn_mean','vw_mean','tnn','tnw','tww', ...
                     'env_mean','colat_deg','lon_deg','spread_deg','has_pole','excluded'});
                T = [table(repmat("s"+s,4,1), repmat("norm",4,1), repmat("omega",4,1), ...
                     'VariableNames', {'subject','cohort','dataset'}), T];
                writetable(T, fullfile(d, "s"+s+"_gaugeflow.csv"));
            end
            G = rheome.scale.reducegaugeflow(d, fullfile(d, 'out'));
            t.verifyEqual(G.n_subjects, [2; 2]);
            t.verifyEqual(G.chart, ["z"; "z"]);
            t.verifyEqual(G.heading_deg(1), 0, 'AbsTol', 1e-9);
            t.verifyEqual(G.R, [1; 0], 'AbsTol', 1e-9);
            t.verifyEqual(G.axis_deg, [0; 0], 'AbsTol', 1e-9);
            t.verifyEqual(G.anisotropy, [1; 1], 'AbsTol', 1e-9);
        end
    end
end

% Author: Diellor Basha, 2026
