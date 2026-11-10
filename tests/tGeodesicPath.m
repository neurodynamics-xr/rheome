classdef tGeodesicPath < matlab.unittest.TestCase
% rheome.geom.draw.line (FlipOut, rheome.geom.draw.flipOut): a great circle on a sphere (its length and its points on
% the surface), and the reused solver reversed. Point-for-point parity with the viewer's FlipOut is tFlipOut.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function aGreatCircleOnTheSphere(t)
            [V, F] = rheome.geom.icosphere(5);  S = struct('Vertices', V, 'Faces', F);
            a = 1;  [~, b] = min(V * V(a,:)' + 0.3*abs(V(:,3) - V(a,3)));     % far from a, not antipodal
            [P, info] = rheome.geom.draw.line(S, a, b);
            truth = acos(V(a,:)*V(b,:)');
            t.verifyEqual(info.length, truth, 'RelTol', 0.02);
            t.verifyGreaterThan(info.length, truth*(1 - 0.005), 'a path on the mesh is no shorter than the arc');
            t.verifyEqual(P([1 end],:), V([a b],:), 'AbsTol', 1e-12);
            n = cross(V(F(info.faces,2),:) - V(F(info.faces,1),:), V(F(info.faces,3),:) - V(F(info.faces,1),:), 2);
            off = abs(sum((P(1:end-1,:) - V(F(info.faces,1),:)).*n, 2))./vecnorm(n, 2, 2);
            t.verifyLessThan(max(off), 1e-10, 'every point lies in its face');
            [~, info2] = rheome.geom.draw.line(S, b, a, info.mesh);     % the reused flip mesh, reversed
            t.verifyEqual(info2.length, info.length, 'RelTol', 0.02);
        end
        function aHeatSolverIsIgnoredAndWarnsOnce(t)
            % the old call: [path, info] = geodesicpath(S, A, B, solver); solver = info.solver (reused for heatDistance)
            [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', V, 'Faces', F);
            [d, solver] = rheome.geom.measure.heatDistance(S, 1);
            clear functions %#ok<CLFUNC>   line's persistent "warned"
            [P, info] = t.verifyWarning(@() rheome.geom.draw.line(S, 1, 200, solver), 'rheome:deprecated');
            t.verifyWarningFree(@() rheome.geom.draw.line(S, 1, 200, solver));          % once per session
            t.verifyEqual(info.solver, solver);                                          % handed back unchanged
            [P0, i0] = rheome.geom.draw.line(S, 1, 200);
            t.verifyEqual(P, P0);  t.verifyEqual(info.length, i0.length);
            t.verifyEqual(rheome.geom.measure.heatDistance(S, 1, info.solver), d);      % the caller's chain still works
            clear functions %#ok<CLFUNC>
            ws = warning('error', 'rheome:deprecated');  c = onCleanup(@() warning(ws));  msg = '';
            try, rheome.geom.draw.line(S, 1, 200, solver); catch e, msg = e.message; end
            t.verifyTrue(contains(msg, 'heat distance') && contains(msg, 'flipMesh'), msg);
            warning(ws);  clear functions %#ok<CLFUNC>
            [Pg, ig] = t.verifyWarning(@() rheome.geom.geodesicpath(S, 1, 200, solver), 'rheome:deprecated');
            t.verifyEqual(Pg, P0);  t.verifyEqual(ig.solver, solver);
        end
    end
end
