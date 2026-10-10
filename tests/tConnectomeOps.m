classdef tConnectomeOps < matlab.unittest.TestCase
% MD16, the connectome operators, against closed forms on two ico3 spheres joined homotopically
% (W_ij = a_i between the same vertex of the two spheres). connectome_validate.m is the full study.
%
% Author: Diellor Basha, 2026

    properties
        S;  W;  a;  n;  R = 0.1;
    end

    methods (TestClassSetup)
        function twoSpheres(tc)
            [Vu, F] = rheome.geom.icosphere(3);  Vu = Vu ./ vecnorm(Vu,2,2);  m = size(Vu,1);
            [~, M1] = rheome.operators.laplace_beltrami(tc.R*Vu, F, 'galerkin');
            tc.n = m;  tc.a = full(sum(M1,2));
            tc.S = struct('Vertices',[tc.R*Vu; tc.R*Vu + [3*tc.R 0 0]], 'Faces',[F; F+m], 'nV',2*m, ...
                'nF',2*size(F,1), 'VertNormals',[Vu; Vu], 'Comment','two spheres', ...
                'Hemi',{{(1:m)', (m+1:2*m)'}}, 'HemiLabel',{{'Cortex L','Cortex R'}}, 'SurfaceFile','');
            tc.W = sparse([1:m, m+1:2*m], [m+1:2*m, 1:m], [tc.a; tc.a], 2*m, 2*m);
        end
    end

    methods (Test)
        function endpointsSnapToVertices(tc)
            % SmoothHops = 0: W is exactly the fibre count between the snapped end vertices.
            V = tc.S.Vertices;  m = tc.n;  ep = zeros(m, 2, 3);
            ep(:,1,:) = V(1:m,:) + 1e-6;  ep(:,2,:) = V(m+1:end,:) - 1e-6;
            [Wc, keep] = rheome.operators.connectome(V, tc.S.Faces, ep, struct('SmoothHops',0, 'EdgeThr',0));
            tc.verifyEqual(full(Wc), full(double(tc.W > 0)));
            tc.verifyNumElements(keep, 2);          % homotopic pairs only: every component is a pair
        end

        function smoothingConnectsTheComponent(tc)
            V = tc.S.Vertices;  m = tc.n;  ep = zeros(m, 2, 3);
            ep(:,1,:) = V(1:m,:);  ep(:,2,:) = V(m+1:end,:);
            [Wc, keep] = rheome.operators.connectome(V, tc.S.Faces, ep);
            tc.verifyEqual(Wc, Wc', 'AbsTol', 1e-12);
            tc.verifyEqual(full(diag(Wc)), zeros(2*m,1));
            tc.verifyNumElements(keep, 2*m);
        end

        function normalizedLaplacianTrick(tc)
            % rheome.eigen.connectome takes the largest of N; lambda_L = 1 - mu must equal eig(I - N).
            Wc = tc.W + 0.01*(tc.S.Vertices*tc.S.Vertices' > 0) .* ~speye(2*tc.n);   % connect it
            [A, B, N] = rheome.operators.connectome_laplacian(Wc, (1:2*tc.n)');
            tc.verifyEqual(full(A), full(B - N), 'AbsTol', 1e-14);
            b = rheome.eigen.connectome(N, 12);
            e = sort(eig(full(A)));
            tc.verifyEqual(b.Lambda, e(1:12), 'AbsTol', 1e-8);
        end

        function antisymmetricBranchShiftsByTwoGamma(tc)
            % Closed form: symmetric modes are blind to the connectome; the l = 0 antisymmetric mode
            % (+1 on one sphere, -1 on the other) sits at exactly 2*gamma (rowsum M = a).
            g = 3/(2*tc.R^2);  m = tc.n;
            [A, B, gOut] = rheome.operators.lb_connectome(tc.S, tc.W, (1:2*m)', g);
            tc.verifyEqual(gOut, g);
            % 24 modes cover degrees 0-2 of both branches; compare those only (eigs can miss part of a
            % degenerate cluster when K cuts at it)
            b = rheome.eigen.modes(A, B, 24);  P = b.Phi;
            par = sum(P(1:m,:).*P(m+1:end,:)) ./ (0.5*sum(P.^2));
            tc.verifyEqual(abs(par), ones(1,24), 'AbsTol', 1e-6);
            la = b.Lambda(par < 0);
            tc.verifyEqual(la(1), 2*g, 'RelTol', 1e-8);
            [A0, B0] = rheome.operators.lb_connectome(tc.S, tc.W, (1:2*m)', 0);
            % decoupled spheres double every eigenvalue, and eigs may return only some copies of a
            % cluster, so match each symmetric eigenvalue to its nearest LB eigenvalue
            b0 = rheome.eigen.modes(A0, B0, 24);  ls = b.Lambda(par > 0);
            tc.verifyLessThan(min(abs(ls(1:9) - b0.Lambda'), [], 2) ./ max(ls(1:9), 1), 1e-8);
        end

        function autoGammaMatchesDiagonals(tc)
            [A, ~, g] = rheome.operators.lb_connectome(tc.S, tc.W, (1:2*tc.n)');
            Lc = spdiags([tc.a; tc.a], 0, 2*tc.n, 2*tc.n) - tc.W;  K = A - g*Lc;
            tc.verifyEqual(g, full(mean(abs(diag(K)))) / mean([tc.a; tc.a]), 'RelTol', 1e-12);
        end

        function liftTriplesEachMode(tc)
            [A, B] = rheome.operators.lb_connectome(tc.S, tc.W, (1:2*tc.n)', 1/tc.R^2);
            b = rheome.eigen.modes(A, B, 6);  q = rheome.eigen.lift(b.Phi, b.Lambda, B);
            tc.verifySize(q.Phi, [4*2*tc.n, 18]);
            tc.verifyEqual(q.Lambda, repelem(b.Lambda(:), 3), 'AbsTol', 0);
            tc.verifyEqual(full(q.Phi'*q.Mass*q.Phi), eye(18), 'AbsTol', 1e-8);
            tc.verifyEqual(q.Phi(1:4:end,:), zeros(2*tc.n,18));       % real (w) slot stays empty
        end
    end
end

% Author: Diellor Basha, 2026
