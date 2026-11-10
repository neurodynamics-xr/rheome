classdef tEigenRouting < matlab.unittest.TestCase
% The "smoothest" gauge and the tree's Fiedler vector solve through rheome.eigen.smallest, not a direct
% eigs(..., 'smallestabs') on a (near-)singular pencil. Synthetic meshes only.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function smallestSolvesWhereSmallestabsWasPoor(tc)
            % two disconnected spheres: an exactly singular pencil with a doubled null mode. The direct
            % 'smallestabs' solve returns eigenpairs with a residual of order 1e-1; smallest, 1e-8 or better.
            [V, F] = rheome.geom.icosphere(3);  V = V .* [1 1.3 2];
            [A, B] = rheome.operators.laplace_beltrami(V, F);
            A = blkdiag(A, A);  B = blkdiag(B, B);
            res = @(v, d) arrayfun(@(i) norm(A*v(:,i) - B*v(:,i)*d(i)) / norm(B*v(:,i)), 1:numel(d));
            wr = warning('off', 'all');  cl = onCleanup(@() warning(wr));
            [v0, d0] = eigs(A, B, 4, 'smallestabs', struct('tol', 1e-6, 'maxit', 500));
            r0 = max(res(v0, real(diag(d0))));
            [v1, d1] = rheome.eigen.smallest(A, B, 4);
            r1 = max(res(v1, diag(d1)));
            tc.log(1, sprintf('residual before (smallestabs) %.2g, after (eigen.smallest) %.2g', r0, r1));
            tc.verifyLessThan(r1, 1e-8);
            tc.verifyGreaterThan(r0, 1e-4, 'the pencil no longer shows smallestabs failing; pick a harder one');
        end

        function gaugeSmoothestEqualsTheDirectSolveWhereThatWorks(tc)
            [V, F] = rheome.geom.icosphere(3);  V = V .* [1 1.3 2];
            C = rheome.operators.connection_laplacian(V, F);
            g = rheome.operators.gauge(V, F, Method="smoothest", Connection=C);
            wr = warning('off', 'all');  cl = onCleanup(@() warning(wr));
            [u, d] = eigs(C.A, C.B, 1, 'smallestabs', struct('tol', 1e-10));
            e1 = cos(angle(u)).*C.e1 + sin(angle(u)).*C.e2;
            % the same field up to a global rotation of the frame: compare the smoothness measure and
            % the relative angles between neighbours
            [ii, jj] = find(triu(C.A ~= 0, 1));
            ang = @(x) acosd(max(-1, min(1, sum(x(ii,:).*x(jj,:), 2))));
            tc.verifyEqual(ang(g.e1), ang(e1), AbsTol=1e-4);
            [w, dn] = rheome.eigen.smallest(C.A, C.B, 1);
            tc.verifyEqual(dn, d, AbsTol=1e-8);
            tc.verifyLessThan(norm(C.A*w - C.B*w*dn) / norm(C.B*w), 1e-8);
        end

        function treeFiedlerSplitIsTheSecondModeNotTheConstant(tc)
            % an ellipsoid has a simple Fiedler vector; the root split must equal the sign of the dense
            % generalized eigenvector of the second-smallest eigenvalue
            [V, F] = rheome.geom.icosphere(3);  V = V .* [1 1.3 2] * 0.07;
            S = struct('Vertices', V, 'Faces', F, 'nV', size(V,1), 'nF', size(F,1));
            [L, M] = rheome.operators.laplace_beltrami(V, F);
            T = rheome.geom.tree(S, L=L, M=M, MaxDepth=1);
            [W, D] = eig(full(L), full(M));  [~, o] = sort(real(diag(D)));
            f = real(W(:, o(2)));
            kids = T(T.parent_id == T.node_id(1), :);
            % vertices on the nodal plane (|f| ~ 0 by symmetry) have no definite sign: leave them out
            sure = abs(f) > 1e-6 * max(abs(f));
            side = zeros(size(f));  side(kids.members{1}) = 1;  side(kids.members{2}) = -1;
            agree = sign(f(sure)) .* side(sure);
            ok = all(agree == agree(1));            % one child is exactly the positive side
            tc.verifyTrue(ok, 'the root split is not the Fiedler sign partition');
        end
    end
end
