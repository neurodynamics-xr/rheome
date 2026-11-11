classdef tOneDelaunay < matlab.unittest.TestCase
% One mesh and mass (PI, 2026-11-10): the intrinsic, mollified Delaunay triangulation with the mixed-Voronoi mass is the
% default of the Laplace-Beltrami operator, the connection Laplacian and the heat method; the embedded mesh and the
% Galerkin / barycentric masses are named options that reproduce the previous operators.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function laplaceBeltramiDefaultIsTheReferenceDelaunay(tc)
            [V, F] = i_mesh();
            [L, M, info] = rheome.operators.laplace_beltrami(V, F);
            S = rheome.operators.intrinsic_delaunay(V, F);
            tc.verifyEqual(full(L), full((S.L + S.L')/2), AbsTol=1e-12);
            tc.verifyEqual(full(diag(M)), S.mass, AbsTol=1e-14);
            tc.verifyEqual(nnz(M - spdiags(diag(M), 0, size(M,1), size(M,1))), 0);
            tc.verifyEqual([info.mesh info.mass], ["intrinsic" "mixedVoronoi"]);
            tc.verifyEqual([info.mollify info.flipEps info.flips], [1e-6 1e-6 S.flips]);
            tc.verifyGreaterThan(S.flips, 0);                    % the test mesh is not Delaunay: this is not the embedded operator
            tc.verifyNotEqual(norm(L - i_old_cotan(V, F), 'fro'), 0);
        end
        function embeddedGalerkinReproducesThePreviousOperator(tc)
            [V, F] = i_mesh();
            [L, M, info] = rheome.operators.laplace_beltrami(V, F, mesh='embedded', mass='galerkin');
            tc.verifyEqual(full(L), full(i_old_cotan(V, F)), AbsTol=1e-13);
            tc.verifyEqual(full(M), full(rheome.operators.mass(V, F, 'galerkin')), AbsTol=1e-15);
            tc.verifyEqual(full(sum(M(:))), i_area(V, F), RelTol=1e-12);
            tc.verifyEqual(info.mesh, "embedded");
            [~, B] = rheome.operators.laplace_beltrami(V, F, mass='barycentric');
            tc.verifyEqual(full(diag(B)), accumarray(F(:), repmat(i_faceArea(V, F)/3, 3, 1)), AbsTol=1e-14);
        end
        function optionsArePassedThrough(tc)
            [V, F] = i_mesh();
            [~, ~, a] = rheome.operators.laplace_beltrami(V, F, mollify=1e-3, flipEps=1e-2);
            S = rheome.operators.intrinsic_delaunay(V, F, mollify=1e-3, flipEps=1e-2);
            tc.verifyEqual([a.mollify a.flipEps a.flips a.mollifyEps], [1e-3 1e-2 S.flips S.mollifyEps]);
            tc.verifyEqual([S.mollify S.flipEps], [1e-3 1e-2]);
        end
        function massTypes(tc)
            [V, F] = i_mesh();
            S = rheome.operators.intrinsic_delaunay(V, F);
            tc.verifyEqual(full(diag(rheome.operators.mass(V, F))), S.mass, AbsTol=1e-14);   % default = mixed Voronoi
            for t = ["mixedVoronoi" "galerkin" "barycentric"]
                tc.verifyEqual(full(sum(rheome.operators.mass(V, F, t), 'all')), i_area(V, F), RelTol=1e-9);
            end
            tc.verifyError(@() rheome.operators.mass(V, F, 'lumped'), 'operators:mass:type');
        end
        function connectionLaplacianDefaultIsTheIntrinsicConstruction(tc)
            [V, F] = i_mesh();
            C = rheome.operators.connection_laplacian(V, F);
            [L, M, S] = rheome.operators.connection_laplacian_idt(V, F);
            tc.verifyEqual(C.mesh, "intrinsic");
            tc.verifyEqual(C.flips, S.flips);
            nV = size(V, 1);
            tc.verifyEqual(norm(C.A - C.A', 'fro'), 0, AbsTol=1e-12);
            % the complex A is the real 2n form: same spectrum, each eigenvalue twice
            ev = sort(real(eig(full(C.A)))); el = sort(eig(full(L)));
            tc.verifyLessThan(max(abs(el(1:2:end) - ev)), 1e-9 * max(1, max(ev)));
            tc.verifyEqual(full(diag(C.B)), M(1:2:end), AbsTol=1e-14);
            tc.verifyEqual(full(diag(C.B)), S.mass, AbsTol=1e-14);
            tc.verifyEqual(size(C.Rt), [nV nV]);
            tc.verifyEqual(abs(nonzeros(C.Rt)), ones(nnz(C.Rt), 1), AbsTol=1e-12);
            tc.verifyEqual(full(C.Rt), full(C.Rt')', AbsTol=1e-12);
            tc.verifyEqual(full(C.Rt .* C.Rt.'), full(spones(C.Rt)), AbsTol=1e-12);          % transport back is the inverse
            E = rheome.operators.connection_laplacian(V, F, mesh='embedded', mass='galerkin');
            tc.verifyEqual(E.mesh, "embedded");
            tc.verifyEqual(full(sum(E.B(:))), i_area(V, F), RelTol=1e-12);
        end
        function heatDistanceUsesTheReferenceDelaunay(tc)
            [V, F] = i_mesh();
            [d, solver] = rheome.geom.measure.heatDistance(struct('Vertices', V, 'Faces', F), 1);
            S = rheome.operators.intrinsic_delaunay(V, F);
            tc.verifyEqual(solver.mesh.flips, S.flips);
            tc.verifyEqual(solver.mesh.len, S.len, AbsTol=1e-14);
            tc.verifyEqual(solver.mesh.mass, S.mass, AbsTol=1e-14);
            tc.verifyEqual([solver.mollify solver.flipEps], [1e-6 1e-6]);
            tc.verifyEqual(d(1), 0, AbsTol=1e-12);
            tc.verifyTrue(all(isfinite(d)) && all(d >= -1e-9));
        end
        function diracFrameMassDefaultIsMixedVoronoi(tc)
            [V, F] = rheome.geom.icosphere(1);
            [~, B] = rheome.operators.dirac_frame(V, F, 0.5);
            S = rheome.operators.intrinsic_delaunay(V, F);
            tc.verifyEqual(full(diag(B)), repelem(S.mass, 4), AbsTol=1e-14);
            [~, Bg] = rheome.operators.dirac_frame(V, F, 0.5, [], true, mass='galerkin');
            tc.verifyEqual(full(Bg), full(kron(rheome.operators.mass(V, F, 'galerkin'), speye(4))), AbsTol=1e-15);
        end
    end
end

function [V, F] = i_mesh()
% a sphere with its vertices jittered radially: obtuse triangles, so the intrinsic Delaunay flips some edges
[V, F] = rheome.geom.icosphere(2);
rng(7);  V = V .* (1 + 0.25*(rand(size(V,1), 1) - 0.5)) .* [1 1.3 0.8];
end
function a = i_faceArea(V, F)
a = 0.5 * vecnorm(cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2), 2, 2);
end
function a = i_area(V, F), a = sum(i_faceArea(V, F)); end
function L = i_old_cotan(V, F)
% the embedded cotan stiffness as it was before 2026-11-10
nV = size(V, 1);  W = sparse(nV, nV);
for k = 1:3
    a = F(:,k); b = F(:,mod(k,3)+1); c = F(:,mod(k+1,3)+1);          % cot at c weights edge (a,b)
    u = V(a,:) - V(c,:); v = V(b,:) - V(c,:);
    cot = sum(u.*v, 2) ./ vecnorm(cross(u, v, 2), 2, 2);
    W = W + sparse([a; b], [b; a], 0.5*[cot; cot], nV, nV);
end
L = spdiags(full(sum(W, 2)), 0, nV, nV) - W;
end
