classdef tGlueHodgeParity < matlab.unittest.TestCase
% rheome.geom.glue, rheome.operators.intrinsic_delaunay / hodge_dirac and rheome.eigen.hodge_dirac against the editor and nxr-convert:
% the editor's golden files (nxr-convert tests/fixtures/editor-parity: a two-piece mesh with a medial wall on each, written by the
% editor's own code) and nxr-convert's solve_hodge on it (tests/fixtures/hodge-solve, written by nxr-convert). The same checks as
% nxr-convert's test_eigen_parity.py. Skipped when the editor fixtures are absent.
%
% Author: Diellor Basha, 2026

    properties
        G   % the editor-parity folder
        P, F, wall
    end

    methods (TestClassSetup)
        function load(tc)
            tc.G = fullfile(getenv('HOME'), 'workspace/research/code/nxr-cortical-flow/backend/convert/nxr-convert/tests/fixtures/editor-parity');
            tc.assumeTrue(isfile(fullfile(tc.G, 'mesh.wall.u8')), 'the editor-parity fixtures are absent');
            tc.P = reshape(i_read(fullfile(tc.G, 'mesh.P.f64'), 'double'), 3, [])';
            tc.F = reshape(i_read(fullfile(tc.G, 'mesh.F.i32'), 'int32'), 3, [])' + 1;
            tc.wall = i_read(fullfile(tc.G, 'mesh.wall.u8'), 'uint8') > 0;
        end
    end

    methods (Test)
        function glueIsTheEditors(tc)
            meta = jsondecode(fileread(fullfile(tc.G, 'meta.json')));
            [X, Gf, keep, bridge] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            tc.verifyEqual(reshape((Gf - 1)', [], 1), i_read(fullfile(tc.G, 'glued.F.i32'), 'int32'), 'the same triangles');
            tc.verifyEqual(keep - 1, i_read(fullfile(tc.G, 'glued.keep.i32'), 'int32'), 'the same kept vertices');
            tc.verifyEqual(reshape(X', [], 1), i_read(fullfile(tc.G, 'glued.P.f64'), 'double'), 'AbsTol', 1e-15, 'the bridge''s points');
            tc.verifyEqual([bridge.vertices bridge.rings], [meta.bridge.vertices meta.bridge.rings]);
            tc.verifyEqual(rheome.geom.meshkey(Gf), meta.meshKey, 'the fingerprint the editor checks stored products by');
        end

        function glueDoesNotDependOnTheWinding(tc)
            [X1, F1, k1] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            [X2, F2, k2] = rheome.geom.glue(tc.P, tc.F(:,[1 3 2]), tc.wall);
            tc.verifyEqual(k1, k2); tc.verifyEqual(X1, X2);
            tc.verifyEqual(sortrows(sort(F1, 2)), sortrows(sort(F2, 2)));
        end

        function flipsAndPencilAreTheEditors(tc)
            meta = jsondecode(fileread(fullfile(tc.G, 'meta.json')));
            [X, Gf] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            H = rheome.operators.hodge_dirac(X, Gf);
            tc.verifyEqual(H.mesh.F - 1, i_read(fullfile(tc.G, 'intrinsic.F.i32'), 'int32'), 'the same flips');
            tc.verifyEqual([H.mesh.flips H.mesh.mollifyEps H.star1Clamped], [meta.flips meta.mollify 0]);
            ip = i_read(fullfile(tc.G, 'hodge.indptr.i32'), 'int32'); n = numel(ip) - 1;
            R = sparse(repelem((1:n)', diff(ip)), i_read(fullfile(tc.G, 'hodge.indices.i32'), 'int32') + 1, ...
                i_read(fullfile(tc.G, 'hodge.data.f64'), 'double'), n, n);
            tc.verifyEqual(H.N, meta.N);
            tc.verifyLessThanOrEqual(full(max(abs(H.K - R), [], 'all')), 1e-12 * full(max(abs(R), [], 'all')), 'K');
            tc.verifyEqual(H.M, i_read(fullfile(tc.G, 'hodge.mass.f64'), 'double'), 'RelTol', 1e-9, 'M');
        end

        function solveIsNxrConverts(tc)
            S = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'hodge-solve');
            meta = jsondecode(fileread(fullfile(S, 'meta.json')));
            [X, Gf] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            H = rheome.operators.hodge_dirac(X, Gf);
            sol = rheome.eigen.hodge_dirac(H, meta.n_modes);
            tc.verifyEqual([numel(sol.exactLambda) numel(sol.coexactLambda)], [meta.exact meta.coexact], 'the same band');
            tc.verifyEqual(sol.exactLambda, i_read(fullfile(S, 'exact_lambda.f64'), 'double'), 'RelTol', 1e-9);
            tc.verifyEqual(sol.coexactLambda, i_read(fullfile(S, 'coexact_lambda.f64'), 'double'), 'RelTol', 1e-9);
            tc.verifyLessThan(i_modegap(sol.exactModes, reshape(i_read(fullfile(S, 'exact_modes.f64'), 'double'), [], meta.exact), H.star0), 1e-6);
            tc.verifyLessThan(i_modegap(sol.coexactModes, reshape(i_read(fullfile(S, 'coexact_modes.f64'), 'double'), [], meta.coexact), H.star2), 1e-6);
        end

        function fusedWeightsReadBackTheirSynthesis(tc)
            % a kernel fused into the modes (carried onto the bridge), its weights synthesised into a current, the current read back:
            % the same weights (rheome.inverse.hodge_weights is the left inverse of rheome.inverse.hodge_synthesis)
            [X, Gf, keep] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            H = rheome.operators.hodge_dirac(X, Gf);
            pkg = rheome.eigen.hodge_package(rheome.eigen.hodge_dirac(H, 60), H);
            tc.verifyTrue(issorted(pkg.mu) && isa(pkg.P, 'single') && all(ismember(pkg.kind, [0 2])));
            W = randn(RandStream('mt19937ar', 'Seed', 3), 3 * size(tc.P, 1), 4);
            B = rheome.inverse.hodge_weights(W, pkg, X, Gf, keep, tc.P);
            K = numel(pkg.mu);
            tc.verifyEqual(B(pkg.kind == 0, :), zeros(nnz(pkg.kind == 0), 4), 'an exact mode carries no normal current');
            J = rheome.inverse.hodge_synthesis(B, pkg, X, Gf);
            back = rheome.inverse.hodge_weights(reshape(permute(J, [2 1 3]), [], 4), pkg, X, Gf);
            tc.verifyLessThan(norm(back - B, 'fro') / norm(B, 'fro'), 1e-10);
            tc.verifyEqual(size(B), [2 * K 4]);
        end

        function everyModeIsAnEigenpairOfThePencil(tc)
            [X, Gf] = rheome.geom.glue(tc.P, tc.F, tc.wall);
            H = rheome.operators.hodge_dirac(X, Gf);
            sol = rheome.eigen.hodge_dirac(H, 60);
            nV = numel(H.star0); nE = numel(H.star1); nF = numel(H.star2);
            for k = 1:numel(sol.exactLambda)                                        % (psi, d0 psi / mu, 0)/sqrt(2)
                mu = sqrt(sol.exactLambda(k)); p = sol.exactModes(:,k);
                x = [p; H.d0 * p / mu; zeros(nF,1)] / sqrt(2);
                tc.verifyLessThanOrEqual(norm(H.K * x - mu * H.M .* x), 1e-6 * mu * norm(H.M .* x));
            end
            for k = 1:numel(sol.coexactLambda)                                      % (0, delta chi / mu, chi)/sqrt(2)
                mu = sqrt(sol.coexactLambda(k)); c = sol.coexactModes(:,k);
                x = [zeros(nV,1); (H.d1' * (H.star2 .* c)) ./ H.star1 / mu; c] / sqrt(2);
                tc.verifyLessThanOrEqual(norm(H.K * x - mu * H.M .* x), 1e-6 * mu * norm(H.M .* x));
            end
            psi = sol.exactModes;
            tc.verifyLessThan(max(abs(psi' * (H.star0 .* psi) - eye(size(psi,2))), [], 'all'), 1e-8, 'star0-orthonormal');
            ne = numel(sol.exactLambda); nc = numel(sol.coexactLambda);
            tc.verifyTrue(abs(ne - nc) <= 3 && ne + nc >= 50, 'a complete band, half and half (Weyl)');
            tc.verifyTrue(issorted(sol.mu) && numel(sol.mu) == 2 * (ne + nc));
            tc.assertEqual(nE, size(H.d0, 1));
        end
    end
end

function g = i_modegap(A, B, w)
% the largest mass-norm difference between two sets of unit modes, each sign aligned
s = sign(sum(A .* w .* B, 1)); D = A - B .* s;
g = max(sqrt(sum(D .* w .* D, 1)));
end

function x = i_read(p, t)
fid = fopen(p); x = fread(fid, inf, [char(t) '=>double']); fclose(fid);       % integers too: exact in double
end
