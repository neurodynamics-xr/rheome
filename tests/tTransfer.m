classdef tTransfer < matlab.unittest.TestCase
% rheome.geom.hemispheres / sphereregister, rheome.eigen.transfer / transfercheck, rheome.inverse.hodge_register on a synthetic
% two-piece surface (two icospheres side by side, a medial cap cut on each, each its own registered sphere), registered onto
% itself: what must hold exactly. The parity with the editor's stored registration is in tNxrStoreParity.
%
% Author: Diellor Basha, 2026

    properties
        side
    end

    methods (TestClassSetup)
        function build(tc)
            [V0, F0] = rheome.geom.icosphere(3); n0 = size(V0, 1);
            P = [V0 - [1.5 0 0]; V0 + [1.5 0 0]]; F = [F0; F0 + n0];
            wall = [V0(:, 1) > 0.75; -V0(:, 1) > 0.75];
            [X, G, keep] = rheome.geom.glue(P, F, wall);
            H = rheome.operators.hodge_dirac(X, G);
            pkg = rheome.eigen.hodge_package(rheome.eigen.hodge_dirac(H, 40), H);
            tc.side = struct('pkg', pkg, 'keep', keep, 'X', X, 'Pin', P, 'Fin', F, 'Sphere', [V0; V0], 'Left', "-x");
        end
    end

    methods (Test)
        function hemispheresByTheAxis(tc)
            h = rheome.geom.hemispheres(tc.side.Pin, tc.side.Fin, "-x");
            tc.verifyTrue(all(tc.side.Pin(h == 'L', 1) < 0) && all(tc.side.Pin(h == 'R', 1) > 0) && ~any(h == '?'));
            tc.verifyEqual(rheome.geom.hemispheres(tc.side.Pin, tc.side.Fin, "+x"), char('L' + 'R' - h), 'the other axis swaps them');
        end

        function registrationIsAPartitionOfUnity(tc)
            [T, R] = rheome.eigen.transfer(tc.side, tc.side);
            m = nnz(tc.side.keep > 0);
            tc.verifyEqual(full(sum(R, 2)), ones(size(R, 1), 1), 'AbsTol', 1e-12, 'every template vertex a blend of subject vertices');
            tc.verifyLessThan(full(max(abs(R(1:m, 1:m) - speye(m)), [], 'all')), 1e-12, 'cortex onto itself: the identity');
            tc.verifyEqual(full(max(abs(R(:, m+1:end)), [], 'all')), 0, 'the subject''s bridge is never sampled');
            k = tc.side.pkg.kind(:);
            tc.verifyEqual(max(abs(T(k == 0, k == 2)), [], 'all') + max(abs(T(k == 2, k == 0)), [], 'all'), 0, 'block-diagonal by kind');
        end

        function registeringWeightsRegistersTheMaps(tc)
            [T, R] = rheome.eigen.transfer(tc.side, tc.side);
            pkg = tc.side.pkg; K = numel(pkg.mu);
            B = randn(RandStream('mt19937ar', 'Seed', 5), 2 * K, 6); x = randn(6, 3);
            Bt = rheome.inverse.hodge_register(B, T, pkg.mu, pkg.mu); w = B * x;
            tc.verifyEqual(Bt * x, [T * w(1:K, :); (T * (pkg.mu .* w(K+1:end, :))) ./ pkg.mu], 'AbsTol', 1e-10);
            tc.verifyEqual(rheome.inverse.hodge_register(B, T), [T * B(1:K, :); T * B(K+1:end, :)], 'AbsTol', 1e-12);
            tc.verifyError(@() rheome.inverse.hodge_register(B, T, pkg.mu), 'rheome:hodge_register');
            C = rheome.eigen.transfercheck(T, R, tc.side, tc.side);
            tc.verifyEqual(height(C), 15); tc.verifyTrue(all(C.l == repelem((1:3)', 2 * (1:3)' + 1)));
        end
    end
end
