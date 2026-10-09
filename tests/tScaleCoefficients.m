classdef tScaleCoefficients < matlab.unittest.TestCase
% TSCALECOEFFICIENTS  rheome.scale.coefficients on two synthetic sphere "hemispheres" (no cached data):
% the Prognome contract (shapes, ladder codes), a planted source found in its own tile and hemisphere,
% a flat analytic envelope, linearity, and the refusals.
%
% Author: Diellor Basha, 2026

    properties
        B
    end

    methods (TestClassSetup)
        function makeBases(tc)
            [V, Fc] = rheome.geom.icosphere(3);  V = 0.05 * V;  nV = size(V, 1);   % 642 vertices, 5 cm
            [L, M] = rheome.operators.laplace_beltrami(V, Fc);
            [Phi, D] = eigs(L, M, 80, 'smallestabs');
            [lam, o] = sort(max(real(diag(D)), 0));  Phi = real(Phi(:, o));
            S = struct('Vertices', V, 'Faces', Fc, 'nV', nV);
            lbo = struct('L', L, 'M', M, 'Phi', Phi, 'Lambda', lam, 'Mass', M);
            tc.B.L = struct('gv', (1:nV)', 'S', S, 'lbo', lbo);
            tc.B.R = struct('gv', nV + (1:nV)', 'S', S, 'lbo', lbo);
        end
    end

    methods (Test)
        function contractShapesAndCodes(tc)
            [K, F] = tc.planted(5, 1);
            C = rheome.scale.coefficients(tc.B, K, F, 200, Depth=3, Band=[5 20], Rate=50);
            nS = numel(C.scales);
            tc.verifyGreaterThan(nS, 1);
            tc.verifyClass(C.W, 'single');
            tc.verifySize(C.W, [16 nS 200]);                       % 2 x 2^3 tiles, 4 s at 50 Hz
            tc.verifyEqual(C.codes, (16:31)');                      % one full ladder level under the root
            tc.verifyEqual(bitshift(C.codes, -3), [2 * ones(8, 1); 3 * ones(8, 1)]);   % L under 2, R under 3
            tc.verifyEqual(C.fs, 50);
            tc.verifyEqual(C.depth, 4);
            tc.verifyTrue(all(isfinite(C.scales) & C.scales > 0), 'wavelet scales only, no lowpass NaN');
            tc.verifyTrue(all(diff(C.scales) > 0));
        end

        function plantedSourceLandsInItsTile(tc)
            v = 5;  [K, F] = tc.planted(v, 1);
            C = rheome.scale.coefficients(tc.B, K, F, 200, Depth=3, Band=[5 20], Rate=50);
            T = rheome.geom.tree(tc.B.L.S, L=tc.B.L.lbo.L, M=tc.B.L.lbo.M, MaxDepth=3);
            G = rheome.geom.tiles(T, tc.B.L.S, 3);
            want = G.node_id(G.tileOf(v)) + 2^3;                    % its ladder code
            mid = 51:150;                                           % away from the filter's edges
            E = median(C.W(:, :, mid), 3);                          % [tile x scale]
            [~, imax] = max(E, [], 1);
            tc.verifyEqual(C.codes(imax(1)), want, 'finest scale: the source tile');
            tc.verifyLessThan(max(E(C.hemisphere == 2, :), [], 'all'), 1e-9 * max(E, [], 'all'), ...
                'a left-hemisphere source leaves the right hemisphere empty');
            e = squeeze(C.W(imax(1), 1, mid));
            tc.verifyLessThan(std(e) / mean(e), 0.05, 'the envelope of a steady 10 Hz source is flat');
        end

        function linearInTheData(tc)
            [K, F] = tc.planted(40, 1);
            C1 = rheome.scale.coefficients(tc.B, K, F, 200, Depth=2, Band=[5 20], Rate=50);
            C2 = rheome.scale.coefficients(tc.B, K, 2 * F, 200, Depth=2, Band=[5 20], Rate=50);
            tc.verifyEqual(C2.W, 2 * C1.W, 'RelTol', 1e-4);
        end

        function refusals(tc)
            [K, F] = tc.planted(5, 1);
            tc.verifyError(@() rheome.scale.coefficients(tc.B, K, F, 200, Depth=12, Band=[5 20], Rate=50), ...
                'scale:coefficients:ladder');
            tc.verifyError(@() rheome.scale.coefficients(tc.B, K, F, 200, Rate=30, Band=[5 10]), ...
                'scale:coefficients:rate');
            tc.verifyError(@() rheome.scale.coefficients(tc.B, K, F, 200, Rate=50, Band=[5 30]), ...
                'scale:coefficients:band');
            tc.verifyError(@() rheome.scale.coefficients(tc.B, K, F(1:end-1, :), 200, Depth=2, Band=[5 20], Rate=50), ...
                'scale:coefficients:channels');
        end
    end

    methods
        function [K, F] = planted(tc, v, amp)
            % identity kernel (channels = the 3 nV current components); a 10 Hz x-current at left vertex v
            nV = numel(tc.B.L.gv) + numel(tc.B.R.gv);
            K = speye(3 * nV);  t = (0:799) / 200;
            F = sparse(3 * nV, 800);  F(3 * (v - 1) + 1, :) = amp * sin(2 * pi * 10 * t);
            K = full(K);  F = full(F);
        end
    end
end

% Author: Diellor Basha, 2026
