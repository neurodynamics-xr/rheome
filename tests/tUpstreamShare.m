classdef tUpstreamShare < matlab.unittest.TestCase
% rheome.flow.upstreamshare -- the share of tile-to-tile travel that runs up a gradient.
%
% Ground truth is planted: a plane wave along a chain of tiles whose gradient value is their position, carried
% by a 3-component analytic current with a different oscillation axis on every tile, on the outward side of the
% tile's normal (eig returns that axis with an arbitrary sign: the normal is what fixes it).
%
% Author: Diellor Basha, 2026

    methods (Static)
        function [Z, g, A, N] = chain(lagPerTile, opts)
            arguments
                lagPerTile (1,1) double                 % phase lag between neighbours, rad (> 0: travel up g)
                opts.nTile (1,1) double = 8
                opts.nT (1,1) double = 2000
                opts.Noise (1,1) double = 0
            end
            rng(3);  n = opts.nTile;  t = (0:opts.nT-1) / 200;
            g = (1:n)';  A = diag(ones(n-1, 1), 1);  A = A + A';
            N = randn(n, 3);  N = N ./ vecnorm(N, 2, 2);                      % tile mean outward normals
            u = N + 0.3 * randn(n, 3);  u = u ./ vecnorm(u, 2, 2);            % tilted axis, same side as N
            s = exp(1i * (2*pi*4*t - lagPerTile * g));                        % 4 Hz wave travelling along g
            Z = permute(u, [1 2 3]) .* reshape(s, n, 1, []);
            Z = Z + opts.Noise * (randn(size(Z)) + 1i * randn(size(Z)));
        end
    end

    methods (Test)
        function upwardTravelIsAllUpstream(tc)
            [Z, g, A, N] = tUpstreamShare.chain(0.6);
            out = rheome.flow.upstreamshare(Z, g, A, Normals=N);
            tc.verifyEqual(out.U, 1, AbsTol=1e-12);
            tc.verifyEqual(size(out.edges, 1), 7);
            tc.verifyGreaterThan(min(out.linearity), 1 - 1e-12);         % one axis per tile, planted
        end
        function downwardTravelIsNone(tc)
            [Z, g, A, N] = tUpstreamShare.chain(-0.6);
            tc.verifyEqual(rheome.flow.upstreamshare(Z, g, A, Normals=N).U, 0, AbsTol=1e-12);
        end
        function reversingTheGradientMirrorsU(tc)
            [Z, g, A, N] = tUpstreamShare.chain(0.6, Noise=0.5);
            a = rheome.flow.upstreamshare(Z, g, A, Normals=N).U;
            b = rheome.flow.upstreamshare(Z, -g, A, Normals=N).U;
            tc.verifyEqual(a + b, 1, AbsTol=1e-12);
        end
        function zeroLagReadsHalfUnderNoise(tc)
            % a standing (zero-lag) pattern -- what the point spread makes -- has no direction: U -> 0.5
            [Z, g, A, N] = tUpstreamShare.chain(0, Noise=0.3, nT=20000);
            tc.verifyEqual(rheome.flow.upstreamshare(Z, g, A, Normals=N).U, 0.5, AbsTol=0.02);
        end
        function aWrongAnchorInvertsThatTilesVotes(tc)
            % the polarity limit: one tile's anchor on the wrong side flips its two edges (2 of 7, equal weights)
            [Z, g, A, N] = tUpstreamShare.chain(0.6);
            N(4, :) = -N(4, :);
            tc.verifyEqual(rheome.flow.upstreamshare(Z, g, A, Normals=N).U, 5/7, AbsTol=1e-12);
        end
        function lagPastHalfACycleReadsBackwards(tc)
            % the documented aliasing limit: a lag of 1.2 pi is read as -0.8 pi
            [Z, g, A, N] = tUpstreamShare.chain(1.2*pi);
            tc.verifyEqual(rheome.flow.upstreamshare(Z, g, A, Normals=N).U, 0, AbsTol=1e-12);
        end
        function maskedTilesCastNoVote(tc)
            [Z, g, A, N] = tUpstreamShare.chain(0.6);
            m = true(8, 1);  m(4) = false;
            out = rheome.flow.upstreamshare(Z, g, A, Normals=N, Mask=m);
            tc.verifyEqual(size(out.edges, 1), 5);
            tc.verifyEqual(out.keptWeight, 5/7, AbsTol=1e-12);
            tc.verifyEqual(out.U, 1, AbsTol=1e-12);
        end
        function scalarFieldIsAccepted(tc)
            [~, g, A] = tUpstreamShare.chain(0.6);
            s = exp(1i * (2*pi*4*(0:1999)/200 - 0.6 * g));                  % any complex scalar per tile
            tc.verifyEqual(rheome.flow.upstreamshare(s, g, A).U, 1, AbsTol=1e-12);
        end
    end
end
