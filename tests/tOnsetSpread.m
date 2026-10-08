classdef tOnsetSpread < matlab.unittest.TestCase
% rheome.detect.onsetspread: burst-onset slowness from half-rise times across tiles.
%
% What is asserted: a burst whose onset is delayed by distance/c reads slowness 1/c; a SIMULTANEOUS
% rise reads 0 however different the tiles' amplitudes are (the reason onsets are half-rise times and
% not threshold crossings); converging reads negative; a planar front reads its slowness vector; and
% too few participating tiles returns ok = false.
%
% Author: Diellor Basha, 2026

    properties
        G; S; fs = 300
    end

    methods (TestClassSetup)
        function tiles(tc)
            [V, F] = rheome.geom.icosphere(4);
            tc.S = struct('Vertices', V * 0.07, 'Faces', F);
            [L, M] = rheome.operators.laplace_beltrami(tc.S.Vertices, tc.S.Faces);
            T = rheome.geom.tree(tc.S, L=L, M=M, MaxDepth=6);
            tc.G = rheome.geom.tiles(T, tc.S, 6, Ruler=rheome.geom.edgegraph(tc.S));
        end
    end

    methods (Test)

        function spreadingReadsItsSlowness(tc)
            s = 2;                                              % 0.5 m/s
            [Y, tp] = tc.burst(@(d) s * d, @(d) ones(size(d)));
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, tp);
            tc.verifyTrue(R.ok);
            tc.verifyEqual(R.slowness, s, 'RelTol', 0.1);
            tc.verifyGreaterThan(R.r2, 0.9);
        end

        function aSimultaneousRiseReadsZeroWhateverTheAmplitudes(tc)
            [Y, tp] = tc.burst(@(d) 0 * d, @(d) 1 + 3 * rand(size(d)));
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, tp);
            tc.verifyLessThan(abs(R.slowness), 0.05);
            tc.verifyLessThan(R.spreadS, 0.01);
        end

        function aThresholdCrossingWouldNotHave(tc)
            % the trap the half-rise avoids: same simultaneous burst, fixed threshold -> fake spread
            [Y, tp] = tc.burst(@(d) 0 * d, @(d) 1 + 3 * d / max(d));
            thrT = arrayfun(@(k) find(Y(k,:) > 0.5, 1) / tc.fs, 1:size(Y,1));
            tc.verifyGreaterThan(range(thrT), 0.05);
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, tp);
            tc.verifyLessThan(abs(R.slowness), 0.05);
        end

        function convergingReadsNegative(tc)
            [Y, tp] = tc.burst(@(d) 0.2 - 2 * d, @(d) ones(size(d)));
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, tp, Pre=1.5);
            tc.verifyLessThan(R.slowness, -1);
        end

        function aPlanarFrontReadsItsSlownessVector(tc)
            % the front runs along the neighbourhood's own first tangent direction: a global axis can
            % point along the surface NORMAL at tile 1, where no front on the surface can lie
            c = tc.S.Vertices(tc.G.centre, :);  nb = tc.G.D(1, :) <= 0.05;
            [~, ~, V] = svd(c(nb,:) - mean(c(nb,:), 1), 'econ');  x = (c - c(1,:)) * V(:, 1);
            s = 3;  tp = 1.5;  nT = round(3 * tc.fs);  t = (0:nT-1) / tc.fs;
            Y = zeros(numel(tc.G.node_id), nT);
            for k = 1:size(Y, 1), Y(k,:) = tc.hann(t - (tp - 0.3) - s * x(k)); end
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, tp, Centres=c, Radius=0.05);
            tc.verifyEqual(R.planar, s, 'RelTol', 0.25);
            tc.verifyGreaterThan(R.planarR2, 0.8);
        end

        function tooFewTilesIsNotAFit(tc)
            Y = zeros(numel(tc.G.node_id), 900);  Y(1, 300:400) = 1;
            R = rheome.detect.onsetspread(Y, tc.G, tc.fs, 1, 350 / tc.fs);
            tc.verifyFalse(R.ok);
        end
    end

    methods
        % Tile envelopes: a 0.6 s Hann burst delayed by delayFn(distance from tile 1), scaled by ampFn.
        function [Y, tp] = burst(tc, delayFn, ampFn)
            d = tc.G.D(1, :)';  tp = 1.5;  nT = round(3 * tc.fs);  t = (0:nT-1) / tc.fs;
            dl = delayFn(d);  a = ampFn(d);  Y = zeros(numel(d), nT);
            for k = 1:numel(d), Y(k,:) = a(k) * tc.hann(t - (tp - 0.3) - dl(k)); end
        end
        function w = hann(~, u)
            w = zeros(size(u));  m = u >= 0 & u <= 0.6;  w(m) = sin(pi * u(m) / 0.6).^2;
        end
    end
end

% Author: Diellor Basha, 2026
