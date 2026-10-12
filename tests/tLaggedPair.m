classdef tLaggedPair < matlab.unittest.TestCase
% rheome.forward.laggedpair on sensor data built from known generators: it finds the planted pair, its lag and
% amplitude ratio, keeps them in the surrogate, and handles one generator per hemisphere (Group).
%
% Author: Diellor Basha, 2026

    properties
        fs = 100
        t
        w
    end

    methods (TestMethodSetup)
        function signal(tc)
            tc.t = (-240:240) / tc.fs;  tc.w = -cos(2*pi*0.8*tc.t) .* exp(-tc.t.^2 / (2 * 0.5^2));
        end
    end

    methods (Test)
        function findsThePlantedPairLagAndRatio(tc)
            rng(1);  Gp = randn(60, 6);  front = logical([1 1 1 0 0 0]);
            Y = Gp(:, 2) * tc.w + 0.7 * Gp(:, 5) * tc.d(0.15) + 0.01 * randn(60, numel(tc.t));
            [Ys, f] = rheome.forward.laggedpair(Y, Gp, front, tc.fs);
            tc.verifyEqual([f.front f.back], [2 5]);
            tc.verifyEqual(f.lagS, 0.15, 'AbsTol', 0.011);
            tc.verifyEqual(f.ratio, 0.7, 'RelTol', 0.02);
            tc.verifyGreaterThan(f.r2, 0.99);
            [~, g] = rheome.forward.laggedpair(Ys, Gp, front, tc.fs);           % the null keeps what it was fitted to
            tc.verifyEqual([g.lagS g.ratio], [f.lagS f.ratio], 'AbsTol', 0.02);
        end

        function oneGeneratorAcrossHemispheres(tc)
            rng(2);  Gp = randn(60, 8);  front = logical([1 1 0 0 1 1 0 0]);  grp = [1 1 1 1 2 2 2 2];
            Y = (Gp(:, 1) + 0.8 * Gp(:, 6)) * tc.w + 0.5 * (Gp(:, 4) + Gp(:, 7)) * tc.d(-0.1) + 0.01 * randn(60, numel(tc.t));
            [~, f] = rheome.forward.laggedpair(Y, Gp, front, tc.fs, Group=grp);
            tc.verifyEqual(f.front, [1 6]);  tc.verifyEqual(f.back, [4 7]);
            tc.verifyEqual(f.lagS, -0.1, 'AbsTol', 0.011);                     % the back leads: negative lag
            tc.verifyEqual(f.ratio, norm([0.5 0.5]) / norm([1 0.8]), 'RelTol', 0.02);   % unit-norm patch currents
            tc.verifyGreaterThan(f.r2, 0.99);
        end

        function ascentFindsWhatEveryCombinationFinds(tc)
            rng(1);  Gp = randn(60, 6);  front = logical([1 1 1 0 0 0]);
            Y = Gp(:, 2) * tc.w + 0.7 * Gp(:, 5) * tc.d(0.15) + 0.01 * randn(60, numel(tc.t));
            [~, f] = rheome.forward.laggedpair(Y, Gp, front, tc.fs, MaxCombos=1);
            tc.verifyEqual([f.front f.back], [2 5]);
        end

        function keepWithZeroLagResidualKeepsTheData(tc)
            % Keep: the fitted generators' own time courses; zerolag: a residual whose two channels are a quarter
            % cycle apart comes back with no lag between them (in phase or in antiphase)
            rng(3);  Gp = randn(60, 6);  front = logical([1 1 1 0 0 0]);
            Y = Gp(:, 3) * tc.w + 0.7 * Gp(:, 4) * tc.d(0.2);
            Ys = rheome.forward.laggedpair(Y, Gp, front, tc.fs, Keep=true, Residual="zerolag");
            tc.verifyEqual(Ys, Y, 'AbsTol', 1e-9 * max(abs(Y(:))));
            n = 1000;  ph = 2*pi*0.05*(0:n-1);  R = [sin(ph); cos(ph); 1e-3 * [ph; -ph] / n];   % the generators take only rows 3-4
            rng(4);  Rz = rheome.forward.laggedpair(R, [0 0; 0 0; 1 0; 0 1], [true false], tc.fs, Residual="zerolag");
            c = corrcoef(Rz(1:2, :)');  tc.verifyEqual(abs(c(1, 2)), 1, 'AbsTol', 1e-6);
            tc.verifyEqual(norm(Rz(1:2, :), 'fro'), norm(R, 'fro'), 'RelTol', 0.5);   % same power, in expectation
        end

        function needsBothSidesInEveryGroup(tc)
            tc.verifyError(@() rheome.forward.laggedpair(randn(10, 50), randn(10, 2), [true true], tc.fs), ?MException);
        end
    end

    methods
        function y = d(tc, lagS)
            y = interp1(tc.t, tc.w, tc.t - lagS, 'linear', 0);
        end
    end
end

% Author: Diellor Basha, 2026
