classdef tIngestBound < matlab.unittest.TestCase
% Stored extrema are BOUNDS, not rounded values.
%
% min, max, absMax and envMax are kept in single. Round-to-nearest can put a minimum one ulp
% above the samples it summarises, and an extremum that is not a bound breaks two things at
% once: the min/max envelope drawn from it excludes samples it covers, and a pruning descent
% that trusts "parent >= child" may drop a subtree it should keep. Rounding outward costs a
% relative 6e-8 and restores both guarantees.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200;  nT = 2000;  C = 2;  X; cfg; fb; bands; frame; g; T0
    end

    methods (TestClassSetup)
        function reduce(tc)
            rng(7);
            t = (0:tc.nT-1)' / tc.fs;
            % values around 1e-13, where single's ulp is ~1e-20: the rounding is invisible
            % in magnitude and decisive in direction.
            tc.X = 1e-13 * randn(tc.nT, tc.C) + 3e-13 * sin(2*pi*10*t);
            tc.cfg = rheome.ingest.config();
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            tc.T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
        end
    end

    methods (Test)

        function everyTilesPairContainsItsSamples(tc)
            % ⭐ The property the envelope rests on, at level 0 and at every rolled-up level.
            T = rheome.ingest.rollup(tc.T0, tc.g);
            F = tc.g.F;
            for L = [0 2 5]
                P = T{L+1};  n = F * 2^L;
                for k = 1:tc.g.K(L+1)
                    a = (k-1)*n + 1;  b = min(k*n, tc.nT);
                    x = tc.X(a:b, :);
                    tc.verifyLessThanOrEqual(double(P.min(k, :)), min(x, [], 1), sprintf('L%d k%d min', L, k));
                    tc.verifyGreaterThanOrEqual(double(P.max(k, :)), max(x, [], 1), sprintf('L%d k%d max', L, k));
                    tc.verifyGreaterThanOrEqual(double(P.absMax(k, :)), max(abs(x), [], 1), sprintf('L%d k%d absMax', L, k));
                end
            end
        end

        function theBoundIsTightToOneUlp(tc)
            % Outward rounding must not become slack: the bound is the nearest single on the
            % safe side, never further.
            F = tc.g.F;
            for k = 1:tc.g.K(1)
                x = tc.X((k-1)*F + (1:F), :);
                mn = min(x, [], 1);  mx = max(x, [], 1);
                tc.verifyLessThanOrEqual(mn - double(tc.T0.min(k, :)), eps(single(mn)));
                tc.verifyLessThanOrEqual(double(tc.T0.max(k, :)) - mx, eps(single(mx)));
            end
        end

        function roundingOutwardKeepsTheRollUpExact(tc)
            % up() and down() are monotone, so max over children of up(x) is up of the max:
            % a parent computed by roll-up still equals one computed directly.
            T = rheome.ingest.rollup(tc.T0, tc.g);
            for L = [1 3]
                gL = rheome.ingest.grid(tc.nT, tc.fs, rheome.ingest.config(FrameFloor = tc.cfg.FrameFloor * 2^L));
                D = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, gL, tc.frame);
                tc.verifyEqual(T{L+1}.min, D.min, sprintf('level %d min', L));
                tc.verifyEqual(T{L+1}.max, D.max, sprintf('level %d max', L));
                tc.verifyEqual(T{L+1}.absMax, D.absMax, sprintf('level %d absMax', L));
                tc.verifyEqual(T{L+1}.envMax, D.envMax, sprintf('level %d envMax', L));
            end
        end

        function theHelperRoundsTheWayItSays(tc)
            x = [1e-13, -1e-13, 0, 1, -1, Inf, -Inf, realmin];
            d = rheome.ingest.bound(x, 'down');  u = rheome.ingest.bound(x, 'up');
            tc.verifyClass(d, 'single');
            tc.verifyTrue(all(double(d) <= x));
            tc.verifyTrue(all(double(u) >= x));
            tc.verifyEqual(d([6 7]), single([Inf -Inf]));          % infinities untouched, not NaN
            tc.verifyEqual(u([6 7]), single([Inf -Inf]));
            tc.verifyError(@() rheome.ingest.bound(1, 'sideways'), 'ingest:bound:dir');
        end

    end
end

% Author: Diellor Basha, 2026
