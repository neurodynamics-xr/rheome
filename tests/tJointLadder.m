classdef tJointLadder < matlab.unittest.TestCase
% The joint (space x time) bookkeeping ladder and its picker: rheome.geom.jointladder, rheome.geom.jointcell.
%
% A row is an observation window -- tile size, octave window, frame rate -- and its speed columns
% are the rate written in metres per second. What is asserted: the speeds are D*rate, D*rate/MinDwell
% and D/window exactly; alpha's octave gets the atlas's 2 s window and a carrier near 30 frames per
% cycle; and the picker takes the coarsest tile that holds the feature, the coarsest rate that
% follows it, the carrier for a phase, and refuses what no row can hold.
%
% Author: Diellor Basha, 2026

    properties
        L; T
    end

    methods (TestClassSetup)
        function ladder(tc)
            [V, F] = rheome.geom.icosphere(4);
            S = struct('Vertices', V * 0.07, 'Faces', F);
            [Lp, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);
            tc.T = rheome.geom.tree(S, L=Lp, M=M, MaxDepth=8);          % tiles down to 17 mm
            tc.L = rheome.geom.jointladder(tc.T, Fs=600);
        end
    end

    methods (Test)

        function speedsAreTheRateInMetresPerSecond(tc)
            L = tc.L;  D = L.tileMM / 1e3;
            tc.verifyEqual(L.speedMaxMS, D .* L.frameHz, 'RelTol', 1e-12);
            tc.verifyEqual(L.speedTrackMS, D .* L.frameHz / 2.5, 'RelTol', 1e-12);
            tc.verifyEqual(L.speedMinMS, D ./ L.windowS, 'RelTol', 1e-12);
            tc.verifyEqual(L.frameHz, 600 ./ 2.^L.level, 'RelTol', 1e-12);
        end

        function alphaGetsTheAtlasWindowAndAThirtyPerCycleCarrier(tc)
            A = tc.L(tc.L.fLo == 8, :);
            tc.verifyEqual(unique(A.windowS), 2);                 % k = 15 at 4 voices: 15/8 -> 2 s
            c = A(A.isCarrier, :);
            tc.verifyEqual(unique(c.frameHz), 300);
            tc.verifyEqual(unique(A.framesPerWindow(A.level == max(A.level))), ...
                2 * 600 / 2^max(A.level), 'RelTol', 1e-12);
            tc.verifyGreaterThanOrEqual(min(A.framesPerWindow), 1);
        end

        function eigenBandIsOneOctaveOfWavenumber(tc)
            tc.verifyEqual(tc.L.lambdaHi ./ tc.L.lambdaLo, 4 * ones(height(tc.L), 1), 'RelTol', 1e-12);
            ts = tc.L.Properties.UserData.TileSigma;
            tc.verifyEqual(tc.L.sigmaMinMM, tc.L.tileMM / ts, 'RelTol', 1e-12);
            tc.verifyEqual(tc.L.lambdaHi, 2 ./ (tc.L.sigmaMinMM / 1e3).^2, 'RelTol', 1e-12);
        end

        function pickerTakesTheCoarsestTileThatHoldsTheFeature(tc)
            sg = 30;  [row, why] = rheome.geom.jointcell(tc.L, SigmaMM=sg, Hz=10, SpeedMS=0.1);
            tc.verifyNotEmpty(row, why);
            ts = tc.L.Properties.UserData.TileSigma;
            tc.verifyLessThanOrEqual(row.tileMM, ts * sg);
            coarser = unique(tc.L.tileMM(tc.L.tileMM > row.tileMM));
            tc.verifyTrue(all(coarser > ts * sg), 'a coarser tile would also have held it');
        end

        function tileSigmaMovesThePickOneDepthPerSqrtTwo(tc)
            L2 = rheome.geom.jointladder(tc.T, Fs=600, TileSigma=sqrt(2));
            L1 = rheome.geom.jointladder(tc.T, Fs=600, TileSigma=1);
            r2 = rheome.geom.jointcell(L2, SigmaMM=30, Hz=10, SpeedMS=0.1);
            r1 = rheome.geom.jointcell(L1, SigmaMM=30, Hz=10, SpeedMS=0.1);
            tc.verifyEqual(r1.depth, r2.depth + 1);
        end

        function envelopeRollsUpToTheCoarsestRateThatFollows(tc)
            row = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=0.1);
            tc.verifyGreaterThanOrEqual(row.speedTrackMS, 0.1);
            same = tc.L(tc.L.depth == row.depth & tc.L.fLo == row.fLo & tc.L.level == row.level + 1, :);
            if ~isempty(same), tc.verifyLessThan(same.speedTrackMS, 0.1); end
            fast = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=1);
            tc.verifyGreaterThan(fast.frameHz, row.frameHz);
        end

        function phaseIsPinnedToTheCarrier(tc)
            row = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=0.01, Measure="phase");
            tc.verifyTrue(row.isCarrier);
            tc.verifyEqual(row.frameHz, 300);
        end

        function detectabilityIsEvidenceOverTheEvent(tc)
            % a 120 mm path: slow is long enough at 10 dB, fast is not; 20 dB rescues the fast one
            slow = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=0.1, SnrDB=10, PathMM=120);
            fast = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=1, SnrDB=10, PathMM=120);
            loud = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=1, SnrDB=20, PathMM=120);
            tc.verifyEqual(slow.eventS, 1.2, 'AbsTol', 1e-12);
            tc.verifyEqual(slow.nTau, 1.2 * 8, 'AbsTol', 1e-12);             % tau = 1/(16-8)
            tc.verifyEqual(slow.evidence, 10 * 9.6, 'AbsTol', 1e-9);
            tc.verifyTrue(slow.detectable);  tc.verifyFalse(fast.detectable);  tc.verifyTrue(loud.detectable);
            still = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=0, SnrDB=10);
            tc.verifyEqual(still.eventS, still.windowS);
        end

        function refusalsSayWhy(tc)
            [row, why] = rheome.geom.jointcell(tc.L, SigmaMM=1, Hz=10);
            tc.verifyEmpty(row);  tc.verifySubstring(char(why), 'No depth is fine enough');
            [row, why] = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=10, SpeedMS=1e4);
            tc.verifyEmpty(row);  tc.verifySubstring(char(why), 'Too fast');
            [row, why] = rheome.geom.jointcell(tc.L, SigmaMM=30, Hz=500);
            tc.verifyEmpty(row);  tc.verifySubstring(char(why), 'outside');
        end
    end
end

% Author: Diellor Basha, 2026
