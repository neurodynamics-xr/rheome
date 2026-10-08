classdef tIngestRollup < matlab.unittest.TestCase
% The store's central property: a level computed by merging equals the same level computed
% directly from the samples. Sums to round-off, extrema and counts exactly. A record with
% an odd frame count at several levels exercises the neutral element.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 4050          % K0 = 81: odd at level 0, then 41, 21, 11, 6, 3, 2, 1
        X
        cfg
        fb
        bands
        frame
    end

    methods (TestClassSetup)
        function makeSignal(tc)
            rng(11);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, 2) + 4e-13*sin(2*pi*7*t);
            tc.cfg = rheome.ingest.config(Bank="morse");
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
        end
    end

    methods (Test)

        function everyLevelIsPresentAndTheTopHasOneRow(tc)
            g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            T = rheome.ingest.rollup(rheome.ingest.reduce(tc.X, tc.fb, tc.bands, g, tc.frame), g);
            tc.verifyNumElements(T, g.Lmax + 1);
            tc.verifyEqual(g.K, [81 41 21 11 6 3 2 1]);
            for L = 0:g.Lmax
                tc.verifyEqual(T{L+1}.level, L);
                tc.verifyEqual(size(T{L+1}.sumX, 1), g.K(L+1));
            end
            tc.verifyEqual(T{end}.n, uint32(tc.nT));
        end

        function rolledUpLevelsEqualDirectComputation(tc)
            g0 = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            T  = rheome.ingest.rollup(rheome.ingest.reduce(tc.X, tc.fb, tc.bands, g0, tc.frame), g0);
            for L = [1 2 4 7]
                cfgL = rheome.ingest.config(Bank="morse", FrameFloor = tc.cfg.FrameFloor * 2^L);
                gL = rheome.ingest.grid(tc.nT, tc.fs, cfgL);
                D  = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, gL, tc.frame);
                R  = T{L+1};
                tc.verifyEqual(R.n,      D.n);
                tc.verifyEqual(R.nCoi,   D.nCoi);
                tc.verifyEqual(R.sumX,   D.sumX,   'RelTol', 1e-12, 'AbsTol', 1e-25);
                tc.verifyEqual(R.sumX2,  D.sumX2,  'RelTol', 1e-12);
                tc.verifyEqual(R.energy, D.energy, 'RelTol', 1e-12);
                tc.verifyEqual(R.absMax, D.absMax);
                tc.verifyEqual(R.envMax, D.envMax);
                tc.verifyEqual(R.min,    D.min);
                tc.verifyEqual(R.max,    D.max);
            end
        end

        function anAbsentChildIsNeutral(tc)
            % Level 0 has 81 rows; level 1 row 41 has one real child. Its statistics must
            % equal that child's, not be corrupted by the padding.
            g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, g, tc.frame);
            T  = rheome.ingest.rollup(T0, g);
            tc.verifyEqual(T{2}.n(41),        T0.n(81));
            tc.verifyEqual(T{2}.sumX(41, :),  T0.sumX(81, :));
            tc.verifyEqual(T{2}.min(41, :),   T0.min(81, :));
            tc.verifyEqual(T{2}.max(41, :),   T0.max(81, :));
            tc.verifyEqual(T{2}.energy(41, :, :), T0.energy(81, :, :));
        end

        function aMaskRollsUpLikeAnyOtherCount(tc)
            mask = true(tc.nT, 1);  mask(101:180) = false;      % 80 samples across frames 3-4
            g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            T = rheome.ingest.rollup(rheome.ingest.reduce(tc.X, tc.fb, tc.bands, g, tc.frame, Mask=mask), g);
            tc.verifyEqual(T{end}.n, uint32(tc.nT - 80));
            tc.verifyEqual(T{2}.n(2), uint32(100 - 80));         % level-1 frame 2 = frames 3,4
        end

    end
end
