classdef tIngestPages < matlab.unittest.TestCase
% The page plan: one level per band from its support, the floor and the cap, contiguous
% jobs, and the byte split. A page is a tile at a level; this pins which level.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 4000
        g
        bands
        frame
    end

    methods (TestClassSetup)
        function makeBank(tc)
            cfg = rheome.ingest.config(Bank="morse");
            [~, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, cfg);
        end
    end

    methods (Test)

        function everyBandIsInExactlyOneJob(tc)
            plan = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1));
            all_ = sort(cell2mat(plan.bands'));
            tc.verifyEqual(all_, 1:height(tc.bands));
        end

        function theLevelIsTheSmallestFrameThatAffordsTheHalo(tc)
            cfg = rheome.ingest.config(Bank="morse", MinPageLength=0.25, PageOverhead=0.5);
            plan = rheome.ingest.pages(tc.g, tc.bands, cfg);
            for r = 1:height(plan)
                for b = plan.bands{r}
                    L = plan.level(r);
                    need = 2 * tc.bands.tSupport(b) / cfg.PageOverhead;
                    if L < tc.g.Lmax
                        tc.verifyGreaterThanOrEqual(tc.g.tExtent(L+1), need);
                        if L > 0, tc.verifyLessThan(tc.g.tExtent(L), need); end
                    end
                end
            end
        end

        function slowBandsRunOnTheWholeRecordWithNoHalo(tc)
            plan = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1));
            top = plan(plan.level == tc.g.Lmax, :);
            tc.verifyNotEmpty(top);
            tc.verifyTrue(all(top.haloSamples == 0));
            tc.verifyTrue(all(top.nPages == 1));
        end

        function theFloorIsHonoured(tc)
            plan = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=4));
            Lmin = ceil(log2(4 / 0.25));
            tc.verifyTrue(all(plan.level >= Lmin));
        end

        function haloIsTheSupportInSamplesAndSpanIsCorePlusTwoHalos(tc)
            plan = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1));
            for r = 1:height(plan)
                if plan.level(r) < tc.g.Lmax
                    tc.verifyEqual(plan.haloSamples(r), max(ceil(tc.bands.tSupport(plan.bands{r}) * tc.fs)));
                    tc.verifyEqual(plan.spanMax(r), min(tc.nT, plan.coreSamples(r) + 2*plan.haloSamples(r)));
                end
            end
        end

        function levelsAreNonDecreasingTowardsLowBands(tc)
            plan = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1));
            tc.verifyTrue(all(diff(plan.level) >= 0));
        end

        function aTightByteBudgetSplitsJobsAndAnImpossibleOneErrors(tc)
            loose = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1));
            tight = rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1, MaxBytes=6e6));
            tc.verifyGreaterThan(height(tight), height(loose));
            tc.verifyTrue(all(tight.bufferBytes <= 6e6));
            tc.verifyError(@() rheome.ingest.pages(tc.g, tc.bands, rheome.ingest.config(Bank="morse", MinPageLength=1, MaxBytes=1e3)), 'ingest:pages:bytes');
        end

    end
end
