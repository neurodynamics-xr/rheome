classdef tIngestReduce < matlab.unittest.TestCase
% Level 0 is the only place samples are read; everything above it is a merge. So what is
% pinned here is (a) equality with an independent oracle on full frames, (b) the mask's
% arithmetic, (c) the band energies partitioning the coefficient energy exactly, and (d)
% the frame-bound identity for a tone. Synthetic, 20 s at 200 Hz.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 4000
        C  = 3
        X
        cfg
        fb
        bands
        frame
        g
    end

    methods (TestClassSetup)
        function makeSignal(tc)
            rng(7);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, tc.C) + 5e-13*sin(2*pi*10*t) .* (1 + 0.5*sin(2*pi*0.2*t));
            tc.cfg = rheome.ingest.config(Bank="morse");
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);        % F = 50, K0 = 80
        end
    end

    methods (Test)

        function allPassStatisticsEqualTheMathWorksExtractorOnFullFrames(tc)
            % signalTimeFeatureExtractor is the oracle: on full, unmasked frames the
            % hand-rolled reduce must reproduce it. Mean*F and RMS^2*F are sums up to
            % round-off; PeakValue is max|x| up to the OUTWARD rounding the store applies to
            % every extremum (rheome.ingest.bound), which is at most one single ulp above.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            F = tc.g.F;
            ex = signalTimeFeatureExtractor("FrameSize", F, "FrameRate", F, "SampleRate", tc.fs, ...
                    "IncompleteFrameRule", "drop", "FeatureFormat", "matrix", ...
                    Mean=true, RMS=true, PeakValue=true);
            M = extract(ex, tc.X);                              % [K0 x 3 x C]
            tc.verifyEqual(T0.sumX,  squeeze(M(:,1,:)) * F,    'RelTol', 1e-12, 'AbsTol', 1e-25);
            tc.verifyEqual(T0.sumX2, squeeze(M(:,2,:)).^2 * F, 'RelTol', 1e-12);
            pk = squeeze(M(:,3,:));
            tc.verifyGreaterThanOrEqual(double(T0.absMax), pk);                  % a bound
            tc.verifyLessThanOrEqual(double(T0.absMax) - pk, eps(single(pk)));   % tight
            tc.verifyEqual(T0.n, uint32(F * ones(tc.g.K0, 1)));
        end

        function minAndMaxBoundEveryFrameToOneUlp(tc)
            % ⚠ THE STORED PAIR IS A BOUND, NOT A ROUNDED VALUE. single() to nearest can put
            % a minimum one ulp above the samples it summarises; the store rounds outward
            % (rheome.ingest.bound) so the pair always contains the frame, which is what the min/max
            % envelope and the pruning descent both rest on. Tight to one ulp either way.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            Xb = reshape(tc.X, tc.g.F, tc.g.K0, tc.C);
            mn = squeeze(min(Xb, [], 1));  mx = squeeze(max(Xb, [], 1));
            tc.verifyLessThanOrEqual(double(T0.min), mn);
            tc.verifyGreaterThanOrEqual(double(T0.max), mx);
            tc.verifyLessThanOrEqual(mn - double(T0.min), eps(single(mn)));
            tc.verifyLessThanOrEqual(double(T0.max) - mx, eps(single(mx)));
        end

        function aMaskedRunIsCountedOutAndOtherFramesAreUntouched(tc)
            mask = true(tc.nT, 1);  mask(1001:1025) = false;     % half of frame 21
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, Mask=mask);
            U0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            tc.verifyEqual(T0.n(21), uint32(25));
            tc.verifyEqual(T0.sumX(21, :), sum(tc.X(1026:1050, :), 1), 'RelTol', 1e-12);
            mn = min(tc.X(1026:1050, :), [], 1);
            tc.verifyLessThanOrEqual(double(T0.min(21, :)), mn);                 % outward-rounded bound
            tc.verifyLessThanOrEqual(mn - double(T0.min(21, :)), eps(single(mn)));
            others = [1:20, 22:tc.g.K0];
            tc.verifyEqual(T0.sumX(others, :),  U0.sumX(others, :));
            tc.verifyEqual(T0.energy(others, :, :), U0.energy(others, :, :));
            tc.verifyEqual(T0.nCoi(others, :), U0.nCoi(others, :));
        end

        function aFullyMaskedFrameHoldsTheNeutralElements(tc)
            mask = true(tc.nT, 1);  mask(1001:1050) = false;     % all of frame 21
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, Mask=mask);
            tc.verifyEqual(T0.n(21), uint32(0));
            tc.verifyEqual(T0.sumX(21, :), zeros(1, tc.C));
            tc.verifyEqual(T0.absMax(21, :), zeros(1, tc.C, 'single'));
            tc.verifyEqual(T0.min(21, :), single(+Inf(1, tc.C)));
            tc.verifyEqual(T0.max(21, :), single(-Inf(1, tc.C)));
            tc.verifyEqual(squeeze(T0.energy(21, :, :)), zeros(tc.C, height(tc.bands)));
        end

        function aPartialLastFrameKeepsItsRealCount(tc)
            nT = tc.nT + 10;
            X = [tc.X; tc.X(1:10, :)];
            [fb, bands, frame] = rheome.ingest.bank(nT, tc.fs, tc.cfg);
            g = rheome.ingest.grid(nT, tc.fs, tc.cfg);
            T0 = rheome.ingest.reduce(X, fb, bands, g, frame);
            tc.verifyEqual(g.K0, 81);
            tc.verifyEqual(T0.n(end), uint32(10));
            tc.verifyEqual(T0.sumX(end, :), sum(X(4001:4010, :), 1), 'RelTol', 1e-12);
        end

        function coneCountIsTheSupportConeFromTheRecordEdges(tc)
            % nCoi is defined analytically (rheome.ingest.cone) with the SUPPORT constant, wider
            % than MATLAB's coi on purpose: it must cover every tile where the boundary
            % extension shows (tIngestPaged pins that). MATLAB's cone is contained in it.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            fc = centerFrequencies(tc.fb);
            lowest = cellfun(@(s) fc(s(end)), tc.bands.scales);
            inside = rheome.ingest.cone(tc.nT, tc.fs, tc.frame, lowest);
            for b = 1:height(tc.bands)
                tc.verifyEqual(T0.nCoi(:, b), uint32(sum(reshape(inside(b, :), tc.g.F, tc.g.K0), 1))');
            end
            [~, ~, coi] = wt(tc.fb, single(tc.X(:, 1)));
            matlabCone = lowest(:) < coi(:)';
            tc.verifyTrue(all(inside(matlabCone)));                 % ours contains MATLAB's
            tc.verifyGreaterThan(tc.frame.support.timeTimesFc, 3 * tc.frame.coiConst);
            % a band's cone reaches one support in from each edge
            b = 5;  d = (0:tc.nT-1)/tc.fs;  d = min(d, (tc.nT-1)/tc.fs - d);
            tc.verifyEqual(inside(b, :), d < tc.frame.support.timeTimesFc / lowest(b));
            % the slowest band is entirely inside; the top band mostly outside
            tc.verifyEqual(double(sum(T0.nCoi(:, end))), tc.nT);
            tc.verifyLessThan(double(sum(T0.nCoi(:, 1))), tc.nT / 10);
        end

        function bandEnergiesPartitionTheCoefficientEnergyExactly(tc)
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            for c = 1:tc.C
                W = wt(tc.fb, single(tc.X(:, c)));
                tot = sum(double(abs(W)).^2, 1);
                ref = sum(reshape(tot, tc.g.F, tc.g.K0), 1)';
                tc.verifyEqual(sum(T0.energy(:, c, :), 3), ref, 'RelTol', 1e-10);
            end
        end

        function envelopeMaxBoundsEverySingleScale(tc)
            % The pruning bound of the design: no scale in a band exceeds the band's
            % envMax on any frame. Exact by construction; pinned so it stays so.
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            W = abs(wt(tc.fb, single(tc.X(:, 2))));
            for b = 1:height(tc.bands)
                for m = tc.bands.scales{b}
                    per = max(reshape(W(m, :), tc.g.F, tc.g.K0), [], 1)';
                    tc.verifyTrue(all(per <= T0.envMax(:, 2, b)));
                end
            end
        end

        function aToneObeysTheFrameBoundIdentity(tc)
            % For a real tone at f0 inside the bank, sum_j energy_j / sumX2 is S(f0)/2 on
            % interior frames, and S/2 lies in [A, B] on the interior (frame.A/B are the
            % bounds of S/2 for exactly this reason). The tone is 8 Hz: two whole cycles
            % per 0.25 s frame, so sumX2 is F/2 on every frame and nothing wobbles with
            % phase. At a non-integer number of cycles per frame the per-frame ratio moves
            % by ~0.7 % -- the frame's doing, not the bank's.
            f0 = 8;
            t = (0:tc.nT-1)' / tc.fs;
            x = sin(2*pi*f0*t);
            [fb, bands, frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);   % own bank: one precision per bank
            T0 = rheome.ingest.reduce(x, fb, bands, tc.g, frame, Precision="double");
            b0 = find(bands.fLo <= f0 & f0 < bands.fHi);
            edge = 2 * ceil(bands.tSupport(b0) / tc.g.tExtent(1)) + 1;   % two supports in from the record edge
            k = (edge+1):(tc.g.K0-edge);
            tc.verifyEqual(T0.sumX2(k, 1), tc.g.F/2 * ones(numel(k), 1), 'RelTol', 1e-12);
            r = sum(T0.energy(k, 1, :), 3) ./ T0.sumX2(k, 1);
            tc.verifyGreaterThan(min(r), frame.A * (1 - 1e-3));
            tc.verifyLessThan(max(r), frame.B * (1 + 1e-3));
        end

        function theByteGuardFailsWithTheNumber(tc)
            tc.verifyError(@() rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, MaxBytes=1), 'ingest:reduce:bytes');
        end

        function lengthAndMaskMismatchesError(tc)
            tc.verifyError(@() rheome.ingest.reduce(tc.X(1:100, :), tc.fb, tc.bands, tc.g, tc.frame), 'ingest:reduce:length');
            tc.verifyError(@() rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame, Mask=true(5, 1)), 'ingest:reduce:mask');
        end

    end
end
