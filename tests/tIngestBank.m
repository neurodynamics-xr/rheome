classdef tIngestBank < matlab.unittest.TestCase
% The MORSE band table (Bank="morse"; the frame bank has tests/tIngestFrame.m).
% The band table is the frequency axis of every stored row. What is pinned: every scale
% belongs to exactly one band (so band energies partition the coefficient energy), full
% bands are exact contiguous octaves, measured extents contain nominal ones, and the bank
% is tight on its interior. Small record so the suite stays fast.
%
% Author: Diellor Basha, 2026

    properties (Constant)
        nT = 4000
        fs = 200
    end

    methods (Test)

        function everyScaleIsInExactlyOneBand(tc)
            [fb, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            all_ = sort(cell2mat(bands.scales'));
            tc.verifyEqual(all_, 1:numel(centerFrequencies(fb)));
        end

        function fullBandsHaveVoicesPerOctaveScalesAndSpanOneOctave(tc)
            cfg = rheome.ingest.config(Bank="morse");
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, cfg);
            full = ~bands.partial;
            tc.verifyTrue(any(full));
            tc.verifyTrue(all(cellfun(@numel, bands.scales(full)) == cfg.VoicesPerOctave));
            tc.verifyEqual(bands.fExtent(full), ones(nnz(full), 1), 'AbsTol', 1e-9);
        end

        function onlyTheBottomBandCanBePartialAndItsExtentSaysSo(tc)
            cfg = rheome.ingest.config(Bank="morse");
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, cfg);
            tc.verifyFalse(any(bands.partial(1:end-1)));
            if bands.partial(end)
                n = numel(bands.scales{end});
                tc.verifyEqual(bands.fExtent(end), n / cfg.VoicesPerOctave, 'AbsTol', 1e-9);
            end
        end

        function bandsAreContiguousInLogFrequency(tc)
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            tc.verifyEqual(bands.fLo(1:end-1), bands.fHi(2:end), 'RelTol', 1e-9);
        end

        function centreIsTheGeometricMean(tc)
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            tc.verifyEqual(bands.fCenter, sqrt(bands.fLo .* bands.fHi), 'RelTol', 1e-12);
        end

        function measuredExtentsContainNominalOnes(tc)
            % The wavelets respond beyond the octave they are assigned to. The table must
            % say so rather than pretend the octave edge is a wall.
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            tc.verifyTrue(all(bands.fLoMeasured < bands.fLo));
            tc.verifyTrue(all(bands.fHiMeasured > bands.fHi));
            tc.verifyTrue(all(bands.tSupport > 0));
        end

        function timeSupportGrowsTowardsLowBandsUntilTheRecordTruncatesIt(tc)
            % waveletsupport TRUNCATES at the record boundary (paged-io notes), so the
            % slowest bands of a short record share one saturated support. Below that the
            % support must grow strictly towards low bands; overall it never shrinks.
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            s = bands.tSupport;                                   % top band first
            tc.verifyTrue(all(diff(s) >= 0));
            free = s < 0.25 * tc.nT / tc.fs;                      % well inside the record
            tc.verifyGreaterThan(nnz(free), 3);
            tc.verifyTrue(all(diff(s(free)) > 0));
        end

        function theBankIsTightOnItsInterior(tc)
            [~, ~, frame] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            tc.verifyLessThan(frame.B / frame.A, 1.01);
            tc.verifyGreaterThan(frame.A, 0);
        end

        function theBankIsCwtfilterbankCappedBelowNyquist(tc)
            % Same scales as a plain cwtfilterbank up to the cap; the cap is where the
            % top filter's gain at Nyquist is 1e-6, about 0.7 of the default top.
            cfg = rheome.ingest.config(Bank="morse");
            [fb, ~, frame] = rheome.ingest.bank(tc.nT, tc.fs, cfg);
            ref = cwtfilterbank('SignalLength', tc.nT, 'SamplingFrequency', tc.fs, ...
                                'Wavelet', 'morse', 'VoicesPerOctave', 10);
            fcr = centerFrequencies(ref);  fcb = centerFrequencies(fb);
            [~, fmaxD] = cwtfreqbounds(tc.nT, tc.fs, 'Wavelet', 'morse');
            tc.verifyLessThan(frame.fmax, 0.75 * fmaxD);
            tc.verifyGreaterThan(frame.fmax, 0.6 * fmaxD);
            [H, f] = freqz(fb);
            tc.verifyLessThan(H(1, end) / max(H(1, :)), 1e-5);
            % the analytic grid: the bank's scales are fmax * 2^(-k/10)
            k = 0:numel(fcb)-1;
            tc.verifyEqual(fcb(:)', frame.fmax * 2.^(-k/10), 'RelTol', 1e-12);
            tc.verifyEqual(numel(fcr) - numel(fcb), 5);       % half an octave fewer, top only
        end

        function maxSupportRaisesTheFloorAndKeepsTheScales(tc)
            % MaxSupport = 2 s -> floor 6.46/2 = 3.2 Hz. Fewer scales, same top, same
            % scales where both banks have them, and every carried support <= 2 s.
            [~, bAll, fAll] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            [~, bCut, fCut] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse", MaxSupport=2));
            tc.verifyEqual(fCut.fmax, fAll.fmax, 'RelTol', 1e-12);
            tc.verifyLessThan(fCut.nF, fAll.nF);
            tc.verifyEqual(fCut.fc, fAll.fc(1:fCut.nF), 'RelTol', 1e-12);
            tc.verifyTrue(all(bCut.tSupport <= 2 * 1.02));
            tc.verifyGreaterThan(fCut.fmin, 3);
            tc.verifyLessThan(fCut.fmin, 3.5);
        end

        function naturalLevelIsTheFirstTileAtLeastTheSupport(tc)
            [~, bands] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            g = rheome.ingest.grid(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            for b = 1:height(bands)
                L = bands.naturalLevel(b);
                tc.verifyGreaterThanOrEqual(g.tExtent(L+1), bands.tSupport(b) * (1 - 1e-12));
                if L > 0, tc.verifyLessThan(g.tExtent(L), bands.tSupport(b)); end
            end
            tc.verifyEqual(bands.naturalLevel(1), 0);                 % top band: below the floor
            tc.verifyTrue(all(diff(bands.naturalLevel) >= 0));
        end

        function explicitLimitsAboveTheCapWarn(tc)
            [~, fmaxD] = cwtfreqbounds(tc.nT, tc.fs, 'Wavelet', 'morse');
            tc.verifyWarning(@() rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse", FrequencyLimits=[2 fmaxD])), 'ingest:bank:nyquist');
        end

        function explicitLimitsAreHonoured(tc)
            cfg = rheome.ingest.config(Bank="morse", FrequencyLimits=[2 40]);
            [fb, bands] = rheome.ingest.bank(tc.nT, tc.fs, cfg);
            fc = centerFrequencies(fb);
            tc.verifyGreaterThanOrEqual(min(fc), 2 * (1 - 1e-9));
            tc.verifyLessThanOrEqual(max(fc), 40 * (1 + 1e-9));
            tc.verifyEqual(bands.fHi(1), 40 * 2^(1/20), 'RelTol', 1e-9);
        end

        function fewerVoicesGiveFewerScalesSameOctaves(tc)
            [~, b10] = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse"));
            [~, b6]  = rheome.ingest.bank(tc.nT, tc.fs, rheome.ingest.config(Bank="morse", VoicesPerOctave=6));
            tc.verifyEqual(b6.fHi(1), b10.fHi(1), 'RelTol', 0.15);   % same top, own half-voice
            tc.verifyEqual(nnz(~b6.partial), nnz(~b10.partial));
        end

    end
end
