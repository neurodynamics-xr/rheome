classdef tIngestGrid < matlab.unittest.TestCase
% The dyadic grid is arithmetic that every stored row depends on. An off-by-one here
% misplaces every tile silently, so the frame counts, the top level and the parent/child
% pairing are pinned against hand-computed values.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function frameCountsPerLevelAreCeilOfHalving(tc)
            cfg = rheome.ingest.config();                      % 0.25 s
            g = rheome.ingest.grid(1000, 8, cfg);              % F = 2 samples, K0 = 500
            tc.verifyEqual(g.F, 2);
            tc.verifyEqual(g.K0, 500);
            tc.verifyEqual(g.K, ceil(500 ./ 2.^(0:g.Lmax)));
        end

        function topLevelIsTheFirstWithOneFrame(tc)
            g = rheome.ingest.grid(1000, 8, rheome.ingest.config());
            tc.verifyEqual(g.K(end), 1);
            tc.verifyGreaterThan(g.K(end-1), 1);
            tc.verifyEqual(g.Lmax, 9);                  % 2^9 = 512 >= 500
        end

        function everyParentHasExactlyTwoChildren(tc)
            % K(L) = ceil(K(L-1)/2): the last parent may have one real child, never zero
            % and never three.
            g = rheome.ingest.grid(1000, 8, rheome.ingest.config());
            for L = 1:g.Lmax
                tc.verifyEqual(g.K(L+1), ceil(g.K(L)/2));
            end
        end

        function aSingleFrameRecordHasLevelZeroOnly(tc)
            g = rheome.ingest.grid(2, 8, rheome.ingest.config());
            tc.verifyEqual(g.Lmax, 0);
            tc.verifyEqual(g.K, 1);
        end

        function partialLastFrameKeepsItsNominalCentre(tc)
            % 7 samples at 4 Hz with 0.25 s frames: F = 1, K0 = 7. Level 1 has 4 frames of
            % 0.5 s; the fourth holds one sample but is centred at 1.75 s like a full one.
            g = rheome.ingest.grid(7, 4, rheome.ingest.config());
            tc.verifyEqual(g.K(2), 4);
            tc.verifyEqual(g.tExtent(2), 0.5);
            tc.verifyEqual(g.tCenter{2}, [0.25 0.75 1.25 1.75], 'AbsTol', 1e-12);
        end

        function extentsDoubleAndCentresTileTheAxis(tc)
            g = rheome.ingest.grid(1000, 8, rheome.ingest.config());
            tc.verifyEqual(g.tExtent, 0.25 * 2.^(0:g.Lmax), 'RelTol', 1e-12);
            tc.verifyEqual(diff(g.tCenter{1}), 0.25 * ones(1, g.K0-1), 'AbsTol', 1e-12);
            tc.verifyEqual(g.tCenter{1}(1), 0.125, 'AbsTol', 1e-12);
        end

        function nonIntegerFrameFloorTimesFsErrors(tc)
            tc.verifyError(@() rheome.ingest.grid(1000, 10, rheome.ingest.config()), 'ingest:grid:frameFloor');
        end

        function configFieldOrderIsCanonical(tc)
            cfg = rheome.ingest.config(MaxBytes=1e9, FrameFloor=0.5);
            tc.verifyEqual(fieldnames(cfg)', ...
                {'FrameFloor','Wavelet','VoicesPerOctave','FrequencyLimits','Precision','MaxBytes','Paged','PageOverhead','MinPageLength','MaxSupport','Bank','Anchor','Oversample','Space','SpaceVoices','ChannelMinTile','ChannelEnvelope','Moments','Peaks','Couple','PhaseBins','SpaceDiagonal','hash'});
        end

        function hashDependsOnParametersNotOnOptionOrder(tc)
            a = rheome.ingest.config();
            b = rheome.ingest.config(Wavelet="morse", FrameFloor=0.25);     % defaults, other order
            c = rheome.ingest.config(VoicesPerOctave=6);
            d = rheome.ingest.config(FrequencyLimits=[1 60]);
            e = rheome.ingest.config(MaxBytes=1e9);                        % a budget, not in the hash
            tc.verifyEqual(a.hash, b.hash);
            tc.verifyEqual(a.hash, e.hash);
            tc.verifyNotEqual(a.hash, c.hash);
            tc.verifyNotEqual(a.hash, d.hash);
            tc.verifyLength(a.hash, 8);
        end

        function badFrequencyLimitsError(tc)
            tc.verifyError(@() rheome.ingest.config(FrequencyLimits=[60 1]), 'ingest:config:limits');
            tc.verifyError(@() rheome.ingest.config(FrequencyLimits=5), 'ingest:config:limits');
        end

    end
end
