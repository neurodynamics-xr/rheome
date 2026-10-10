classdef tWaveletTile < matlab.unittest.TestCase
% rheome.selection.wavelettile: wavelets as a tiling. The level structure is identical to the hard dyadic
% tiling; only the within-level density differs, and it differs by one constant.
%
% Author: Diellor Basha, 2026
    properties
        T; oc
    end
    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                C = rheome.select.catalog();
                h = find(C.recording_id == string(rheomeTestSubject()) & C.bank == "frame" & ...
                         C.store == "default", 1);
                t.T = rheome.selection.wavelettile(rheome.select.open(char(C.file(h))));
            catch
                t.assumeFail('no default frame store for test subject');
            end
            % the true octave bands: the partial ones at the spectrum's edges obey none of this
            t.oc = t.T.supportPerTile > 0.6 & t.T.supportPerTile < 0.7;
        end
    end
    methods (Test)
        function oneSupportIsExactlyTwoOfTheNextLevel(t)
            % ⭐ the claim that makes this a dyadic tiling and not just a bank
            % ⚠ only between pairs where BOTH bands are full octaves; band 2 is 0.75 octaves, so its
            % neighbour's ratio is 0.595 and that is correct rather than a failure.
            r = t.T.childrenPerSupport;
            pair = t.oc & [false; t.oc(1:end-1)] & isfinite(r);
            t.verifyGreaterThan(sum(pair), 4, 'too few octave pairs to check the nesting');
            t.verifyEqual(r(pair), repmat(2, sum(pair), 1), 'RelTol', 1e-9);
        end
        function supportIsScaleInvariantInCycles(t)
            t.verifyEqual(t.T.supportCycles, repmat(t.T.supportCycles(1), height(t.T), 1), ...
                'RelTol', 1e-12, 'the support must be a fixed number of cycles');
        end
        function theSelectionStrideIsTheSupport(t)
            t.verifyEqual(t.T.strideSelect, t.T.supportSec, 'RelTol', 1e-12);
        end
        function theTileIsOneAndAHalfSupportsLong(t)
            % ⭐ the exchange rate with the hard tiling: the tile is the support rounded UP to the next
            % power of two, so one wavelet per tile covers two thirds of it
            t.verifyEqual(t.T.positionsPerTileSelect(t.oc), ...
                repmat(1.5033, sum(t.oc), 1), 'RelTol', 1e-3);
            t.verifyEqual(t.T.supportPerTile(t.oc), repmat(0.6652, sum(t.oc), 1), 'RelTol', 1e-3);
        end
        function theFrameRedundancyIsConstantAcrossOctaves(t)
            % ⚠ this is what separates a wavelet tiling from a hard one: 12.5, not 1
            r = t.T.redundancyFrame(t.oc);
            t.verifyLessThan(std(r)/mean(r), 0.02);
            t.verifyGreaterThan(mean(r), 10);
        end
        function theSelectionLatticeIsASubLatticeOfTheFrameLattice(t)
            C = rheome.select.catalog();
            h = find(C.recording_id == string(rheomeTestSubject()) & C.bank == "frame" & C.store == "default", 1);
            db = rheome.select.open(char(C.file(h)));
            a = rheome.selection.wavelettile(db, Bands=6, Positions=true);
            f = rheome.selection.wavelettile(db, Bands=6, Positions=true, Stride="frame");
            t.verifyEqual(height(f)/height(a), a.redundancyFrame(1), 'RelTol', 0.02);
        end
        function itIsRegisteredAsAProductKernelThatDoesNotMerge(t)
            S = rheome.selection.registry();  r = S(S.id == "wavelet_tile", :);
            t.verifyEqual(r.domain_kind, "product");
            t.verifyEqual(r.axes{1}, {'time','frequency'});
            t.verifyEqual(r.kind, "kernel");
            t.verifyFalse(r.merges);
            t.verifyTrue(r.nested);          % the two-to-one level nesting above
        end
        function theSupportIsAnEnergyThresholdNotAWidth(t)
            % ⚠ 5.92 / 8.31 / 15.05 cycles at 0.95 / 0.99 / 0.999 -- any support figure is a
            % statement about a threshold, so the threshold travels with it
            C = rheome.select.catalog();
            h = find(C.recording_id == string(rheomeTestSubject()) & C.bank == "frame" & C.store == "default", 1);
            db = rheome.select.open(char(C.file(h)));
            a = rheome.selection.wavelettile(db, Bands=6, Energy=0.95);
            b = rheome.selection.wavelettile(db, Bands=6, Energy=0.999);
            t.verifyLessThan(a.supportCycles, b.supportCycles);
            t.verifyEqual(a.energy, 0.95);
        end
    end
end
