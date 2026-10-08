classdef tSelectLadder < matlab.unittest.TestCase
% SELECT.LADDER -- the constant-Q ladder: tile, cycles and rate per band.
%
% The two invariants are the reason the store is built this way, so they are asserted rather
% than described: cycles per tile is the same in every octave, and so is the number of the
% band's OWN samples per tile, while the record's samples per tile double with every level.
%
% Author: Diellor Basha, 2026

    properties
        fx; db; T
    end

    methods (TestClassSetup)
        function open(tc)
            tc.fx = selFixture();
            tc.db = tc.fx.db;
            tc.T  = rheome.select.ladder(tc.db, Scope="channel");
        end
    end

    methods (Test)

        function itHasOneRowPerBandKeyedLikeTheBandRelation(tc)
            tc.verifyEqual(height(tc.T), height(tc.db.bands));
            tc.verifyEqual(tc.T.band_id, tc.db.bands.j);
            tc.verifyEqual(tc.T.f_lo, tc.db.bands.fLo);
            tc.verifyEqual(tc.T.natural_level, tc.db.bands.naturalLevel);
            tc.verifyTrue(all(ismember(["band_id","tile_s","cycles_per_tile","rate_hz", ...
                "samples_per_tile","samples_full","floor_level"], string(tc.T.Properties.VariableNames))));
        end

        function theTileIsTheFirstLevelAtLeastAsLongAsTheSupport(tc)
            g = tc.db.grid;
            tc.verifyTrue(all(tc.T.tile_s >= tc.T.support_s - 1e-12));
            fine = max(tc.T.natural_level - 1, 0);
            tc.verifyTrue(all(g.tExtent(fine + 1)' < tc.T.support_s | tc.T.natural_level == 0));
            tc.verifyEqual(tc.T.tile_s, g.tExtent(tc.T.natural_level + 1)');
        end

        function cyclesPerTileIsTheSameInEveryOctave(tc)
            % ⭐ Constant Q: halve the frequency, double the support, double the tile. The
            % product tile * f_center is therefore an invariant of the bank, not of the band.
            oct = tc.T.kind == "octave" & ~tc.db.bands.partial;
            tc.assumeGreaterThan(sum(oct), 2);
            c = tc.T.cycles_per_tile(oct);
            tc.verifyEqual(max(c) / min(c), 1, 'RelTol', 1e-9);
            tc.verifyEqual(c, tc.T.tile_s(oct) .* tc.T.f_center(oct), 'RelTol', 1e-12);
        end

        function theBandsOwnCostPerTileIsFlatWhileTheRecordsDoubles(tc)
            % ⭐ The multirate claim. An octave costs the same to describe wherever it sits;
            % the record's samples for the same tile double with every level.
            oct = tc.T.kind == "octave" & ~tc.db.bands.partial;
            tc.assumeGreaterThan(sum(oct), 2);
            s = tc.T.samples_per_tile(oct);
            tc.verifyLessThan(max(s) / min(s), 1.3);                      % flat to rounding
            full = tc.T.samples_full(oct);
            tc.verifyEqual(full, tc.T.tile_s(oct) * tc.db.meta.fs, 'RelTol', 1e-12);
            [~, i] = sort(tc.T.natural_level(oct));
            r = full(i);
            tc.verifyEqual(r(2:end) ./ r(1:end-1), 2 * ones(numel(r)-1, 1), 'RelTol', 1e-12);
            tc.verifyGreaterThan(max(full) / max(s), 50);                 % the whole point
        end

        function theRateIsTheOneTheBankActuallyUses(tc)
            % Not a nominal bandwidth: the rate a member's coefficients arrive at, which is
            % what wt returns. Checked against the bank itself.
            m = tc.db.meta;  cfg = m.cfg;
            fb = rheome.timefilterbank(m.nT, 'SamplingFrequency', m.fs, 'VoicesPerOctave', cfg.VoicesPerOctave, ...
                                'Anchor', cfg.Anchor, 'FrequencyLimits', cfg.FrequencyLimits, 'Oversample', cfg.Oversample);
            r = rates(fb);
            C = wt(fb, tc.fx.X(:, 1));
            tc.verifyEqual(r(:), [C.rate]', 'RelTol', 1e-12);             % rates == wt's rates
            for j = 1:height(tc.T)
                mem = tc.db.bands.scales{j};
                tc.verifyEqual(tc.T.rate_hz(j), max(r(mem)), 'RelTol', 1e-12);
            end
            tc.verifyTrue(all(tc.T.rate_hz > 0));
        end

        function aMorseStoreReportsNoRateRatherThanAWrongOne(tc)
            % The Morse path does not evaluate a band at its own rate, so the column is NaN
            % instead of a number that would mean something else.
            db2 = tc.db;  db2.meta.cfg.Bank = 'morse';
            T2 = rheome.select.ladder(db2);
            tc.verifyTrue(all(isnan(T2.rate_hz)));
            tc.verifyTrue(all(isnan(T2.samples_per_tile)));
            tc.verifyEqual(T2.tile_s, tc.T.tile_s);                       % the tiling is unchanged
        end

        function theChannelFloorColumnSaysWhatOneSensorCanActuallyBeRead(tc)
            tc.verifyTrue(all(tc.T.floor_level >= tc.T.natural_level));
            tc.verifyEqual(tc.T.floor_level, max(tc.T.natural_level, tc.db.channelLevel));
            T3 = rheome.select.ladder(tc.db);                                    % Scope none: no floor columns
            tc.verifyFalse(ismember('floor_level', T3.Properties.VariableNames));
        end

    end
end

% Author: Diellor Basha, 2026
