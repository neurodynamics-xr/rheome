classdef tSelectDerive < matlab.unittest.TestCase
% SELECT.DERIVE -- the scalarization layer: every column a function of the stored moments.
%
% The point of these is that nothing here is stored, so each column must be reproducible from
% the relations alone (rheome.select.rows), and the identities the constant-Q construction promises
% (shares sum to one, the tight frame partitions energy, a parent's power is its children's
% n-weighted mean) must hold to floating point.
%
% Author: Diellor Basha, 2026

    properties
        fx; db
    end

    methods (TestClassSetup)
        function open(tc)
            tc.fx = selFixture();
            tc.db = tc.fx.db;
        end
    end

    methods (Test)

        function theVocabularyDescribesItself(tc)
            % A caller (and the browser's menu) must be able to ask what exists without
            % reading the source: name, kind, the moments it needs, the formula.
            V = rheome.select.derive();
            tc.verifyClass(V, 'table');
            tc.verifyEqual(V.Properties.VariableNames, {'name','kind','needs','formula'});
            tc.verifyTrue(all(ismember(["rms","crest","spectralEntropy","share"], V.name)));
            tc.verifyTrue(all(ismember(V.kind, ["time","spectral","band","fft"])));
            tc.verifyEqual(numel(unique(V.name)), height(V));
        end

        function theShorthandsExpandToKinds(tc)
            V = rheome.select.derive();
            R = rheome.select.derive(tc.db, Level=4, Channels=1, Stats="time");
            tc.verifyTrue(all(ismember(V.name(V.kind == "time"), string(R.Properties.VariableNames))));
            tc.verifyFalse(ismember("spectralEntropy", string(R.Properties.VariableNames)));
            A = rheome.select.derive(tc.db, Level=4, Channels=1, Stats="all");
            tc.verifyFalse(ismember("share", string(A.Properties.VariableNames)));   % "all" is the SCALARS
            tc.verifyTrue(ismember("spectralEntropy", string(A.Properties.VariableNames)));
            B = rheome.select.derive(tc.db, Level=4, Channels=1, Stats="band");
            tc.verifyTrue(all(ismember(["bandPower","share","bandEnergy","bandEnvMax"], string(B.Properties.VariableNames))));
        end

        function anUnknownStatisticIsRefusedByName(tc)
            tc.verifyError(@() rheome.select.derive(tc.db, Level=2, Stats="spectralWhatsit"), 'select:derive:stat');
            tc.verifyError(@() rheome.select.derive(tc.db, Stats="rms"), 'select:derive:level');
            tc.verifyError(@() rheome.select.derive(tc.db, Level=99), 'select:derive:level');
        end

        function theTimeStatisticsComeOutOfTheStoredMoments(tc)
            % The reference is rheome.select.rows: the relations a database would hold.
            L = 4;
            R = rheome.select.derive(tc.db, Level=L, Stats=["mean","rms","std","power","peak","ptp","crest"]);
            F = rheome.select.rows(tc.db, 'feature', Level=L);
            T = rheome.select.rows(tc.db, 'tile', Level=L);
            n = double(T.n);
            for c = 1:tc.fx.C
                r = R(R.unit_id == c, :);
                f = F(F.channel_id == c, :);
                tc.verifyEqual(r.mean, f.sum_x ./ n, 'RelTol', 1e-12, 'mean');
                tc.verifyEqual(r.rms, sqrt(f.sum_x2 ./ n), 'RelTol', 1e-12, 'rms');
                tc.verifyEqual(r.power, f.sum_x2 ./ n, 'RelTol', 1e-12, 'power');
                tc.verifyEqual(r.peak, double(f.abs_max), 'RelTol', 1e-12, 'peak');
                tc.verifyEqual(r.ptp, double(f.max) - double(f.min), 'RelTol', 1e-12, 'ptp');
                tc.verifyEqual(r.crest, double(f.abs_max) ./ sqrt(f.sum_x2 ./ n), 'RelTol', 1e-12, 'crest');
                sd = sqrt(max(f.sum_x2 ./ n - (f.sum_x ./ n).^2, 0)) .* sqrt(n ./ max(n - 1, 1));
                tc.verifyEqual(r.std, sd, 'RelTol', 1e-9, 'std');
            end
        end

        function aParentsPowerIsItsChildrensWeightedMeanButItsRmsIsNot(tc)
            % ⭐ THE REASON MOMENTS ARE STORED AND NOT SUMMARIES. Power merges exactly as an
            % n-weighted mean; rms does not merge at all (sqrt of a mean is not a mean of
            % sqrts), so a level-of-detail view that averaged fine rms values would show a
            % number the data never had. Both facts are checked here.
            L = 3;
            P = rheome.select.derive(tc.db, Level=L,   Channels=1, Stats=["power","rms"]);
            C = rheome.select.derive(tc.db, Level=L-1, Channels=1, Stats=["power","rms"]);
            K = height(P);
            for k = 1:K
                i = [2*k-1, 2*k];  i = i(i <= height(C));
                nk = C.n(i);  w = nk / sum(nk);
                tc.verifyEqual(P.power(k), sum(w .* C.power(i)), 'RelTol', 1e-12);
            end
            avgRms = arrayfun(@(k) mean(C.rms(min([2*k-1, 2*k], height(C)))), (1:K)');
            tc.verifyGreaterThan(max(abs(P.rms - avgRms) ./ P.rms), 1e-6);   % genuinely different
        end

        function sharesPartitionTheTileAndBandPowerIsEnergyOverN(tc)
            L = 4;
            R = rheome.select.derive(tc.db, Level=L, Channels=2, Stats=["share","bandPower","bandEnergy","energy"]);
            tc.verifyEqual(sum(R.share, 2), ones(height(R), 1), 'AbsTol', 1e-12);
            tc.verifyEqual(R.bandPower, R.bandEnergy ./ R.n, 'RelTol', 1e-12);
            tc.verifyEqual(R.energy, sum(R.bandEnergy, 2), 'RelTol', 1e-12);
            b = R.Properties.CustomProperties.bands;
            tc.verifyEqual(b, tc.db.grid.bandsAt{L+1});
            F = rheome.select.rows(tc.db, 'feature_band', Level=L);
            F = F(F.channel_id == 2, :);
            for i = 1:numel(b)
                tc.verifyEqual(R.bandEnergy(:, i), F.energy(F.band_id == b(i)), 'RelTol', 1e-12);
            end
        end

        function theTightFrameLeavesNoResidualAtTheTop(tc)
            % At Lmax every band is carried and the frame is tight over the record, so the
            % bands account for sumX2. Anywhere below, the residual is signed: a transient's
            % response crosses tile edges (measured here at level Lmax-1: up to 0.22).
            % ⚠ THE TOLERANCE IS THE FIXTURE'S, NOT THE CONSTRUCTION'S. On this 20 s, 200 Hz
            % toy the top-level residual is 1.3e-4, set by the bank's edge members over a very
            % short record; on the reference store (600 s, 600 Hz, 15 bands) the same quantity is
            % 1e-11. The tolerance below pins the fixture, not the claim.
            L = tc.db.grid.Lmax;
            R = rheome.select.derive(tc.db, Level=L, Stats="residual");
            tc.verifyEqual(numel(tc.db.grid.bandsAt{L+1}), height(tc.db.bands));
            tc.verifyEqual(R.residual, zeros(height(R), 1), 'AbsTol', 5e-4);
        end

        function theSpectralStatisticsStayInsideTheirDefinitions(tc)
            L = 4;
            R = rheome.select.derive(tc.db, Level=L, Stats="spectral");
            b = R.Properties.CustomProperties.bands;
            j = arrayfun(@(x) find(tc.db.bands.j == x, 1), b);
            tc.verifyTrue(all(R.spectralEntropy >= 0 & R.spectralEntropy <= 1));
            tc.verifyTrue(all(R.spectralFlatness > 0 & R.spectralFlatness <= 1 + 1e-12));
            tc.verifyTrue(all(R.spectralCrest >= 1 - 1e-12 & R.spectralCrest <= numel(b) + 1e-9));
            tc.verifyTrue(all(R.spectralCentroid >= min(tc.db.bands.fCenter(j)) - 1e-9));
            tc.verifyTrue(all(R.spectralCentroid <= max(tc.db.bands.fCenter(j)) + 1e-9));
            tc.verifyTrue(all(ismember(R.dominantBand, b)));
            tc.verifyTrue(all(R.coneFrac >= 0 & R.coneFrac <= 1));
        end

        function theDominantBandFindsTheFixturesRhythms(tc)
            % Channel 1 carries 10 Hz bursts, channel 2 a continuous 3 Hz rhythm: the
            % dominant band over the record must be the octave each one lives in.
            % ⚠ AT THE TOP LEVEL, WHERE EVERY BAND IS CARRIED. At 4 s tiles the 2-4 Hz octave
            % is below the diagonal, so the dominant band there is the widest noise octave --
            % correct, and not a statement about the rhythm.
            L = tc.db.grid.Lmax;
            R = rheome.select.derive(tc.db, Level=L, Stats=["dominantBand","dominantFreq"]);
            for pair = {{1, 10}, {2, 3}}
                c = pair{1}{1};  f = pair{1}{2};
                r = R(R.unit_id == c, :);
                j = mode(r.dominantBand);
                i = find(tc.db.bands.j == j, 1);
                tc.verifyTrue(tc.db.bands.fLo(i) <= f && f < tc.db.bands.fHi(i), ...
                    sprintf('channel %d: dominant band %d is %.3g-%.3g Hz, not %g Hz', c, j, tc.db.bands.fLo(i), tc.db.bands.fHi(i), f));
            end
        end

        function theWindowKeepsEveryTileThatIntersectsIt(tc)
            % A coarse tile overlapping a short window IS the tile covering that window;
            % keeping only tiles inside the window would return nothing.
            L = 5;  ext = tc.db.grid.tExtent(L+1);
            R = rheome.select.derive(tc.db, Level=L, Channels=1, Stats="rms", Window=[9 10]);
            tc.verifyGreaterThanOrEqual(height(R), 1);
            tc.verifyTrue(all(R.t_lo < 10 & R.t_hi > 9));
            tc.verifyLessThanOrEqual(R.t_hi(1) - R.t_lo(1), ext + 1e-12);
            W = rheome.select.derive(tc.db, Level=L, Channels=1, Stats="rms");
            tc.verifyEqual(height(W), tc.db.grid.K(L+1));
            past = rheome.select.derive(tc.db, Level=L, Channels=1, Stats="rms", Window=[1e6 1e6+1]);
            tc.verifyEqual(height(past), 1);                       % never empty: the nearest tile
        end

        function aNodeCarriesItsSensorCountBecauseItsMomentsAreSums(tc)
            % ⚠ g_sumX2 adds every member channel's squares while n counts SAMPLES, so a
            % node's rms is sqrt(nSensors) times a channel's and its "crest" is not a crest.
            % The count travels with the rows so a caller can normalise.
            L = 4;
            G = rheome.select.derive(tc.db, Level=L, Scope="group", Stats=["rms","crest"]);
            tc.verifyTrue(ismember("n_sensors", string(G.Properties.VariableNames)));
            root = min(tc.db.groupNodes);
            r = G(G.unit_id == root, :);
            tc.verifyEqual(unique(r.n_sensors), tc.fx.C);              % the root holds them all
            Ch = rheome.select.derive(tc.db, Level=L, Stats="rms");
            per = zeros(height(r), 1);
            for c = 1:tc.fx.C, per = per + Ch.rms(Ch.unit_id == c).^2; end
            tc.verifyEqual(r.rms.^2, per, 'RelTol', 1e-10);            % the sum of squares
            % the numerator is one channel's peak while the denominator sums every channel's
            % power, so a node's crest is strictly below the largest crest among its members
            % -- on a 64-sensor node it lands under 1, which no real crest can.
            cc = zeros(height(r), tc.fx.C);
            Ch2 = rheome.select.derive(tc.db, Level=L, Stats="crest");
            for c = 1:tc.fx.C, cc(:, c) = Ch2.crest(Ch2.unit_id == c); end
            tc.verifyLessThan(r.crest, max(cc, [], 2));
            tc.verifyLessThan(median(r.crest), max(cc(:)) / sqrt(tc.fx.C) * 1.5);
            C1 = rheome.select.derive(tc.db, Level=L, Channels=1, Stats="rms");
            tc.verifyEqual(unique(C1.n_sensors), 1);                   % a channel is one sensor
        end

        function groupScopeAddressesTheSensorTreeAndSumsItsMembers(tc)
            % A node's moments are the exact sums over its member channels (the position
            % pyramid), so a node's energy is the sum of its sensors' energies.
            L = 4;
            G = rheome.select.derive(tc.db, Level=L, Scope="group", Stats=["energy","rms"]);
            tc.verifyEqual(unique(G.unit_id)', tc.db.groupNodes(:)');
            tc.verifyEqual(unique(G.scope), "group");
            root = min(tc.db.groupNodes);
            Ch = rheome.select.derive(tc.db, Level=L, Stats="energy");
            tot = zeros(tc.db.grid.K(L+1), 1);
            for c = 1:tc.fx.C, tot = tot + Ch.energy(Ch.unit_id == c); end
            tc.verifyEqual(G.energy(G.unit_id == root), tot, 'RelTol', 1e-10);
        end

        function derivingCostsOnlyTheMomentsTheStatisticNeeds(tc)
            % ⭐ The reason the browser stays cheap: asking for rms must not read the band
            % arrays. Counted through the select handle's cost meter.
            db1 = rheome.select.open(tc.fx.store);
            rheome.select.derive(db1, Level=4, Channels=1, Stats="rms");
            db2 = rheome.select.open(tc.fx.store);
            rheome.select.derive(db2, Level=4, Channels=1, Stats="spectralEntropy");
            tc.verifyLessThan(db1.cost('bytes'), db2.cost('bytes'));     % no band arrays for rms
            tc.verifyLessThanOrEqual(db1.cost('reads'), 2);              % n and sumX2, nothing else
            db3 = rheome.select.open(tc.fx.store);
            rheome.select.derive(db3, Level=4, Channels=1, Stats="all");
            tc.verifyGreaterThan(db3.cost('bytes'), db2.cost('bytes'));
        end

    end
end

% Author: Diellor Basha, 2026
