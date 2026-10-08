classdef tIngestCouple < matlab.unittest.TestCase
% Phase-amplitude coupling as a FIRST-CLASS MERGEABLE feature.
%
% The peak had to live at one level because it does not merge. Coupling does: the complex
% sum of fast amplitude times the slow unit phasor, and the plain sum of that amplitude, are
% both sums, so the vector length and the preferred phase come back exactly at any coarser
% level. What is asserted is that the estimator finds a planted coupling and not a control,
% that the phase it reports is the one that was planted, and that merging is exact rather
% than approximate.
%
% Author: Diellor Basha, 2026

    properties
        fs = 600; nT; t; x; xc; cfg; fb; bands; frame; g; s; q
        fslow = 6;  ffast = 60;  atPhase = pi/2      % the slow ANALYTIC phase of the burst
    end

    methods (TestClassSetup)
        function signal(tc)
            tc.nT = 120 * tc.fs;
            tc.t = (0:tc.nT-1)' / tc.fs;
            rng(29);
            slow = sin(2*pi*tc.fslow*tc.t);
            % ⚠ the burst peaks where the slow ARGUMENT is pi; the analytic phase of a sine
            % leads its argument by pi/2, so the preferred phase to expect is pi/2.
            amp = 1 + cos(2*pi*tc.fslow*tc.t - pi);
            tc.x  = 2e-13*slow + 0.5e-13*amp.*sin(2*pi*tc.ffast*tc.t) + 0.2e-13*randn(tc.nT,1);
            tc.xc = 2e-13*slow + 0.5e-13*sin(2*pi*tc.ffast*tc.t)      + 0.2e-13*randn(tc.nT,1);
            tc.cfg = rheome.ingest.config();
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            tc.s = find(tc.bands.fLo == 4 & strcmp(tc.bands.kind, 'octave'));
            tc.q = find(tc.bands.fLo == 32 & strcmp(tc.bands.kind, 'octave'));
        end
    end

    methods (Test)

        function aPairIsAccumulatedAtTheSlowBandsOwnLevel(tc)
            P = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            tc.verifyNumElements(P, 1);
            tc.verifyEqual(P.level, tc.bands.naturalLevel(tc.s));
            tc.verifySize(P.vec, [tc.g.K(P.level+1), 1]);
            tc.verifyTrue(~isreal(P.vec));
            tc.verifyTrue(all(P.amp > 0));
            % ⭐ the same number of slow cycles in every band: that is the point
            L = rheome.select.ladder(struct('bands', tc.bands, 'grid', tc.g, 'meta', ...
                struct('fs', tc.fs, 'nT', tc.nT, 'duration', tc.nT/tc.fs, 'cfg', tc.cfg, 'C', 1)));
            oct = L.kind == "octave" & ~tc.bands.partial;
            tc.verifyEqual(max(L.cycles_per_tile(oct)) / min(L.cycles_per_tile(oct)), 1, 'RelTol', 1e-9);
        end

        function autoPairsAreEveryOctavePairFarEnoughApart(tc)
            [~, pr] = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            tc.verifyEqual(pr, [tc.s tc.q]);
            [P, pa] = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g);
            oc = log2(tc.bands.fCenter(pa(:,2)) ./ tc.bands.fCenter(pa(:,1)));
            tc.verifyTrue(all(oc >= 2));                       % slow first, far enough apart
            tc.verifyEqual(size(pa, 1), sum(arrayfun(@(e) size(e.pairs,1), P)));
            tc.verifyEqual(numel(unique(pa, 'rows')), numel(pa));
        end

        function anAdjacentOctavePairIsRefusedByName(tc)
            % ⚠ the fast band's rate cannot carry twice the slow band's top edge there, so
            % the product would alias; the refusal says the octaves and the requirement.
            near = find(tc.bands.fLo == 8 & strcmp(tc.bands.kind, 'octave'));
            tc.verifyError(@() rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s near]), ...
                           'ingest:couple:pairs');
        end

        function thePlantedCouplingIsFoundAndAControlIsNot(tc)
            P  = rheome.ingest.couple(tc.x,  tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            Pc = rheome.ingest.couple(tc.xc, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            mvl  = median(abs(P.vec) ./ P.amp);
            mvlc = median(abs(Pc.vec) ./ Pc.amp);
            tc.verifyGreaterThan(mvl, 0.1);
            tc.verifyLessThan(mvlc, 0.06);
            tc.verifyGreaterThan(mvl / mvlc, 3);
        end

        function thePreferredPhaseIsTheOneThatWasPlanted(tc)
            P = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            ph = angle(sum(P.vec));                            % the whole record's resultant
            tc.verifyEqual(ph, tc.atPhase, 'AbsTol', 0.25);
        end

        function theAccumulatorsMergeExactly(tc)
            % ⭐ THE PROPERTY THE PEAK DOES NOT HAVE. A parent's sums are its children's, so a
            % coarser estimate is exact -- and rests on twice as many slow cycles.
            P = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q]);
            L = P.level;
            % ⚠ a coarser accumulation is the SAME level index on a doubled floor, which is
            % twice the tile; keeping the level and halving it would be the same tile again.
            gL = rheome.ingest.grid(tc.nT, tc.fs, rheome.ingest.config(FrameFloor = tc.cfg.FrameFloor * 2));
            D = rheome.ingest.couple(tc.x, tc.fb, tc.bands, gL, Pairs=[tc.s tc.q]);
            tc.verifyEqual(gL.tExtent(D.level+1), 2 * tc.g.tExtent(L+1), 'RelTol', 1e-12);
            v = P.vec(1:2:end-1) + P.vec(2:2:end);
            a = P.amp(1:2:end-1) + P.amp(2:2:end);
            n = min(numel(v), numel(D.vec));
            tc.verifyEqual(v(1:n), D.vec(1:n), 'RelTol', 1e-9);
            tc.verifyEqual(a(1:n), D.amp(1:n), 'RelTol', 1e-9);
            tc.verifyEqual(L, P.level);
        end

        function theModulogramIsSumsToo(tc)
            P = rheome.ingest.couple(tc.x, tc.fb, tc.bands, tc.g, Pairs=[tc.s tc.q], PhaseBins=18);
            tc.verifySize(P.hist, [tc.g.K(P.level+1), 1, 18]);
            h = squeeze(sum(P.hist, 1));
            tc.verifyEqual(sum(h), sum(P.amp), 'RelTol', 1e-9);   % the bins partition the sum
            [~, im] = max(h);
            phase = -pi + (im - 0.5) * 2*pi/18;
            tc.verifyEqual(phase, tc.atPhase, 'AbsTol', 0.4);
            p = h / sum(h);
            mi = (log(18) + sum(p .* log(max(p, realmin)))) / log(18);   % Tort's index
            tc.verifyGreaterThan(mi, 0.005);
        end

        function theStoreKeepsTheSumsAndSelectReadsThemBack(tc)
            fx = selFixture();  db = fx.db;
            T = rheome.select.coupling(db, Channels=1);
            tc.verifyGreaterThan(height(T), 0);
            tc.verifyTrue(all(T.octaves >= 2));
            tc.verifyTrue(all(T.mvl >= 0 & T.mvl <= 1));
            tc.verifyTrue(all(abs(T.pref_phase) <= pi));
            p = [T.slow_band(1) T.fast_band(1)];
            fbx = rheome.ingest.bank(fx.nT, fx.fs, fx.cfg);          % the fixture's own bank
            P = rheome.ingest.couple(fx.X(:, 1), fbx, db.bands, db.grid, Pairs=p);
            R = rheome.select.coupling(db, Channels=1, Pairs=p);
            tc.verifyEqual(R.mvl, abs(P.vec) ./ P.amp, 'RelTol', 1e-5);
            tc.verifyEqual(R.pref_phase, angle(P.vec), 'AbsTol', 1e-5);
            tc.verifyEqual(unique(R.n_slow_cycles), db.grid.tExtent(R.level(1)+1) * db.bands.fCenter(p(1)), 'RelTol', 1e-9);
        end

        function readingItAtACoarserLevelSumsTheSums(tc)
            fx = selFixture();  db = fx.db;
            T = rheome.select.coupling(db, Channels=1);
            p = [T.slow_band(1) T.fast_band(1)];
            home = db.bands.naturalLevel(p(1));
            A = rheome.select.coupling(db, Channels=1, Pairs=p, Level=home);
            B = rheome.select.coupling(db, Channels=1, Pairs=p, Level=home+1);
            v = A.amp .* A.mvl .* exp(1i * A.pref_phase);
            vm = v(1:2:end-1) + v(2:2:end);  am = A.amp(1:2:end-1) + A.amp(2:2:end);
            n = min(numel(vm), height(B));
            tc.verifyEqual(B.mvl(1:n), abs(vm(1:n)) ./ am(1:n), 'AbsTol', 1e-6);
            tc.verifyEqual(B.n_slow_cycles(1), 2 * A.n_slow_cycles(1), 'RelTol', 1e-9);
            tc.verifyError(@() rheome.select.coupling(db, Pairs=p, Level=home-1), 'select:coupling:level');
        end

        function theReadFiltersOnChannelsPairsAndWindow(tc)
            fx = selFixture();  db = fx.db;
            all_ = rheome.select.coupling(db);
            one = rheome.select.coupling(db, Channels=2);
            tc.verifyEqual(unique(one.channel_id), 2);
            tc.verifyLessThan(height(one), height(all_));
            w = rheome.select.coupling(db, Channels=2, Window=[5 9]);
            % a tile is kept when it INTERSECTS the window, so a 16 s tile covering 5-9 s
            % comes back with a centre outside it -- which is the tile that covers the window.
            tc.verifyTrue(all(w.t_center - w.t_extent/2 < 9 & w.t_center + w.t_extent/2 > 5));
            tc.verifyLessThanOrEqual(height(w), height(one));
            tc.verifyGreaterThan(height(w), 0);
        end

    end
end

% Author: Diellor Basha, 2026
