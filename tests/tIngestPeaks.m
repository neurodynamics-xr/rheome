classdef tIngestPeaks < matlab.unittest.TestCase
% The features that are NOT closed form: the per-tile spectral peak, and the shape moments.
%
% Two families arrive together here. The shape moments (sum|x|, sum sqrt|x|, sum x^3,
% sum x^4) are sums, so they merge like every other moment and give skewness, kurtosis and
% the shape/impulse/clearance factors at every level. The spectral peak does NOT merge -- and
% does not have to, because the diagonal gives every band exactly one home level, so it is
% computed there once. What is asserted is that both are right against planted signals.
%
% Author: Diellor Basha, 2026

    properties
        fs = 600; nT; t; x; cfg; fb; bands; frame; g
        f1 = 10.3;  a1 = 3e-13;         % an alpha tone
        f2 = 22.7;  a2 = 1e-13;         % a beta tone
    end

    methods (TestClassSetup)
        function signal(tc)
            tc.nT = 120 * tc.fs;
            tc.t = (0:tc.nT-1)' / tc.fs;
            rng(17);
            tc.x = 0.4e-13*randn(tc.nT,1) + tc.a1*sin(2*pi*tc.f1*tc.t) + tc.a2*sin(2*pi*tc.f2*tc.t);
            tc.cfg = rheome.ingest.config();
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
        end
    end

    methods (Test)

        function everyBandGetsItsPeakAtItsOwnLevelAndNowhereElse(tc)
            P = rheome.ingest.peaks(tc.x, tc.bands, tc.g);
            lv = [P.level];
            tc.verifyEqual(sort(lv), sort(unique(tc.bands.naturalLevel)'));
            all_ = [];
            for i = 1:numel(P)
                tc.verifyEqual(sort(P(i).bands), sort(find(tc.bands.naturalLevel == P(i).level)'));
                tc.verifySize(P(i).freq, [tc.g.K(P(i).level+1), numel(P(i).bands)]);
                all_ = [all_, P(i).bands];                       %#ok<AGROW>
            end
            tc.verifyEqual(sort(all_), (1:height(tc.bands)));     % each band exactly once
        end

        function thePeakFindsAPlantedToneToAFractionOfABin(tc)
            % ⭐ Why FFT and not the bank: the bank's resolution is one member (1.19x at four
            % voices); a 2 s tile's own spectrum has 0.5 Hz bins and the parabolic fit puts
            % the peak between them, which is what makes an alpha peak quotable.
            P = rheome.ingest.peaks(tc.x, tc.bands, tc.g);
            for pair = {{tc.f1, tc.a1, 8}, {tc.f2, tc.a2, 16}}
                f0 = pair{1}{1};  a0 = pair{1}{2};  lo = pair{1}{3};
                b = find(tc.bands.fLo == lo & strcmp(tc.bands.kind, 'octave'));
                e = P([P.level] == tc.bands.naturalLevel(b));
                j = find(e.bands == b, 1);
                fpk = median(e.freq(:, j), 'omitnan');
                apk = median(e.amp(:, j), 'omitnan');
                bin = tc.fs / (tc.g.F * 2^e.level);
                tc.verifyEqual(fpk, f0, 'AbsTol', 0.25 * bin, sprintf('%g Hz', f0));
                tc.verifyEqual(apk, a0, 'RelTol', 0.10, sprintf('%g amplitude', f0));
            end
        end

        function theWindowedPowerAgreesWithTheBanksExactEnergy(tc)
            % ⚠ These are two different estimates on purpose: `energy` is exact by Parseval
            % over the whole record, fftPower is a windowed tile's own. On stationary content
            % they agree to a few per cent, and the test pins that rather than pretending
            % they are the same number.
            T0 = rheome.ingest.reduce(tc.x, tc.fb, tc.bands, tc.g, tc.frame);
            T = rheome.ingest.rollup(T0, tc.g);
            P = rheome.ingest.peaks(tc.x, tc.bands, tc.g);
            b = find(tc.bands.fLo == 8 & strcmp(tc.bands.kind, 'octave'));
            L = tc.bands.naturalLevel(b);
            e = P([P.level] == L);
            bankPower = squeeze(T{L+1}.energy(:, 1, b)) ./ double(T{L+1}.n);
            fftPower = e.power(:, e.bands == b);
            ok = isfinite(fftPower) & bankPower > 0;
            r = median(fftPower(ok) ./ bankPower(ok));
            tc.verifyEqual(r, 1, 'RelTol', 0.15);
        end

        function theShapeMomentsAreSumsAndThereforeMerge(tc)
            T0 = rheome.ingest.reduce(tc.x, tc.fb, tc.bands, tc.g, tc.frame);
            tc.verifyTrue(all(isfield(T0, {'sumAbs','sumSqrt','sumX3','sumX4'})));
            T = rheome.ingest.rollup(T0, tc.g);
            for L = [1 3 5]
                gL = rheome.ingest.grid(tc.nT, tc.fs, rheome.ingest.config(FrameFloor = tc.cfg.FrameFloor * 2^L));
                D = rheome.ingest.reduce(tc.x, tc.fb, tc.bands, gL, tc.frame);
                for f = ["sumAbs","sumSqrt","sumX3","sumX4"]
                    tc.verifyEqual(T{L+1}.(f), D.(f), 'RelTol', 1e-12, sprintf('%s at level %d', f, L));
                end
            end
        end

        function theShapeFactorsMatchTheirDefinitionsOnKnownSignals(tc)
            % A sine has crest sqrt(2), shape factor pi/(2 sqrt 2) and kurtosis 1.5; Gaussian
            % noise has kurtosis 3 and skewness 0. Checked on one tile of each.
            n = 6000;  tt = (0:n-1)'/tc.fs;
            for pair = {{sin(2*pi*10*tt), sqrt(2), pi/(2*sqrt(2)), 1.5, 0}, ...
                        {randn(n,1),      NaN,     NaN,            3.0, 0}}
                y = pair{1}{1};
                m = struct('n', n, 'sumX', sum(y), 'sumX2', sum(y.^2), 'absMax', max(abs(y)), ...
                           'sumAbs', sum(abs(y)), 'sumSqrt', sum(sqrt(abs(y))), ...
                           'sumX3', sum(y.^3), 'sumX4', sum(y.^4));
                mu = m.sumX/n;  r2 = m.sumX2/n;  sg = sqrt(max(r2 - mu^2, 0));
                m3 = m.sumX3/n - 3*mu*r2 + 2*mu^3;
                m4 = m.sumX4/n - 4*mu*m.sumX3/n + 6*mu^2*r2 - 3*mu^4;
                if ~isnan(pair{1}{2})
                    tc.verifyEqual(m.absMax/sqrt(r2), pair{1}{2}, 'RelTol', 0.01);
                    tc.verifyEqual(sqrt(r2)/(m.sumAbs/n), pair{1}{3}, 'RelTol', 0.01);
                end
                tc.verifyEqual(m4/sg^4, pair{1}{4}, 'RelTol', 0.06);
                tc.verifyEqual(m3/sg^3, pair{1}{5}, 'AbsTol', 0.12);
            end
        end

        function theStoreKeepsThePeaksAtTheHomeLevelAndSelectReadsThemThere(tc)
            % The round trip: what rheome.ingest.peaks computes is what the store holds and what
            % rheome.select.derive gives back, at the level the diagonal assigns and nowhere else.
            fx = selFixture();  db = fx.db;
            P = rheome.ingest.peaks(fx.X(:, 1), db.bands, db.grid);
            for i = 1:numel(P)
                L = P(i).level;
                R = rheome.select.derive(db, Level=L, Channels=1, Stats="fft", Bands=P(i).bands);
                tc.verifyEqual(double(R.peakFreq), P(i).freq, 'RelTol', 1e-6, sprintf('level %d freq', L));
                tc.verifyEqual(double(R.peakAmp),  P(i).amp,  'RelTol', 1e-6, sprintf('level %d amp', L));
                tc.verifyEqual(double(R.fftPower), P(i).power, 'RelTol', 1e-6, sprintf('level %d power', L));
                tc.verifyEqual(R.Properties.CustomProperties.bands, P(i).bands);
            end
        end

        function askingForAPeakAwayFromItsHomeLevelSaysWhereItLives(tc)
            fx = selFixture();  db = fx.db;
            b = find(db.bands.fLo == 8 & strcmp(db.bands.kind, 'octave'));
            home = db.bands.naturalLevel(b);
            tc.verifyNotEmpty(rheome.select.derive(db, Level=home, Channels=1, Stats="peakFreq", Bands=b));
            tc.verifyError(@() rheome.select.derive(db, Level=home+1, Channels=1, Stats="peakFreq", Bands=b), ...
                           'select:peaks:level');
            try
                rheome.select.derive(db, Level=home+1, Channels=1, Stats="peakFreq", Bands=b);
            catch err
                tc.verifySubstring(err.message, sprintf('level %d', home));
                tc.verifySubstring(err.message, 'not mergeable');
            end
            tc.verifyError(@() rheome.select.derive(db, Level=home, Scope="group", Stats="peakFreq", Bands=b), ...
                           'select:peaks:scope');
        end

        function theNamedBandsSayHowWellTheOctavesCoverThem(tc)
            % ⚠ The store's bands are octaves on an absolute grid; delta/theta/alpha/beta/
            % gamma are conventions with ragged edges. The mapping records the approximation
            % rather than hiding it.
            fx = selFixture();
            T = rheome.select.named(fx.db);
            tc.verifyTrue(all(ismember(["delta","theta","alpha"], T.name)));
            a = T(T.name == "alpha", :);
            tc.verifyEqual(a.f_lo, 8);  tc.verifyEqual(a.f_hi, 13);
            tc.verifyEqual(a.band_lo, 8);  tc.verifyEqual(a.band_hi, 16);
            tc.verifyEqual(a.covered, (13-8)/(16-8), 'RelTol', 1e-12);
            tc.verifyTrue(all(T.covered > 0 & T.covered <= 1));
            for i = 1:height(T)
                tc.verifyTrue(all(ismember(T.band_ids{i}, fx.db.bands.j)));
            end
        end

        function aWindowChoiceIsHonouredAndZeroPaddingSharpensNothingItShouldNot(tc)
            Ph = rheome.ingest.peaks(tc.x, tc.bands, tc.g, Window="hann");
            Pr = rheome.ingest.peaks(tc.x, tc.bands, tc.g, Window="rect", Pad=1);
            b = find(tc.bands.fLo == 8 & strcmp(tc.bands.kind, 'octave'));
            L = tc.bands.naturalLevel(b);
            fh = median(Ph([Ph.level] == L).freq(:, Ph([Ph.level] == L).bands == b), 'omitnan');
            fr = median(Pr([Pr.level] == L).freq(:, Pr([Pr.level] == L).bands == b), 'omitnan');
            tc.verifyEqual(fh, tc.f1, 'AbsTol', 0.2);
            tc.verifyEqual(fr, tc.f1, 'AbsTol', 0.5);            % a rectangular window leaks more
        end

    end
end

% Author: Diellor Basha, 2026
