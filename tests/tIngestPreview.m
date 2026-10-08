classdef tIngestPreview < matlab.unittest.TestCase
% rheome.ingest.preview: a min/max pyramid beside a recording, with no transform at all.
%
% The preview exists so a recording can be LOOKED at before anything is loaded: browse the
% channels, pick a stretch, and read off what that stretch can resolve. It is a tile store
% with no bands, so the things asserted here are that it contains the signal at every level,
% that it holds nothing it does not need, and that every reader says what it is rather than
% failing on a missing array.
%
% Author: Diellor Basha, 2026

    properties
        dir; rec; file; db; X; fs = 200; nT = 4000; C = 3
    end

    methods (TestClassSetup)
        function build(tc)
            rng(23);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13 * randn(tc.nT, tc.C) + 5e-13 * sin(2*pi*7*t) .* (t > 4 & t < 8);
            tc.X(1500, 1) = 9e-12;                                % a spike, to be contained
            tc.dir = tempname;  mkdir(tc.dir);
            tc.addTeardown(@() rmdir(tc.dir, 's'));
            tc.rec = fullfile(tc.dir, 'recording.mat');
            F = tc.X.';  Time = t';  sfreq = tc.fs;                                    %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);  ChannelName = {'A','B','C'};                 %#ok<NASGU>
            ChannelType = repmat({'MEG'}, 1, tc.C);                                    %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;  Comment = 'p';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(tc.rec, 'F','Time','sfreq','ChannelFlag','ChannelName','ChannelType', ...
                 'nCh','nT','Comment','Events','nAvg','source', '-v7.3', '-nocompression');
            tc.file = rheome.ingest.preview(tc.rec, Floor=0.05, Verbose=false);
            tc.db = rheome.select.open(tc.file);
        end
    end

    methods (Test)

        function itHoldsMinMaxAndCountsAndNothingElse(tc)
            w = whos('-file', tc.file);
            nm = string({w.name});
            lv = nm(startsWith(nm, 'L'));
            suffix = unique(extractAfter(lv, 4));
            tc.verifyEqual(sort(suffix(:)), ["max"; "min"; "n"]);
            tc.verifyEqual(numel(lv), 3 * (tc.db.grid.Lmax + 1));
            tc.verifyTrue(all(ismember(["meta","bands","frame","grid"], nm)));
            tc.verifyEmpty(tc.db.bands);
            tc.verifyTrue(tc.db.preview);
            tc.verifyEqual(tc.db.kind, 'preview');
        end

        function everyTileContainsItsSamplesAtEveryLevel(tc)
            % ⭐ The one guarantee a preview must keep: what you see bounds what is there.
            for L = 0:tc.db.grid.Lmax
                R = rheome.select.derive(tc.db, Level=L, Stats=["min","max"]);
                for c = 1:tc.C
                    r = R(R.unit_id == c, :);
                    for k = 1:height(r)
                        a = round(r.t_lo(k) * tc.fs) + 1;
                        b = min(round(r.t_hi(k) * tc.fs), tc.nT);
                        if b < a, continue; end
                        x = tc.X(a:b, c);
                        tc.verifyLessThanOrEqual(r.min(k), min(x), sprintf('L%d c%d k%d', L, c, k));
                        tc.verifyGreaterThanOrEqual(r.max(k), max(x), sprintf('L%d c%d k%d', L, c, k));
                    end
                end
            end
        end

        function theSpikeSurvivesToTheTopLevel(tc)
            % A subsampled preview would step over it; an extremum pyramid cannot.
            R = rheome.select.derive(tc.db, Level=tc.db.grid.Lmax, Channels=1, Stats="max");
            tc.verifyGreaterThanOrEqual(R.max(1), 9e-12);
        end

        function aParentIsTheExtremumOfItsChildren(tc)
            for L = 1:4
                P = rheome.select.derive(tc.db, Level=L, Channels=2, Stats=["min","max"]);
                C = rheome.select.derive(tc.db, Level=L-1, Channels=2, Stats=["min","max"]);
                for k = 1:height(P)
                    i = [2*k-1, 2*k];  i = i(i <= height(C));
                    tc.verifyEqual(P.min(k), min(C.min(i)), sprintf('L%d k%d', L, k));
                    tc.verifyEqual(P.max(k), max(C.max(i)), sprintf('L%d k%d', L, k));
                end
            end
        end

        function theCountsAddUpToTheRecord(tc)
            for L = [0 3 tc.db.grid.Lmax]
                R = rheome.select.rows(tc.db, 'tile', Level=L);
                tc.verifyEqual(sum(double(R.n)), tc.nT, sprintf('level %d', L));
            end
        end

        function onlyMinMaxAndPtpCanBeDerivedAndTheRestSaysWhy(tc)
            R = rheome.select.derive(tc.db, Level=3, Channels=1, Stats="all");
            got = intersect(string(R.Properties.VariableNames), ["min","max","ptp"]);
            tc.verifyEqual(sort(got(:)), ["max"; "min"; "ptp"]);
            tc.verifyFalse(ismember("rms", string(R.Properties.VariableNames)));
            tc.verifyError(@() rheome.select.derive(tc.db, Level=3, Stats="rms"), 'select:derive:preview');
            tc.verifyError(@() rheome.select.derive(tc.db, Level=3, Stats="share"), 'select:derive:preview');
            try
                rheome.select.derive(tc.db, Level=3, Stats="spectralEntropy");
            catch err
                tc.verifySubstring(err.message, 'preview store');
            end
        end

        function rebuildingIsANoOpAndADifferentFloorIsRefused(tc)
            f2 = rheome.ingest.preview(tc.rec, Floor=0.05, Verbose=false);
            tc.verifyEqual(f2, tc.file);
            f3 = rheome.ingest.preview(tc.rec, Floor=0.1, Verbose=false);      % its own file
            tc.verifyNotEqual(f3, tc.file);
            d3 = rheome.select.open(f3);
            tc.verifyEqual(d3.grid.tExtent(1), 0.1);
            tc.verifyError(@() rheome.ingest.preview(tc.rec, Floor=0.05, Channels=[1 2], Verbose=false), ...
                           'ingest:preview:hash');
        end

        function theBaseOptionWritesOneLevelOnly(tc)
            f = rheome.ingest.preview(tc.rec, Floor=0.2, Levels="base", Verbose=false, Overwrite=true);
            w = whos('-file', f);
            nm = string({w.name});
            tc.verifyEqual(sum(startsWith(nm, 'L')), 3);                 % n, min, max at level 0
            d = rheome.select.open(f);
            tc.verifyEqual(height(rheome.select.derive(d, Level=0, Channels=1, Stats="ptp")), d.grid.K(1));
        end

        function aWindowLengthNamesItsBandSupport(tc)
            % ⭐ What the explorer is for: a stretch chosen in the preview already says which
            % bands can be resolved in it.
            T = rheome.ingest.support([0.25 1 2 4 16]);
            tc.verifyEqual(T.window_s, [0.25 1 2 4 16]');
            tc.verifyTrue(all(diff(T.f_lowest) < 0));                    % longer window, lower floor
            tc.verifyEqual(T.f_lowest .* T.window_s, repmat(15, 5, 1), 'RelTol', 1e-12);
            tc.verifyTrue(all(T.band_lo >= T.f_lowest));                 % the octave fully inside
            tc.verifyTrue(all(T.band_hi == 2 * T.band_lo));
            tc.verifyEqual(T.cycles_in_window, repmat(T.cycles_in_window(1), 5, 1), 'RelTol', 1e-9);
            t1 = rheome.ingest.support(2, Voices=10);
            tc.verifyGreaterThan(t1.f_lowest, T.f_lowest(3));            % more voices, longer support
        end

        function theSupportMatchesTheStoresOwnLadder(tc)
            % The same relation the bank measured. ⚠ A band's stored support is k over its
            % LOW EDGE, not its centre (rheome.ingest.bank takes the band's slowest member), so a
            % window that just fits a band resolves from that band's low edge up -- which is
            % exactly what the octave column says.
            fx = selFixture();
            L = rheome.select.ladder(fx.db);
            oct = L.kind == "octave" & ~fx.db.bands.partial;
            S = rheome.ingest.support(L.support_s(oct));
            tc.verifyEqual(S.f_lowest, L.f_lo(oct), 'RelTol', 0.02);
            tc.verifyEqual(S.band_lo, L.f_lo(oct), 'RelTol', 1e-9);
        end

    end
end

% Author: Diellor Basha, 2026
