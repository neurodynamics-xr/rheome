classdef tIngestEnvelope < matlab.unittest.TestCase
% The envelope exemption: min and max below the channel floor, and nothing else.
%
% The position diagonal drops per-channel rows below ChannelMinTile because the band arrays
% are what cost. min and max are the two cheapest columns and the only ones that make a
% level-of-detail view of a TRACE possible, so they are kept at every level. This checks
% that the store really holds them there, that the rest is really absent, and that the query
% side says which is which.
%
% Author: Diellor Basha, 2026

    properties
        dir; rec; file; db; X; fs = 200; nT = 4000; C = 4; cfg
    end

    methods (TestClassSetup)
        function build(tc)
            rng(19);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13 * randn(tc.nT, tc.C) + 4e-13 * sin(2*pi*10*t) .* (t > 5 & t < 9);
            tc.dir = tempname;  mkdir(tc.dir);
            tc.addTeardown(@() rmdir(tc.dir, 's'));
            tc.rec = fullfile(tc.dir, 'recording.mat');
            F = tc.X.';  Time = t';  sfreq = tc.fs;                                   %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);  ChannelName = {'A','B','C','D'};            %#ok<NASGU>
            ChannelType = repmat({'MEG'}, 1, tc.C);                                   %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;  Comment = 'env';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(tc.rec, 'F','Time','sfreq','ChannelFlag','ChannelName','ChannelType', ...
                 'nCh','nT','Comment','Events','nAvg','source', '-v7.3', '-nocompression');
            P = [0 0 0; 0.03 0 0; 0 0.03 0; 0.03 0.03 0];
            tc.cfg = rheome.ingest.config(ChannelMinTile=16);          % the diagonal, with the exemption on
            tc.file = rheome.ingest.build(tc.rec, tc.cfg, Verbose=false, Groups=P);
            tc.db = rheome.select.open(tc.file);
        end
    end

    methods (Test)

        function theTwoFloorsDifferAndTheStoreSaysSo(tc)
            tc.verifyGreaterThan(tc.db.channelLevel, 0);
            tc.verifyEqual(tc.db.envelopeLevel, 0);
            tc.verifyEqual(tc.db.grid.tExtent(tc.db.channelLevel + 1), 16);
        end

        function onlyMinAndMaxAreWrittenBelowTheChannelFloor(tc)
            w = whos('-file', tc.file);
            nm = string({w.name});
            for L = 0:tc.db.channelLevel-1
                tc.verifyTrue(ismember(sprintf('L%02d_min', L), nm), sprintf('L%02d_min', L));
                tc.verifyTrue(ismember(sprintf('L%02d_max', L), nm), sprintf('L%02d_max', L));
                for f = ["sumX","sumX2","absMax","energy","envMax"]
                    tc.verifyFalse(ismember(sprintf('L%02d_%s', L, f), nm), sprintf('L%02d_%s', L, f));
                end
            end
            L = tc.db.channelLevel;                               % at the floor, everything
            for f = ["sumX","sumX2","absMax","min","max","energy","envMax"]
                tc.verifyTrue(ismember(sprintf('L%02d_%s', L, f), nm), sprintf('L%02d_%s', L, f));
            end
        end

        function theExemptionCostsWhatItShould(tc)
            % Two extra single columns per channel per tile, and nothing else.
            g = tc.db.grid;
            extra = 0;
            for L = 0:tc.db.channelLevel-1, extra = extra + 2 * tc.C * g.K(L+1) * 4; end
            w = whos('-file', tc.file);
            got = sum([w(startsWith({w.name}, 'L') & (endsWith({w.name}, '_min') | endsWith({w.name}, '_max')) ...
                        & ~contains({w.name}, '_g_')).bytes]);
            tc.verifyGreaterThanOrEqual(got, extra);
            tc.verifyLessThan(got, 3 * extra);                    % the floor's own rows, not more
        end

        function derivingTheEnvelopeBelowTheFloorWorksAndTheRestRefuses(tc)
            R = rheome.select.derive(tc.db, Level=0, Channels=2, Stats=["min","max","ptp"]);
            tc.verifyEqual(height(R), tc.db.grid.K(1));
            tc.verifyEqual(R.ptp, R.max - R.min);
            tc.verifyError(@() rheome.select.derive(tc.db, Level=0, Channels=2, Stats="rms"), 'select:level');
            tc.verifyError(@() rheome.select.derive(tc.db, Level=0, Channels=2, Stats="share"), 'select:level');
            try
                rheome.select.derive(tc.db, Level=0, Channels=2, Stats="rms");
            catch err
                tc.verifySubstring(err.message, 'envelope');       % the message names the way out
            end
        end

        function theStoredPairContainsEverySampleOfItsTile(tc)
            % ⭐ The guarantee the whole envelope rests on, read back off disk.
            for L = [0 1 3 6]
                R = rheome.select.derive(tc.db, Level=L, Stats=["min","max"]);
                ext = tc.db.grid.tExtent(L+1);
                for c = 1:tc.C
                    r = R(R.unit_id == c, :);
                    for k = 1:height(r)
                        a = round(r.t_lo(k) * tc.fs) + 1;  b = min(round(r.t_hi(k) * tc.fs), tc.nT);
                        x = tc.X(a:b, c);
                        tc.verifyLessThanOrEqual(r.min(k), min(x), sprintf('L%d c%d k%d', L, c, k));
                        tc.verifyGreaterThanOrEqual(r.max(k), max(x), sprintf('L%d c%d k%d', L, c, k));
                    end
                end
            end
        end

        function turningTheExemptionOffRestoresTheOldLayout(tc)
            f = rheome.ingest.build(tc.rec, rheome.ingest.config(ChannelMinTile=16, ChannelEnvelope=false), ...
                             Verbose=false, Groups=[0 0 0; 0.03 0 0; 0 0.03 0; 0.03 0.03 0], Overwrite=true);
            d2 = rheome.select.open(f);
            tc.verifyEqual(d2.envelopeLevel, d2.channelLevel);
            w = whos('-file', f);
            tc.verifyFalse(ismember('L00_min', string({w.name})));
            tc.verifyError(@() rheome.select.derive(d2, Level=0, Channels=1, Stats="min"), 'select:level');
        end

    end
end

% Author: Diellor Basha, 2026
