classdef tIngestBuild < matlab.unittest.TestCase
% The store: layout, sliceability, channel selection, the hash guard, determinism, and
% the long-form rows. A synthetic recording store in the +import/recording.m layout, so
% no dataset is needed. 20 s at 100 Hz, four channels of which one is a bad MEG and one
% is not MEG.
%
% Author: Diellor Basha, 2026

    properties
        Dir
        RecFile
        fs = 100
        nT = 2000
        F
        cfg
    end

    methods (TestClassSetup)
        function writeStore(tc)
            tc.Dir = tempname;  mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            tc.RecFile = fullfile(tc.Dir, 'recording.mat');
            rng(5);
            t = (0:tc.nT-1) / tc.fs;
            F = 1e-13*randn(4, tc.nT) + 3e-13*sin(2*pi*9*t);       %#ok<PROP>
            tc.F = F;
            Time = t;  sfreq = tc.fs;                               %#ok<NASGU>
            ChannelFlag = [1; 1; -1; 1];                            %#ok<NASGU>
            ChannelName = {'MLC11', 'MLC12', 'MLC13', 'EOG1'};      %#ok<NASGU>
            ChannelType = {'MEG', 'MEG', 'MEG', 'EOG'};             %#ok<NASGU>
            nCh = 4;  nT = tc.nT;  Comment = 'synthetic';           %#ok<NASGU>
            Events = [];  nAvg = 1;  source = struct('test', true); %#ok<NASGU>
            save(tc.RecFile, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', ...
                 'ChannelType', 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', ...
                 '-v7.3', '-nocompression');
            tc.cfg = rheome.ingest.config();
        end
    end

    methods (Test)

        function theStoreIsWrittenBesideTheRecordingWithParametersInItsName(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            tc.verifyEqual(fileparts(file), tc.Dir);
            [~, nm] = fileparts(file);
            tc.verifyEqual(nm, 'ingest__F0.25__V4__frame');                  % the default bank, four voices
            tc.verifyEqual(exist(file, 'file'), 2);
        end

        function onlyGoodMegChannelsAreSelected(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            m = matfile(file);  meta = m.meta;
            tc.verifyEqual(meta.iChannel, [1 2]);
            tc.verifyEqual(meta.ChannelName, {'MLC11', 'MLC12'});
            tc.verifyEqual(meta.C, 2);
            tc.verifyEqual(size(m.L00_sumX, 2), 2);
        end

        function everyLevelIsAFlatVariableWithTheDiagonalShape(tc)
            % Level L carries the bands whose natural level is <= L (grid.bandsAt).
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            m = matfile(file);  g = m.grid;  bands = m.bands;
            w = whos(m);
            names = {w.name};
            for L = 0:g.Lmax
                K = g.K(L+1);  nb = numel(g.bandsAt{L+1});
                for st = ["n","sumX","sumX2","absMax","min","max","energy","envMax","nCoi"]
                    v = sprintf('L%02d_%s', L, st);
                    tc.verifyTrue(ismember(v, names), v);
                end
                tc.verifyEqual(g.bandsAt{L+1}, find(bands.naturalLevel <= L)');
                tc.verifyEqual(size(m.(sprintf('L%02d_energy', L)), 3), nb);
                tc.verifyEqual(size(m.(sprintf('L%02d_nCoi', L)), 2),   nb);
                tc.verifyEqual(size(m.(sprintf('L%02d_n', L))),         [K 1]);
            end
            tc.verifyLessThan(numel(g.bandsAt{1}), height(bands));          % level 0 is not all bands
            tc.verifyEqual(numel(g.bandsAt{end}), height(bands));           % the top is
        end

        function aSliceThroughMatfileEqualsTheDirectReduce(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            m = matfile(file);  g = m.grid;
            [fb, bands, frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            gg = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
            T0 = rheome.ingest.reduce(tc.F([1 2], :).', fb, bands, gg, frame);
            T  = rheome.ingest.rollup(T0, gg);
            % the first level that carries at least two bands (on this 100 Hz store even
            % the top band's support exceeds the 0.25 s floor, so level 0 may carry none)
            L = find(cellfun(@numel, g.bandsAt) >= 2, 1) - 1;
            bb = g.bandsAt{L+1};
            v = sprintf('L%02d_energy', L);
            tc.verifyEqual(m.(v)(:, 2, 2), T{L+1}.energy(:, 2, bb(2)));
            tc.verifyEqual(m.L00_sumX2(:, :), T0.sumX2);
            v = sprintf('L%02d_envMax', L);
            tc.verifyEqual(m.(v)(:, :, :), T{L+1}.envMax(:, :, bb));
            v = sprintf('L%02d_nCoi', L);
            tc.verifyEqual(m.(v)(:, :), T{L+1}.nCoi(:, bb));
        end

        function metaRecordsWhatARebuildNeeds(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            meta = load(file, 'meta');  meta = meta.meta;
            tc.verifyEqual(meta.hash, tc.cfg.hash);
            tc.verifyEqual(meta.fs, tc.fs);
            tc.verifyEqual(meta.nT, tc.nT);
            tc.verifyEqual(meta.nValid, tc.nT);
            tc.verifyEqual(meta.units, 'T');
            tc.verifyNotEmpty(meta.release);
            tc.verifyGreaterThan(meta.buildSeconds, 0);
        end

        function theSameHashIsNotRebuiltAndADifferentOneRefuses(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            d1 = dir(file);
            file2 = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            d2 = dir(file);
            tc.verifyEqual(file2, file);
            tc.verifyEqual(d1.datenum, d2.datenum);
            other = rheome.ingest.config(FrequencyLimits=[1 20]);       % same file name, other hash
            tc.verifyError(@() rheome.ingest.build(tc.RecFile, other, Verbose=false), 'ingest:build:hash');
            f3 = rheome.ingest.build(tc.RecFile, other, Verbose=false, Overwrite=true);
            meta = load(f3, 'meta');
            tc.verifyEqual(meta.meta.hash, other.hash);
            rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true);   % restore
        end

        function twoBuildsAreIdenticalInEveryLevelVariable(tc)
            f1 = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true);
            A = load(f1);
            f2 = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true, ChannelChunk=1);
            B = load(f2);
            names = fieldnames(A);
            for k = 1:numel(names)
                if startsWith(names{k}, 'L')
                    tc.verifyEqual(A.(names{k}), B.(names{k}), names{k});
                end
            end
        end

        function aMaskIsAppliedAndRecorded(tc)
            mask = true(tc.nT, 1);  mask(101:200) = false;
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true, Mask=mask);
            m = matfile(file);  meta = m.meta;
            tc.verifyTrue(meta.masked);
            tc.verifyEqual(meta.nValid, tc.nT - 100);
            n = m.L00_n;
            tc.verifyEqual(sum(double(n)), tc.nT - 100);
            rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true);   % restore
        end

        function explicitChannelIndicesAreHonoured(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true, Channels=[4 1]);
            meta = load(file, 'meta');  meta = meta.meta;
            tc.verifyEqual(meta.iChannel, [4 1]);
            tc.verifyEqual(meta.ChannelType, {'EOG', 'MEG'});
            rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false, Overwrite=true);   % restore
        end

        function longFormRowsMatchTheStoreAndSkipBandsBelowTheirNaturalLevel(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            m = matfile(file);  g = m.grid;  bands = m.bands;
            R = rheome.ingest.totable(file, Levels=2, Stats=["sumX2","energy","n"]);
            K = g.K(3);  present = g.bandsAt{3};  nb = numel(present);
            tc.verifyEqual(height(R), K*2 + K*2*nb + K);          % sumX2, energy, n
            tc.verifyEqual(unique(R.level), 2);
            tc.verifyEqual(unique(R.tExtent), 1);                 % 0.25 * 2^2
            tc.verifyEqual(unique(R.band(R.stat == "energy"))', present);
            b = present(end);  k = find(present == b);
            e = R(R.stat == "energy" & R.channel == 2 & R.band == b, :);
            tc.verifyEqual(e.value, m.L02_energy(:, 2, k));
            tc.verifyEqual(e.fCenter, repmat(bands.fCenter(b), K, 1));
            tc.verifyEqual(e.tCenter, g.tCenter{3}(:));
            ap = R(R.stat == "sumX2", :);
            tc.verifyTrue(all(ap.band == 0) && all(isnan(ap.fCenter)));
            nn = R(R.stat == "n", :);
            tc.verifyTrue(all(nn.channel == 0));
        end

        function totableDefaultsToTheTopLevel(tc)
            file = rheome.ingest.build(tc.RecFile, tc.cfg, Verbose=false);
            m = matfile(file);  g = m.grid;
            R = rheome.ingest.totable(file, Stats="sumX");
            tc.verifyEqual(height(R), 2);
            tc.verifyEqual(unique(R.level), g.Lmax);
        end

    end
end
