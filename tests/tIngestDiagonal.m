classdef tIngestDiagonal < matlab.unittest.TestCase
% The position diagonal: channel rows only from ChannelMinTile up, wavelength bands only
% on nodes at least one wavelength across (plus the root, and the shortest band on
% sensors); the queries say so below it, and everything above it is unchanged.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 4000
        C  = 12
        arr
        rec
        dir_
        X
    end

    methods (TestClassSetup)
        function make(tc)
            rng(47);
            tc.arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 2e-2);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, tc.C) + 3e-13*sin(2*pi*9*t) .* (1 + 0.2*tc.arr.Pos(:, 1)' / tc.arr.Aperture);
            tc.dir_ = tempname;  mkdir(tc.dir_);
            tc.addTeardown(@() rmdir(tc.dir_, 's'));
            tc.rec = fullfile(tc.dir_, 'recording.mat');
            F = tc.X.';  Time = t';  sfreq = tc.fs;                                    %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);  ChannelName = compose('S%d', 1:tc.C);  ChannelType = repmat({'MEG'}, 1, tc.C); %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;  Comment = '';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(tc.rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
                 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
        end
    end

    methods (Test)

        function channelRowsStartAtTheChosenTile(tc)
            f = rheome.ingest.build(tc.rec, rheome.ingest.config(ChannelMinTile=2), Verbose=false, Groups=tc.arr.Pos, Overwrite=true);   % level 3
            m = matfile(f);  g = m.grid;  names = {whos(m).name};
            tc.verifyEqual(g.channelLevel, 3);
            for L = 0:2, tc.verifyFalse(ismember(sprintf('L%02d_sumX2', L), names)); tc.verifyTrue(ismember(sprintf('L%02d_g_sumX2', L), names)); end
            for L = 3:g.Lmax, tc.verifyTrue(ismember(sprintf('L%02d_sumX2', L), names)); end
            tc.verifyTrue(ismember('L00_n', names));                                     % per-record rows stay
            db = rheome.select.open(f);
            q = rheome.select.query(Stat="absMax", Op=">", Threshold="p50", Level=1);
            tc.verifyError(@() rheome.select.frames(db, q), 'select:level');
            tc.verifyError(@() rheome.select.rows(db, "feature", Level=2), 'select:level');
            qg = rheome.select.query(Stat="absMax", Op=">", Threshold="p50", Level=1, Scope="group");
            [Rf, ~] = rheome.select.frames(db, qg);  Rs = rheome.select.scan(db, qg);                  % groups carry the fine levels
            tc.verifyEqual(height(Rf), height(Rs));  tc.verifyGreaterThan(height(Rf), 0);
            q3 = rheome.select.query(Stat="absMax", Op=">", Threshold="p50", Level=3);
            tc.verifyEqual(height(rheome.select.frames(db, q3)), height(rheome.select.scan(db, q3)));
        end

        function wavelengthBandsSitOnNodesThatCanHoldThem(tc)
            f = rheome.ingest.build(tc.rec, rheome.ingest.config(ChannelMinTile=0), Verbose=false, Groups=tc.arr.Pos, Overwrite=true);
            m = matfile(f);  g = m.grid;  meta = m.meta;  sl = g.spaceSlots;  sb = meta.sbands;  tr = meta.tree;
            tc.assumeTrue(all(isfinite(sb.wavelength_lo)), 'needs a calibrated grid');
            tc.verifyGreaterThan(height(sl), 0);
            for i = 1:height(sl)
                n = sl.group_id(i);  b = sl.sband_id(i);
                tc.verifyTrue(tr.parent_id(n) == 0 || tr.diameter(n) >= sb.wavelength_lo(b));
            end
            tc.verifyEqual(nnz(sl.group_id == 1), height(sb));                             % the root carries every band
            tc.verifyLessThan(height(sl), height(sb) * numel(meta.groupNodes));           % something was dropped
            tc.verifyEqual(g.chanSpace, find(sb.wavelength_lo <= 0)');                      % sensors: the shortest band only
            db = rheome.select.open(f);
            longest = find(strcmp(sb.kind, 'octave'), 1);
            small = tr.node_id(find(~tr.is_leaf & tr.parent_id ~= 0 & tr.diameter < sb.wavelength_lo(longest), 1));
            tc.assumeNotEmpty(small, 'needs a non-root node smaller than the longest octave');
            q = rheome.select.query(Stat="spaceEnergy", Op=">", Threshold=0, Level=2, Band=longest, Scope="group", Groups=small);
            tc.verifyError(@() rheome.select.scan(db, q), 'select:space');
            qr = rheome.select.query(Stat="spaceEnergy", Op=">", Threshold=0, Level=2, Band=longest, Scope="group", Groups=1);
            tc.verifyEqual(height(rheome.select.scan(db, qr)), g.K(3));
            R = rheome.select.rows(db, "feature_group_space", Level=2);
            tc.verifyEqual(height(R), g.K(3) * height(sl));
        end

        function theDiagonalCanBeSwitchedOff(tc)
            f = rheome.ingest.build(tc.rec, rheome.ingest.config(ChannelMinTile=0, SpaceDiagonal=false), Verbose=false, Groups=tc.arr.Pos, Overwrite=true);
            m = matfile(f);  g = m.grid;  meta = m.meta;
            tc.verifyEqual(height(g.spaceSlots), height(meta.sbands) * numel(meta.groupNodes));
            tc.verifyEqual(g.chanSpace, 1:height(meta.sbands));
        end

    end
end
