classdef tIngestGroups < matlab.unittest.TestCase
% Group rows are exact roll-ups of channel rows along the sensor tree: sums add, maxima
% take the max, and rolling up in time then space equals space then time.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 2000
        C  = 12
        X
        P
        cfg
        fb
        bands
        frame
        g
        tree
    end

    methods (TestClassSetup)
        function make(tc)
            rng(41);
            t = (0:tc.nT-1)' / tc.fs;
            tc.X = 1e-13*randn(tc.nT, tc.C) + 3e-13*sin(2*pi*9*t) .* (1 + 0.5*sin(2*pi*0.3*t)) .* (1:tc.C) / tc.C;
            arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 2e-2);
            tc.P = arr.Pos;
            tc.tree = rheome.sensors.tree(rheome.sensors.graph(arr, 'Faces', false));
            tc.cfg = rheome.ingest.config(ChannelMinTile=0);
            [tc.fb, tc.bands, tc.frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
        end
    end

    methods (Test)

        function groupSumsEqualTheSumOfTheirChannels(tc)
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            Q = rheome.ingest.groups(T0, tc.tree, 1:tc.C);
            for i = 1:numel(Q.nodes)
                mem = tc.tree.members{Q.nodes(i)};
                tc.verifyEqual(Q.sumX(:, i),  sum(T0.sumX(:, mem), 2),  'RelTol', 1e-12);
                tc.verifyEqual(Q.sumX2(:, i), sum(T0.sumX2(:, mem), 2), 'RelTol', 1e-12);
                tc.verifyEqual(Q.energy(:, i, :), sum(T0.energy(:, mem, :), 2), 'RelTol', 1e-12);
                tc.verifyEqual(Q.absMax(:, i), max(T0.absMax(:, mem), [], 2));
                tc.verifyEqual(Q.max(:, i),    max(T0.max(:, mem), [], 2));
                tc.verifyEqual(Q.min(:, i),    min(T0.min(:, mem), [], 2));
                tc.verifyEqual(Q.envMax(:, i, :), max(T0.envMax(:, mem, :), [], 2));
            end
            tc.verifyEqual(Q.nodes(1), 1);                                       % the root first
            tc.verifyEqual(Q.sumX2(:, 1), sum(T0.sumX2, 2), 'RelTol', 1e-12);
        end

        function timeAndSpaceRollUpsCommute(tc)
            T0 = rheome.ingest.reduce(tc.X, tc.fb, tc.bands, tc.g, tc.frame);
            T  = rheome.ingest.rollup(T0, tc.g);
            A = rheome.ingest.groups(T, tc.tree, 1:tc.C);                               % space after time
            B = rheome.ingest.rollup(rheome.ingest.groups(T0, tc.tree, 1:tc.C), tc.g);       % time after space
            for L = 1:numel(T)
                tc.verifyEqual(A{L}.energy, B{L}.energy, 'RelTol', 1e-12);
                tc.verifyEqual(A{L}.envMax, B{L}.envMax);
                tc.verifyEqual(A{L}.min, B{L}.min);
            end
        end

        function aChannelSubsetMapsThroughChanRows(tc)
            % the store holds channels 2, 5, 7 of the array; groups sum only those
            T0 = rheome.ingest.reduce(tc.X(:, [2 5 7]), tc.fb, tc.bands, tc.g, tc.frame);
            Q = rheome.ingest.groups(T0, tc.tree, [2 5 7]);
            tc.verifyEqual(Q.sumX2(:, 1), sum(T0.sumX2, 2), 'RelTol', 1e-12);
            leafParent = tc.tree.parent_id(tc.tree.channel_id == 5);
            i = find(Q.nodes == leafParent);
            mem = intersect(tc.tree.members{leafParent}, [2 5 7]);
            [~, loc] = ismember(mem, [2 5 7]);
            tc.verifyEqual(Q.sumX2(:, i), sum(T0.sumX2(:, loc), 2), 'RelTol', 1e-12);
        end

        function theBuildWritesGroupArraysWhenPositionsAreGiven(tc)
            d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            rec = fullfile(d, 'recording.mat');
            F = tc.X.';  Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;                     %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);  ChannelName = compose('S%d', 1:tc.C);  ChannelType = repmat({'MEG'}, 1, tc.C); %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;  Comment = '';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
                 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
            f = rheome.ingest.build(rec, tc.cfg, Verbose=false, Groups=tc.P);
            m = matfile(f);  meta = m.meta;  g = m.grid;
            tc.verifyEqual(height(meta.tree), 2 * tc.C - 1);
            w = whos(m);  names = {w.name};
            tc.verifyTrue(ismember('L00_g_energy', names));
            tc.verifyEqual(size(m.L00_g_sumX2), [g.K0, tc.C - 1]);
            tc.verifyEqual(m.L00_g_sumX2(:, 1), sum(m.L00_sumX2, 2), 'RelTol', 1e-12);    % the root = all channels
            fNo = rheome.ingest.build(rec, tc.cfg, Verbose=false, Groups=false, Overwrite=true);
            wn = whos(matfile(fNo));
            tc.verifyFalse(ismember('L00_g_energy', {wn.name}));
        end

    end
end
