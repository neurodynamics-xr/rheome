classdef tIngestSpatial < matlab.unittest.TestCase
% The wavelength axis in the ingest: at the root of the sensor tree the spatial bands
% partition sumX2 per tile exactly (Parseval on the graph); the fields roll up in time
% and space and commute; the build writes them when positions exist.
%
% Author: Diellor Basha, 2026

    properties
        fs = 200
        nT = 2000
        C  = 12
        X
        arr
        G
        tree
        cfg
        g
    end

    methods (TestClassSetup)
        function make(tc)
            rng(43);
            t = (0:tc.nT-1)' / tc.fs;
            tc.arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 2e-2);
            tc.G = rheome.sensors.graph(tc.arr, 'Faces', false);
            tc.tree = rheome.sensors.tree(tc.G);
            pattern = 1 + 0.3 * tc.arr.Pos(:, 1)' / tc.arr.Aperture;           % a nearly uniform pattern: longest wavelengths
            tc.X = 1e-13*randn(tc.nT, tc.C) + 3e-13*sin(2*pi*9*t) .* pattern;
            tc.cfg = rheome.ingest.config(ChannelMinTile=0);
            tc.g = rheome.ingest.grid(tc.nT, tc.fs, tc.cfg);
        end
    end

    methods (Test)

        function theBandTableCoversTheSpectrumTightly(tc)
            [~, sbands, sframe] = rheome.ingest.spatial(@(a, b) tc.X(a:b, :), tc.G, tc.g, tc.cfg);
            tc.verifyEqual(sframe.A, 1, 'AbsTol', 1e-10);  tc.verifyEqual(sframe.B, 1, 'AbsTol', 1e-10);
            tc.verifyEqual(sbands.kind{1}, 'longer');  tc.verifyEqual(sbands.kind{end}, 'shorter');
            tc.verifyTrue(all(strcmp(sbands.kind(2:end-1), 'octave')));
            tc.verifyEqual(sbands.k_lo(2:end), sbands.k_hi(1:end-1), 'RelTol', 1e-12);     % contiguous in k
            tc.verifyEqual(sum(sbands.n_modes), tc.C, 'RelTol', 1e-10);                     % gains^2 sum to 1 per mode
            tc.verifyEqual(sbands.k_hi(end), sqrt(max(sframe.Lambda)), 'RelTol', 1e-12);
        end

        function theRootPartitionIsExactAndChannelsAreNot(tc)
            [S0, sbands] = rheome.ingest.spatial(@(a, b) tc.X(a:b, :), tc.G, tc.g, tc.cfg);
            [fb, bands, frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            T0 = rheome.ingest.reduce(tc.X, fb, bands, tc.g, frame);
            rootE = squeeze(sum(S0.senergy, [2 3]));                                       % over channels and spatial bands
            tc.verifyEqual(rootE, sum(T0.sumX2, 2), 'RelTol', 1e-10);
            perChan = sum(S0.senergy, 3);
            tc.verifyGreaterThan(max(abs(perChan - T0.sumX2) ./ T0.sumX2, [], 'all'), 1e-3);  % not per channel
            tc.verifyEqual(height(sbands), size(S0.senergy, 3));
        end

        function aNearlyUniformPatternSitsInTheLongestBand(tc)
            [S0, sbands] = rheome.ingest.spatial(@(a, b) tc.X(a:b, :), tc.G, tc.g, tc.cfg);
            E = squeeze(sum(S0.senergy, [1 2]));
            longs = strcmp(sbands.kind, 'longer');                                        % the spatial mean and beyond the array
            tc.verifyGreaterThan(sum(E(longs)) / sum(E), 0.8);
        end

        function timeAndSpaceRollUpsCommuteForTheSpatialFields(tc)
            [S0] = rheome.ingest.spatial(@(a, b) tc.X(a:b, :), tc.G, tc.g, tc.cfg);
            [fb, bands, frame] = rheome.ingest.bank(tc.nT, tc.fs, tc.cfg);
            T0 = rheome.ingest.reduce(tc.X, fb, bands, tc.g, frame);
            T0.senergy = S0.senergy;  T0.senvMax = S0.senvMax;
            T = rheome.ingest.rollup(T0, tc.g);
            A = rheome.ingest.groups(T, tc.tree, 1:tc.C);
            B = rheome.ingest.rollup(rheome.ingest.groups(T0, tc.tree, 1:tc.C), tc.g);
            for L = 1:numel(T)
                tc.verifyEqual(A{L}.senergy, B{L}.senergy, 'RelTol', 1e-12);
                tc.verifyEqual(A{L}.senvMax, B{L}.senvMax);
            end
            tc.verifyEqual(squeeze(sum(A{end}.senergy(1, 1, :))), sum(tc.X.^2, 'all'), 'RelTol', 1e-10);   % root, top level
        end

        function aMaskRemovesSamplesFromTheSpatialEnergies(tc)
            mask = true(tc.nT, 1);  mask(101:150) = false;                                  % frame 3
            S0 = rheome.ingest.spatial(@(a, b) tc.X(a:b, :), tc.G, tc.g, tc.cfg, Mask=mask);
            tc.verifyEqual(squeeze(S0.senergy(3, :, :)), zeros(tc.C, size(S0.senergy, 3)));
            tc.verifyEqual(squeeze(S0.senvMax(3, :, :)), zeros(tc.C, size(S0.senergy, 3), 'single'));
        end

        function theBuildWritesTheSpatialArrays(tc)
            d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            rec = fullfile(d, 'recording.mat');
            F = tc.X.';  Time = (0:tc.nT-1)/tc.fs;  sfreq = tc.fs;                     %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);  ChannelName = compose('S%d', 1:tc.C);  ChannelType = repmat({'MEG'}, 1, tc.C); %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;  Comment = '';  Events = [];  nAvg = 1;  source = struct(); %#ok<NASGU>
            save(rec, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', 'ChannelType', ...
                 'nCh', 'nT', 'Comment', 'Events', 'nAvg', 'source', '-v7.3', '-nocompression');
            f = rheome.ingest.build(rec, tc.cfg, Verbose=false, Groups=tc.arr.Pos);
            m = matfile(f);  meta = m.meta;  g = m.grid;  names = {whos(m).name};
            tc.verifyTrue(all(ismember({'L00_s_energy','L00_s_envMax','L00_gs_energy','L00_gs_envMax'}, names)));
            nS = height(meta.sbands);
            tc.verifyEqual(size(m.L00_s_energy, 1), g.K0);  tc.verifyEqual(size(m.L00_s_energy, 2), tc.C);
            tc.verifyEqual(size(m.L00_s_energy, 3), numel(g.chanSpace));                    % a sensor carries the shortest band only
            sl = g.spaceSlots;
            tc.verifyEqual(size(m.L00_gs_energy), [g.K0, height(sl)]);
            E = m.L00_gs_energy;  s2 = m.L00_sumX2;
            tc.verifyEqual(sum(E(:, sl.group_id == 1), 2), sum(s2, 2), 'RelTol', 1e-10);      % root partition in the store
            f2 = rheome.ingest.build(rec, rheome.ingest.config(Space=false, ChannelMinTile=0), Verbose=false, Groups=tc.arr.Pos, Overwrite=true);
            tc.verifyFalse(ismember('L00_s_energy', {whos(matfile(f2)).name}));
        end

    end
end
