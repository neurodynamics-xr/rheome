classdef tImportRawRecording < matlab.unittest.TestCase
% rheome.import.rawrecording -- native-rate BST-BIN into the pageable v7.3 store.
%
% Author: Diellor Basha, 2026

    properties
        Name = '_tmp_rawrec_test'
        DsDir; RawMat
        nCh = 5; ep = 10; nS = 40; hdr = 64; fs = 100
        Ftrue; Proj
    end

    methods (TestClassSetup)
        function build(tc)
            src = tempname; mkdir(src);
            tc.addTeardown(@() rmdir(src, 's'));
            bstF = fullfile(src, 'synth.bst');
            tc.RawMat = fullfile(src, 'data_0raw_synth.mat');

            tc.Ftrue = (1:tc.nCh)' + (1:tc.nS)/1000;
            fid = fopen(bstF, 'w', 'ieee-le');
            fwrite(fid, zeros(tc.hdr,1), 'uint8');
            for e = 0:(tc.nS/tc.ep - 1)
                blk = tc.Ftrue(:, e*tc.ep + (1:tc.ep));
                fwrite(fid, reshape(blk.', [], 1), 'float32');
            end
            fclose(fid);

            F = struct('filename', bstF, 'format', 'BST-BIN', 'byteorder', 'l', ...
                'prop', struct('times', [0 (tc.nS-1)/tc.fs], 'sfreq', tc.fs, 'nAvg', 1, ...
                               'currCtfComp', 3, 'destCtfComp', 3), ...
                'epochs', [], 'channelflag', ones(tc.nCh,1), ...
                'header', struct('nchannels', tc.nCh, 'nsamples', tc.nS, ...
                                 'epochsize', tc.ep, 'hdrsize', tc.hdr, 'sfreq', tc.fs));
            ChannelFlag = ones(tc.nCh,1); Time = (0:tc.nS-1)/tc.fs; Comment = 'synth raw'; %#ok<NASGU>
            save(tc.RawMat, 'F', 'ChannelFlag', 'Time', 'Comment');

            u = zeros(tc.nCh,1); u(1) = 0.6; u(2) = 0.8;
            Projector = struct('Comment', {'c'}, 'Components', {[u, zeros(tc.nCh,tc.nCh-1)]}, ...
                'CompMask', {[1 zeros(1,tc.nCh-1)]}, 'Status', {1}, ...
                'SingVal', {ones(1,tc.nCh)}, 'Method', {'SSP_pca'});
            tc.Proj = eye(tc.nCh) - u*u';
            Channel = struct('Name', {'A','B','C','D','E'}, ...
                             'Type', {'MEG','MEG','MEG','MEG','MEG REF'}, 'Loc', {[],[],[],[],[]});
            MegRefCoef = []; %#ok<NASGU>
            save(fullfile(src, 'channel_synth.mat'), 'Channel', 'Projector', 'MegRefCoef');

            tc.DsDir = fullfile(rheome.load.root(), tc.Name);
            tc.assumeFalse(exist(tc.DsDir,'dir') == 7, 'scratch dataset exists');
            mkdir(tc.DsDir);
            tc.addTeardown(@() rmdir(tc.DsDir, 's'));
        end
    end

    methods (Test)

        function storesSingleBecauseTheSourceIsFloat32(tc)
            % The .bst is float32, so single is BIT-EXACT, not lossy -- and halves the store.
            f = rheome.import.rawrecording(tc.Name, tc.RawMat);
            w = whos('-file', f);
            iF = strcmp({w.name}, 'F');
            tc.verifyEqual(w(iF).class, 'single');
            tc.verifyEqual(w(iF).size, [tc.nCh tc.nS]);
        end

        function contentEqualsTheReaderOverTheWholeRecord(tc)
            f = rheome.import.rawrecording(tc.Name, tc.RawMat);
            R = rheome.io.read.rawbst(tc.RawMat, [1 tc.nS]);
            m = matfile(f);
            tc.verifyEqual(double(m.F(:, 1:tc.nS)), R.F, 'AbsTol', 1e-4);
            tc.verifyEqual(double(m.F(:, 1:tc.nS)), tc.Proj * tc.Ftrue, 'AbsTol', 1e-4);
        end

        function chunkBoundariesDoNotCorruptTheStore(tc)
            % Written in chunks, so a chunk seam is exactly where samples get dropped or
            % duplicated -- and the result still looks like plausible data.
            f = rheome.import.rawrecording(tc.Name, tc.RawMat, struct('chunkSamples', 7));
            R = rheome.io.read.rawbst(tc.RawMat, [1 tc.nS]);
            m = matfile(f);
            tc.verifyEqual(double(m.F(:, 1:tc.nS)), R.F, 'AbsTol', 1e-4);
        end

        function carriesTheNativeRateNotAResample(tc)
            f = rheome.import.rawrecording(tc.Name, tc.RawMat);
            h = builtin('load', f, 'sfreq', 'nT', 'nCh', 'Time', 'nProjector');
            tc.verifyEqual(h.sfreq, tc.fs, 'AbsTol', 1e-12);
            tc.verifyEqual(h.nT, tc.nS);
            tc.verifyEqual(h.nCh, tc.nCh);
            tc.verifyEqual(h.nProjector, 1);
        end

        function pagedrecordingOpensTheNativeStore(tc)
            rheome.import.rawrecording(tc.Name, tc.RawMat);
            pr = rheome.pagedrecording(tc.Name, 'PageLength', 10, 'Overlap', 3, 'Store', 'native');
            tc.verifyEqual(pr.NumSamples, tc.nS);
            tc.verifyEqual(pr.SamplingFrequency, tc.fs);
            p = page(pr, 2);
            expected = tc.Proj * tc.Ftrue;
            tc.verifyEqual(double(p.F(:, p.core)), expected(:, 11:20), 'AbsTol', 1e-4);
        end

    end
end

% Author: Diellor Basha, 2026
