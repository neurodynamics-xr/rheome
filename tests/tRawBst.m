classdef tRawBst < matlab.unittest.TestCase
% rheome.io.read.rawbst -- the native BST-BIN reader.
%
% Every constant here was established by MEASUREMENT against a real Brainstorm file
% (a reference subject) and against brainstorm3's own in_fread/in_fread_bst:
%
%   * data is float32, little-endian, ALREADY in Tesla and ALREADY CTF-compensated
%   * it is EPOCH-BLOCKED in header.epochsize blocks, CHANNEL-MAJOR inside each block --
%     and that holds even though sFile.epochs is EMPTY (a continuous file). Reading it as
%     one flat channel-major run gives plausible per-channel magnitudes and r ~ 0.
%   * the SSP/ICA Projector must be applied on read. Without it the reconstruction reaches
%     only r = 0.80 against Brainstorm's own import; with it, r = 1.0000 (min 0.9999).
%
% A reader that is WRONG here still returns numbers that look like MEG, which is why these
% are pinned against a synthetic file whose contents are known exactly.
%
% Author: Diellor Basha, 2026

    properties
        Dir; RawMat; BstFile; ChanFile
        nCh = 5; ep = 10; nS = 40; hdr = 64; fs = 100
        Ftrue; Proj; Proj2
    end

    methods (TestClassSetup)
        function writeSynthetic(tc)
            tc.Dir = tempname;  mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
            tc.BstFile  = fullfile(tc.Dir, 'synth.bst');
            tc.RawMat   = fullfile(tc.Dir, 'data_0raw_synth.mat');
            tc.ChanFile = fullfile(tc.Dir, 'channel_synth.mat');

            % F(c,n) = c + n/1000 -- every sample identifies its (channel, index)
            tc.Ftrue = (1:tc.nCh)' + (1:tc.nS)/1000;

            % write EPOCH-BLOCKED, CHANNEL-MAJOR within each block
            fid = fopen(tc.BstFile, 'w', 'ieee-le');
            fwrite(fid, zeros(tc.hdr,1), 'uint8');
            for e = 0:(tc.nS/tc.ep - 1)
                blk = tc.Ftrue(:, e*tc.ep + (1:tc.ep));    % [nCh x ep]
                fwrite(fid, reshape(blk.', [], 1), 'float32');   % channel-major inside
            end
            fclose(fid);

            F = struct('filename', tc.BstFile, 'format', 'BST-BIN', 'byteorder', 'l', ...
                'prop', struct('times', [0 (tc.nS-1)/tc.fs], 'sfreq', tc.fs, 'nAvg', 1, ...
                               'currCtfComp', 3, 'destCtfComp', 3), ...
                'epochs', [], 'channelflag', ones(tc.nCh,1), ...
                'header', struct('nchannels', tc.nCh, 'nsamples', tc.nS, ...
                                 'epochsize', tc.ep, 'hdrsize', tc.hdr, 'sfreq', tc.fs));
            ChannelFlag = ones(tc.nCh,1);  Time = (0:tc.nS-1)/tc.fs;  Comment = 'synth'; %#ok<NASGU>
            save(tc.RawMat, 'F', 'ChannelFlag', 'Time', 'Comment');

            % The DISABLED projector must be a DIFFERENT direction from the active one:
            % a projector is idempotent, so re-applying the same one changes nothing and a
            % test built on it can never fail. This one removes a direction the active one
            % leaves alone, so skipping it is observable.
            u  = zeros(tc.nCh,1); u(1)  = 0.6; u(2) = 0.8;    % active
            u2 = zeros(tc.nCh,1); u2(3) = 1.0;                % disabled -- distinct direction
            Projector = struct('Comment', {'active','disabled'}, ...
                'Components', {[u, zeros(tc.nCh, tc.nCh-1)], [u2, zeros(tc.nCh, tc.nCh-1)]}, ...
                'CompMask',   {[1 zeros(1,tc.nCh-1)], [1 zeros(1,tc.nCh-1)]}, ...
                'Status',     {1, 0}, ...
                'SingVal',    {ones(1,tc.nCh), ones(1,tc.nCh)}, 'Method', {'SSP_pca','SSP_pca'});
            tc.Proj  = eye(tc.nCh) - u*u';
            tc.Proj2 = eye(tc.nCh) - u2*u2';
            Channel = struct('Name', {'A','B','C','D','E'}, ...
                             'Type', {'MEG','MEG','MEG','MEG','MEG REF'}, 'Loc', {[],[],[],[],[]});
            MegRefCoef = [];                                    %#ok<NASGU>
            save(tc.ChanFile, 'Channel', 'Projector', 'MegRefCoef');
        end
    end

    methods (Test)

        function headerReadTouchesNoData(tc)
            R = rheome.io.read.rawbst(tc.RawMat);
            tc.verifyEqual(R.nCh, tc.nCh);
            tc.verifyEqual(R.nT, tc.nS);
            tc.verifyEqual(R.sfreq, tc.fs);
            tc.verifyEmpty(R.F);
        end

        function readsARangeInsideOneEpoch(tc)
            R = rheome.io.read.rawbst(tc.RawMat, [3 7], struct('projector', false));
            tc.verifyEqual(R.F, tc.Ftrue(:, 3:7), 'AbsTol', 1e-5);
        end

        function readsARangeSpanningEpochBoundaries(tc)
            % The block structure is the whole trap: a range crossing a boundary is where a
            % flat-layout reader silently returns the wrong samples.
            R = rheome.io.read.rawbst(tc.RawMat, [8 23], struct('projector', false));
            tc.verifyEqual(R.F, tc.Ftrue(:, 8:23), 'AbsTol', 1e-5);
        end

        function readsTheWholeRecord(tc)
            R = rheome.io.read.rawbst(tc.RawMat, [1 tc.nS], struct('projector', false));
            tc.verifyEqual(R.F, tc.Ftrue, 'AbsTol', 1e-5);
        end

        function appliesActiveProjectorsByDefault(tc)
            % Measured on real data: without this the reconstruction tops out at r = 0.80.
            R = rheome.io.read.rawbst(tc.RawMat, [3 7]);
            tc.verifyEqual(R.F, tc.Proj * tc.Ftrue(:, 3:7), 'AbsTol', 1e-5);
            tc.verifyEqual(R.nProjector, 1);
        end

        function skipsProjectorsThatAreNotActive(tc)
            % Two are stored; only Status==1 counts. Applying the disabled one would square
            % the correction and quietly attenuate real signal.
            R = rheome.io.read.rawbst(tc.RawMat, [3 7]);
            tc.verifyEqual(R.F, tc.Proj * tc.Ftrue(:, 3:7), 'AbsTol', 1e-5);
            % had the disabled one been applied too, channel 3 would have been zeroed
            bothApplied = tc.Proj2 * tc.Proj * tc.Ftrue(:, 3:7);
            tc.verifyGreaterThan(norm(R.F - bothApplied), 1);
            tc.verifyGreaterThan(norm(R.F(3,:)), 1);
        end

        function timeVectorMatchesTheRequestedSamples(tc)
            R = rheome.io.read.rawbst(tc.RawMat, [8 23]);
            tc.verifyEqual(R.Time, (7:22)/tc.fs, 'AbsTol', 1e-12);
            tc.verifyEqual(R.samples, [8 23]);
        end

        function rejectsAnOutOfRangeRequest(tc)
            tc.verifyError(@() rheome.io.read.rawbst(tc.RawMat, [0 5]),  'io:read:rawbst:range');
            tc.verifyError(@() rheome.io.read.rawbst(tc.RawMat, [1 999]),'io:read:rawbst:range');
        end

    end
end

% Author: Diellor Basha, 2026
