classdef tPagedRecording < matlab.unittest.TestCase
% Page arithmetic for @pagedrecording. The core/overlap/clip bookkeeping is the bug
% class that silently corrupts a downstream table -- an off-by-one in the core mask
% double-counts or drops samples without ever raising an error -- so it is pinned here
% against a synthetic store whose contents are known exactly.
%
% Author: Diellor Basha, 2026

    properties
        File            % synthetic recording store
        C   = 6         % channels
        nT  = 100       % samples
        fs  = 10        % Hz
    end

    methods (TestClassSetup)
        function writeStore(tc)
            d = tempname;  mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            tc.File = fullfile(d, 'recording.mat');

            % F(c,n) = c*1000 + n  -- every sample identifies its own (channel, index)
            F = (1:tc.C)' * 1000 + (1:tc.nT);          %#ok<NASGU>
            Time        = (0:tc.nT-1) / tc.fs;          %#ok<NASGU>
            sfreq       = tc.fs;                        %#ok<NASGU>
            ChannelFlag = ones(tc.C, 1);                %#ok<NASGU>
            ChannelName = compose('CH%d', 1:tc.C);      %#ok<NASGU>
            ChannelType = repmat({'MEG'}, 1, tc.C);     %#ok<NASGU>
            nCh = tc.C;  nT = tc.nT;                    %#ok<NASGU>
            Comment = 'synthetic';  source = 'test';    %#ok<NASGU>
            save(tc.File, 'F', 'Time', 'sfreq', 'ChannelFlag', 'ChannelName', ...
                 'ChannelType', 'nCh', 'nT', 'Comment', 'source', '-v7.3', '-nocompression');
        end
    end

    methods (Test)

        function pagesTileTheRecordExactlyOnce(tc)
            % The cores must PARTITION 1:nT -- no gap (samples never measured), no
            % overlap (samples measured twice, inflating every statistic).
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            tc.verifyEqual(pr.NumPages, 4);
            seen = [];
            for k = 1:pr.NumPages
                p = page(pr, k);
                seen = [seen, p.samples(p.core)]; %#ok<AGROW>
            end
            tc.verifyEqual(sort(seen), 1:tc.nT);
        end

        function interiorPageCarriesFullOverlapOnBothSides(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            p  = page(pr, 2);                       % core 31:60, span 24:67
            tc.verifyEqual(p.samples([1 end]), [24 67]);
            tc.verifyEqual(p.samples(p.core), 31:60);
            tc.verifyEqual(p.clipped, [0 0]);
        end

        function firstAndLastPagesClipToTheRecordAndSaySo(tc)
            % An edge page cannot have its full COI margin. Silently returning a short
            % span would make the caller trust cone-contaminated samples.
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            p1 = page(pr, 1);
            tc.verifyEqual(p1.samples([1 end]), [1 37]);
            tc.verifyEqual(p1.clipped, [7 0]);
            p4 = page(pr, 4);                       % core 91:100 (partial), span 84:100
            tc.verifyEqual(p4.samples([1 end]), [84 100]);
            tc.verifyEqual(p4.samples(p4.core), 91:100);
            tc.verifyEqual(p4.clipped, [0 7]);
        end

        function pageDataIsTheStoredSlice(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            p  = page(pr, 2);
            expected = (1:tc.C)' * 1000 + (24:67);
            tc.verifyEqual(double(p.F), expected, 'AbsTol', 1e-9);
        end

        function coreMaskSelectsTheCoreSamplesOfTheData(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            p  = page(pr, 3);
            expected = (1:tc.C)' * 1000 + (61:90);
            tc.verifyEqual(double(p.F(:, p.core)), expected, 'AbsTol', 1e-9);
        end

        function channelSelectionSubsetsRowsAndCounts(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Channels', [2 5]);
            tc.verifyEqual(pr.NumChannels, 2);
            p = page(pr, 1);
            tc.verifyEqual(size(p.F, 1), 2);
            tc.verifyEqual(double(p.F(:, 1)), [2001; 5001], 'AbsTol', 1e-9);
        end

        function readsAsSingleByDefaultAndDoubleOnRequest(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30);
            tc.verifyClass(getfield(page(pr, 1), 'F'), 'single');
            prd = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Precision', 'double');
            tc.verifyClass(getfield(page(prd, 1), 'F'), 'double');
        end

        function timeVectorTracksTheSpan(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30, 'Overlap', 7);
            p  = page(pr, 2);
            tc.verifyEqual(p.t, (23:66) / tc.fs, 'AbsTol', 1e-12);   % sample n -> (n-1)/fs
        end

        function zeroOverlapIsTheDefault(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30);
            p  = page(pr, 2);
            tc.verifyEqual(pr.Overlap, 0);
            tc.verifyEqual(p.samples([1 end]), [31 60]);
            tc.verifyTrue(all(p.core));
        end

        function pageIndexOutOfRangeErrors(tc)
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30);
            tc.verifyError(@() page(pr, 0), 'pagedrecording:pageIndex');
            tc.verifyError(@() page(pr, 5), 'pagedrecording:pageIndex');
        end

        function fullAxisIsExposedWithoutReadingData(tc)
            % The caller declares the FULL time axis from the object; only page() touches F.
            pr = rheome.pagedrecording(tc.File, 'PageLength', 30);
            tc.verifyEqual(pr.NumSamples, tc.nT);
            tc.verifyEqual(pr.SamplingFrequency, tc.fs);
            tc.verifyEqual(numel(pr.Time), tc.nT);
            tc.verifyEqual(pr.Duration, tc.nT / tc.fs, 'AbsTol', 1e-12);
        end

    end
end

% Author: Diellor Basha, 2026
