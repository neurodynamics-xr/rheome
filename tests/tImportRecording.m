classdef tImportRecording < matlab.unittest.TestCase
% rheome.import.recording surfaces the cached recording as a v7.3 store whose F is a TOP-LEVEL
% variable. The whole point is sliceability: a store that round-trips its numbers but
% still nests F inside a struct is useless, because MatFile cannot index into one. That
% is what these tests actually assert.
%
% Runs against a small synthetic dataset folder in the +data cache, removed afterwards.
%
% Author: Diellor Basha, 2026

    properties
        Name = '_tmp_import_recording_test'
        Dir
        C   = 5
        nT  = 40
        fs  = 20
    end

    methods (TestClassSetup)
        function writeSyntheticStudy(tc)
            tc.Dir = fullfile(rheome.load.root(), tc.Name);
            tc.assumeFalse(exist(tc.Dir, 'dir') == 7, 'test scratch dataset already exists');
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));

            rec = struct();
            rec.F           = (1:tc.C)' * 1000 + (1:tc.nT);
            rec.Time        = (0:tc.nT-1) / tc.fs;
            rec.ChannelFlag = [1 1 -1 1 1]';
            rec.sfreq       = tc.fs;
            rec.nCh         = tc.C;
            rec.nT          = tc.nT;
            rec.Comment     = 'synthetic raw block';
            rec.Events      = [];
            rec.nAvg        = 1;

            chan = struct('Name', {compose('CH%d', 1:tc.C)}, ...
                          'Type', {{'MEG','MEG','MEG','MEG REF','EOG'}});
            studyDir = 'synthetic';  dataName = 'synthetic.mat';
            hm = struct(); ncov = struct();                       %#ok<NASGU>
            save(fullfile(tc.Dir, 'study.mat'), 'rec', 'chan', 'hm', 'ncov', ...
                 'studyDir', 'dataName', '-v7.3');
        end
    end

    methods (Test)

        function fIsATopLevelVariableSoMatfileCanSliceIt(tc)
            % The defining requirement. Nested inside a struct this call errors with
            % "MatFile objects only support '()' indexing" -- which is exactly why the
            % store exists separately from study.mat.
            f = rheome.import.recording(tc.Name);
            m = matfile(f);
            tc.verifyTrue(ismember('F', who(m)), 'F must be a top-level variable');
            slice = m.F(:, 11:15);
            tc.verifyEqual(double(slice), (1:tc.C)' * 1000 + (11:15), 'AbsTol', 1e-9);
        end

        function headerRoundTripsWithoutTheData(tc)
            f = rheome.import.recording(tc.Name);
            h = builtin('load', f, 'Time', 'sfreq', 'ChannelFlag', 'nCh', 'nT');
            tc.verifyEqual(h.sfreq, tc.fs, 'AbsTol', 1e-12);
            tc.verifyEqual(h.nCh, tc.C);
            tc.verifyEqual(h.nT, tc.nT);
            tc.verifyEqual(h.Time, (0:tc.nT-1) / tc.fs, 'AbsTol', 1e-12);
            tc.verifyEqual(h.ChannelFlag(:), [1 1 -1 1 1]');
        end

        function channelNamesAndTypesComeFromTheChannelFile(tc)
            % The MEG/good selection downstream is built from these; without them the
            % caller has to reopen study.mat and the paging saves nothing.
            f = rheome.import.recording(tc.Name);
            h = builtin('load', f, 'ChannelName', 'ChannelType');
            tc.verifyEqual(h.ChannelType, {'MEG','MEG','MEG','MEG REF','EOG'});
            tc.verifyEqual(h.ChannelName, compose('CH%d', 1:tc.C));
        end

        function theStoreIsWhatPagedrecordingConsumes(tc)
            rheome.import.recording(tc.Name);
            pr = rheome.pagedrecording(tc.Name, 'PageLength', 10, 'Overlap', 3);
            tc.verifyEqual(pr.NumSamples, tc.nT);
            tc.verifyEqual(pr.NumPages, 4);
            p = page(pr, 2);
            tc.verifyEqual(double(p.F(:, p.core)), (1:tc.C)' * 1000 + (11:20), 'AbsTol', 1e-9);
        end

        function refusesARecordingThatCarriesNoData(tc)
            % rheome.io.read.recording leaves F empty for a raw LINK file. Writing a store with
            % no F would fail later, at page(), far from the cause.
            d = fullfile(rheome.load.root(), '_tmp_import_recording_empty');
            mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            rec = struct('F', [], 'Time', [], 'ChannelFlag', [], 'sfreq', [], ...
                         'nCh', 0, 'nT', 0, 'Comment', '', 'Events', [], 'nAvg', 1);
            chan = struct('Name', {{}}, 'Type', {{}});
            studyDir = ''; dataName = ''; hm = struct(); ncov = struct();  %#ok<NASGU>
            save(fullfile(d, 'study.mat'), 'rec', 'chan', 'hm', 'ncov', ...
                 'studyDir', 'dataName', '-v7.3');
            tc.verifyError(@() rheome.import.recording('_tmp_import_recording_empty'), ...
                'import:recording:noData');
        end

    end
end

% Author: Diellor Basha, 2026
