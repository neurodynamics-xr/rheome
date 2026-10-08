% Author: Diellor Basha, 2026
classdef tQcEmptyRoom < matlab.unittest.TestCase
% TQCEMPTYROOM  Detector v2: which recordings count as empty-room (never labelled), and that
% the detector's rule step leaves an empty-room recording with no flags whatever its metrics.
%
% See also: rheome.qc.isemptyroom, rheome.qc.badchannels

methods (Test)
    function omegaNoiseRunsAreEmptyRoom(tc)
        tc.verifyTrue(rheome.qc.isemptyroom('sub-01_ses-04_task-noise_run-02_meg'));
        tc.verifyTrue(rheome.qc.isemptyroom('/x/sub-02/@rawsub-02_ses-03_task-noise_run-01_meg/data_0raw_sub-02_ses-03_task-noise_run-01_meg.mat'));
        tc.verifyTrue(rheome.qc.isemptyroom('/data/sub-01/ses-04/meg/sub-01_ses-04_task-noise_run-02_meg.ds'));
    end
    function bidsEmptyRoomConventionsCount(tc)
        tc.verifyTrue(rheome.qc.isemptyroom('sub-emptyroom_ses-20190101_task-noise_meg'));
        tc.verifyTrue(rheome.qc.isemptyroom('sub-01_task-emptyroom_meg'));
        tc.verifyTrue(rheome.qc.isemptyroom('SUB-01_TASK-NOISE_MEG'));
    end
    function subjectRecordingsAreNot(tc)
        names = {'sub-01_ses-01_task-rest_run-01_meg', 'sub-03_ses-01_task-restaftertask_run-02_meg', ...
            'sub-01_task-restaftersleep_meg', 'sub-01_task-noisetest_meg', '/noise/sub-01_task-rest_meg', ...
            'sub-04_ses-01_task-rest_run-04_meg'};
        tf = rheome.qc.isemptyroom(names);
        tc.verifySize(tf, size(names));
        tc.verifyFalse(any(tf));
    end
    function stringArrayKeepsShape(tc)
        tf = rheome.qc.isemptyroom(["sub-1_task-noise_meg"; "sub-1_task-rest_meg"]);
        tc.verifyEqual(tf, [true; false]);
    end
    function detectorSourceGatesEveryRule(tc)
        % v2 must apply the empty-room gate AFTER all rules, so no rule (FLAT, NOISY, JUMPY, PSD_*) escapes it
        src = fileread(which('rheome.qc.badchannels'));
        tc.verifyNotEmpty(regexp(src, 'qc\.isemptyroom', 'once'));
        tc.verifyNotEmpty(regexp(src, "P\.version\s*=\s*'v2'", 'once'));
        iRules = strfind(src, '[flag, why] = i_rules(M, P);');
        iGate = strfind(src, 'flag(:) = false;');
        tc.verifyNotEmpty(iRules); tc.verifyNotEmpty(iGate);
        tc.verifyGreaterThan(iGate(1), iRules(1));
    end
end
end
% Author: Diellor Basha, 2026
