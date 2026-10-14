classdef tReadBold < matlab.unittest.TestCase
% rheome.io.read.bold on a hand-made bst-fmri overlay: the run's confounds are paired by tag, the time axis
% gives the TR, and a map on another surface or with the wrong vertex count is refused.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function readsRunsWithConfounds(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = i_run(d, 'task_rest_run_1', 'sub-1/tess_cortex_pial_low.mat', 5, 8);
            i_run(d, 'task_rest_run_2', 'sub-1/tess_cortex_pial_low.mat', 5, 6);
            B = rheome.io.read.bold({f, strrep(f, 'run_1', 'run_2')}, 5);
            tc.verifyEqual({B.name}, {'task_rest_run_1', 'task_rest_run_2'});
            tc.verifyEqual([B.tr], [2 2]);  tc.verifySize(B(1).x, [5 8]);  tc.verifyClass(B(1).x, 'single');
            tc.verifySize(B(2).confounds, [6 2]);  tc.verifyEqual(B(1).confoundNames, ["csf" "white_matter"]);
        end

        function refusesWrongSurfaceOrSize(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            f = i_run(d, 'task_rest', 'sub-1/tess_cortex_pial_high.mat', 5, 8);
            tc.verifyError(@() rheome.io.read.bold({f}), 'io:read:bold:surface');
            f = i_run(d, 'task_rest', 'sub-1/tess_cortex_pial_low.mat', 5, 8);
            tc.verifyError(@() rheome.io.read.bold({f}, 6), 'io:read:bold:vertices');
            delete(fullfile(d, 'matrix_confounds_fmri_task_rest.mat'));
            tc.verifyError(@() rheome.io.read.bold({f}), 'io:read:bold:confounds');
        end
    end
end

function f = i_run(d, tag, surf, nV, nT)
    ImageGridAmp = rand(nV, nT);  Time = (0:nT-1) * 2;  SurfaceFile = surf;  DisplayUnits = 'a.u.'; %#ok<NASGU>
    f = fullfile(d, ['results_surface_fmri_' tag '.mat']);
    save(f, 'ImageGridAmp', 'Time', 'SurfaceFile', 'DisplayUnits');
    Value = rand(2, nT);  Description = {'csf'; 'white_matter'}; %#ok<NASGU>
    save(fullfile(d, ['matrix_confounds_fmri_' tag '.mat']), 'Value', 'Description', 'Time');
end

% Author: Diellor Basha, 2026
