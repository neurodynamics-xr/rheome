classdef tOutPath < matlab.unittest.TestCase
% rheome.load.outroot and rheome.load.outpath: where analysis outputs go.
%
% What is asserted: in a git clone the default root is <repo>/results (an installed .mltbx has no
% .git and defaults under userpath instead, checked by build/check_install.m); RHEOME_OUT overrides it; the layout is
% <root>/[<subject>/]<analysis>/<file> with "/" in the analysis making nested folders; the folder is
% created; and an empty file returns the folder itself.
%
% Author: Diellor Basha, 2026

    properties
        saved
    end

    methods (TestMethodSetup)
        function keepEnv(tc)
            tc.saved = getenv('RHEOME_OUT');
            tc.addTeardown(@() setenv('RHEOME_OUT', tc.saved));
        end
    end

    methods (Test)

        function defaultRootIsResultsInTheClone(tc)
            setenv('RHEOME_OUT', '');
            repo = fileparts(fileparts(fileparts(which('rheome.load.root'))));
            tc.verifyEqual(rheome.load.outroot(), fullfile(repo, 'results'));
            tc.verifyNotEqual(rheome.load.outroot(), rheome.load.root(), 'outputs must not share +data');
        end

        function envOverridesAndLayoutIsSubjectAnalysisFile(tc)
            d = tempname;  setenv('RHEOME_OUT', d);  tc.addTeardown(@() rmdir(d, 's'));
            p = rheome.load.outpath("alpha_occupancy", "x.mat", "sub01");
            tc.verifyEqual(p, fullfile(d, 'sub01', 'alpha_occupancy', 'x.mat'));
            tc.verifyTrue(isfolder(fileparts(p)));
            tc.verifyEqual(rheome.load.outpath("reports", "a"), fullfile(d, 'reports', 'a'));
            tc.verifyEqual(rheome.load.outpath("reports/figures", ""), fullfile(d, 'reports', 'figures'));
            tc.verifyTrue(isfolder(fullfile(d, 'reports', 'figures')));
            tc.verifyEqual(rheome.load.outpath("plant_track", "p.csv", "planted"), fullfile(d, 'planted', 'plant_track', 'p.csv'));
        end
    end
end

% Author: Diellor Basha, 2026
