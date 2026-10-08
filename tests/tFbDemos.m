classdef tFbDemos < matlab.unittest.TestCase
% All three demos run end to end and every self-check passes.

    methods (TestMethodTeardown)
        function shut(~), close all; end
    end

    methods (Test)

        function allThreeRunAndPass(tc)
            names = {'filterbank_cwt','filterbank_graph','filterbank_joint'};
            for i = 1:numel(names)
                close all;
                out = feval(['rheome.demos.' names{i}]);
                if isfield(out,'skipped') && out.skipped
                    continue;                       % toolbox absent: documented skip
                end
                tc.verifyTrue(out.ok, sprintf('%s: every self-check must pass', names{i}));
            end
        end

        function exportWritesFivePngsPerDemo(tc)
            d = fullfile(tempdir, 'fbdemo_png_test');
            if exist(d,'dir'), rmdir(d,'s'); end
            c = onCleanup(@() rmdir(d,'s'));
            close all;
            rheome.demos.filterbank_graph(d);
            f = dir(fullfile(d, '*.png'));
            tc.verifyEqual(numel(f), 5, 'a figure opened but never saved would show up here');
        end

        function graphAndJointNeedNoToolbox(tc)
            % The suite must not become toolbox-dependent by accident.
            for n = {'filterbank_graph','filterbank_joint'}
                src = fileread(which(['rheome.demos.' n{1}]));
                tc.verifyEmpty(regexp(src, 'cwtfilterbank', 'once'), n{1});
            end
        end

    end
end

% Author: Diellor Basha, 2026
