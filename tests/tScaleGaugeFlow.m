classdef tScaleGaugeFlow < matlab.unittest.TestCase
% rheome.scale.run's periodicflow writes gaugeflow.csv (the alpha flow per ico3 sphere patch in the group
% gauge, charts z and x) on a synthetic cached subject with a registration sphere (scaleSynthSubject),
% and rheome.scale.reducegaugeflow reads it back.
%
% Author: Diellor Basha, 2026

    properties
        Name = 'synth_gaugeflow'
        Root
    end

    methods (TestClassSetup)
        function cache(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', tc.Root);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(tc.Root, tc.Name));
        end
    end

    methods (Test)
        function periodicflowWritesTheGaugeFlowPerPatch(tc)
            out = fullfile(tc.Root, 'out');
            rheome.scale.run(tc.Name, Analyses="periodicflow", OutDir=out, FlowTiles=2, Cohort="c", Dataset="d");
            f = fullfile(out, tc.Name, 'gaugeflow.csv');
            tc.assertTrue(isfile(f), 'periodicflow must write gaugeflow.csv when the surface has a sphere');
            G = readtable(f, 'TextType', 'string');
            tc.verifyEqual(sort(unique(G.chart))', ["x" "z"]);
            tc.verifyEqual(height(G), 2 * 2 * 642, 'two charts x two bands x 642 ico3 patches');
            tc.verifyEqual(unique(G.subject), string(tc.Name));
            tc.verifyEqual(sum(G.has_pole(G.chart == "z" & G.band == G.band(1))), 2);
            Gr = rheome.scale.reducegaugeflow(out, fullfile(tc.Root, 'group'));
            tc.verifyTrue(all(Gr.n_subjects == 1));
            tc.verifyTrue(isfile(fullfile(tc.Root, 'group', 'group_gaugeflow.csv')));
        end
    end
end

% Author: Diellor Basha, 2026
