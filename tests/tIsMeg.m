classdef tIsMeg < matlab.unittest.TestCase
% rheome.sensors.ismeg: CTF and MEGIN measurement channels in, references and other modalities out; and
% rheome.scale.sensors on a MEGIN-typed subject keeps the types, so the MNE whitens mags and grads apart.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function selectsCtfAndMegin(tc)
            ty = {'MEG','MEG MAG','MEG GRAD','MEG REF','EEG','EOG','Stim','meg grad'};
            tc.verifyEqual(rheome.sensors.ismeg(ty), logical([1 1 1 0 0 0 0 1]));
        end

        function scaleSensorsSelectsMegin(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            old = getenv('RHEOME_DATA');  setenv('RHEOME_DATA', d);  tc.addTeardown(@setenv, 'RHEOME_DATA', old);
            scaleSynthSubject(fullfile(d, 'synth_megin'));
            f = fullfile(d, 'synth_megin', 'study.mat');  s = load(f);
            ty = repmat({'MEG GRAD', 'MEG GRAD', 'MEG MAG'}, 1, 14);  ty = ty(1:40);  ty{40} = 'MEG REF';
            s.chan.Type = ty;  save(f, '-struct', 's');
            S = rheome.scale.sensors('synth_megin');
            tc.verifyEqual(S.iSel(:)', 1:39);
            tc.verifyEqual(S.chT(:)', ty(1:39));
            tc.verifyEqual(size(S.Res.ImagingKernel, 2), 39);
        end
    end
end
