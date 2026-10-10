classdef tViewBundle < matlab.unittest.TestCase
% TVIEWBUNDLE  rheome.io.write.bundle writes rheome.bundle@1 as the Editor reads it, and rheome.view needs no viewer:
% 'none' writes beneath RHEOME_VIEWER (one bundle at a time), 'show' draws with rheome.show.
%
% Author: Diellor Basha, 2026
    properties
        S = struct('Vertices', [0 0 0; 1 0 0; 0 1 0; 0 0 1], 'Faces', [1 2 3; 1 3 4; 1 4 2; 2 4 3]);
    end
    methods (Test)
        function writesTheFormat(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            I = rand(4, 3);  v = rand(12, 3);
            m = rheome.io.write.bundle(d, tc.S, struct('activity', I, 'velocity', v), 'Fs', 100);
            j = jsondecode(fileread(m));
            tc.verifyEqual(j.format, 'rheome.bundle@1');
            tc.verifyEqual([j.time.samples j.time.fs], [3 100]);
            tc.verifyEqual({j.fields.type}, {'real', 'vector'});
            tc.verifyEqual(reshape(rd(d, 'positions.bin', 'single'), 3, []), single(tc.S.Vertices'));
            tc.verifyEqual(reshape(rd(d, 'faces.bin', 'int32'), 3, []), int32(tc.S.Faces' - 1));
            tc.verifyEqual(reshape(rd(d, j.fields(2).values, 'single'), 12, 3), single(v));
        end
        function refusesABadResult(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));  mkdir(d);
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(5, 2)), 'rheome:bundle:rows');
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, struct('a', rand(4, 2), 'b', rand(4, 3))), 'rheome:bundle:time');
        end
        function noneWritesBeneathTheViewer(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            a = rheome.view(tc.S, rand(4, 2), 'Viewer', 'none');
            pause(0.002);
            b = rheome.view(tc.S, rand(12, 2), 'Viewer', 'none');
            tc.verifyTrue(startsWith(b.manifest, fullfile(d, 'bundles')));
            tc.verifyTrue(isfile(b.manifest));
            tc.verifyFalse(isfile(a.manifest));   % one bundle at a time
            tc.verifyEqual(b.viewer, 'none');
        end
        function showDrawsWithoutTheViewer(tc)
            out = rheome.view(tc.S, struct('speed', rand(4, 2)), 'Viewer', 'show');
            c = onCleanup(@() close(out.figure));
            tc.verifyEqual(out.viewer, 'show');
            tc.verifyClass(out.figure, 'matlab.ui.Figure');
            tc.verifyEmpty(out.manifest);
        end
    end
end

function x = rd(d, f, cls)
    fid = fopen(fullfile(d, f), 'r', 'ieee-le');  x = fread(fid, inf, ['*' cls]);  fclose(fid);
end

function cleanup(d, old)
    setenv('RHEOME_VIEWER', old);
    if isfolder(d), rmdir(d, 's'); end
end
