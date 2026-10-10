classdef tFsSurf < matlab.unittest.TestCase
% TFSSURF  rheome.io.read.fssurf reads the FreeSurfer triangle format: a written mesh round-trips, a bad magic is refused.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function roundTrip(tc)
            V = [0 0 0; 10 0 0; 0 10 0; 0 0 10];  F = [1 2 3; 1 3 4; 1 4 2; 2 4 3];
            f = [tempname '.surf'];  c = onCleanup(@() delete(f));
            fid = fopen(f, 'w', 'ieee-be');
            fwrite(fid, [255 255 254], 'uchar');  fprintf(fid, 'created by test\n\n');
            fwrite(fid, [4 4], 'int32');  fwrite(fid, V.', 'float32');  fwrite(fid, (F-1).', 'int32');  fclose(fid);
            S = rheome.io.read.fssurf(f);
            tc.verifyEqual(S.Vertices, 1e-3*V, 'AbsTol', 1e-12);
            tc.verifyEqual(S.Faces, F);
            tc.verifyEqual([S.nV S.nF], [4 4]);
            tc.verifyEqual(S.Comment, 'created by test');
        end
        function refusesOtherFormats(tc)
            f = [tempname '.surf'];  c = onCleanup(@() delete(f));
            fid = fopen(f, 'w');  fwrite(fid, 'not a surface');  fclose(fid);
            tc.verifyError(@() rheome.io.read.fssurf(f), 'io:read:fssurf:format');
        end
    end
end

% Author: Diellor Basha, 2026
