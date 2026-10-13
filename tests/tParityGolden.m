classdef tParityGolden < matlab.unittest.TestCase
% TPARITYGOLDEN  parityWriteGolden writes rheome.golden@1 as the viewer's parity test reads it: every array of
% golden.json in its file, float64 little-endian, row-major in the shape it names; I is the planted bump.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function writesTheFormat(tc)
            d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
            [V, F] = rheome.geom.icosphere(3);  V = 0.07 * V;     % a cortex-sized sphere, as the viewer's fixture
            stem = fullfile(d, 'fx');
            put([stem '.positions.f32'], V', 'single');
            put([stem '.faces.u32'], F' - 1, 'uint32');
            out = fullfile(d, 'golden');
            parityWriteGolden(stem, out);
            j = jsondecode(fileread(fullfile(out, 'golden.json')));
            tc.verifyEqual(j.format, 'rheome.golden@1');
            nV = size(V, 1);  T = j.T;
            want = struct('Lx', [nV 1], 'Mlumped', [nV 1], 'MgalRowSum', [nV 1], 'I', [T nV], ...
                'flow_alpha1', [T-1 3*nV], 'flow_alpha0_01', [T-1 3*nV], 'div', [T-1 nV], 'curl', [T-1 nV]);
            for k = fieldnames(want)'
                a = j.arrays.(k{1});
                tc.verifyEqual(a.shape(:)', want.(k{1}), k{1});
                tc.verifyEqual(dir(fullfile(out, a.file)).bytes, 8 * prod(want.(k{1})), k{1});
            end
            fid = fopen(fullfile(out, 'I.bin'), 'r', 'ieee-le');  I = reshape(fread(fid, inf, 'double'), nV, T);  fclose(fid);
            Vs = double(single(V));
            I1 = exp(-sum((Vs - j.bump.c0(:)').^2, 2) / (2 * j.bump.sigma^2));
            tc.verifyEqual(I(:, 1), I1, 'AbsTol', 1e-12);       % row-major: row t is frame t
            fid = fopen(fullfile(out, 'Mlumped.bin'), 'r', 'ieee-le');  M = fread(fid, inf, 'double');  fclose(fid);
            tc.verifyEqual(sum(M), 4*pi*0.07^2, 'RelTol', 0.01);   % the lumped mass is the area
        end
    end
end

function put(f, x, cls)
    fid = fopen(f, 'w', 'ieee-le');  fwrite(fid, x, cls);  fclose(fid);
end
