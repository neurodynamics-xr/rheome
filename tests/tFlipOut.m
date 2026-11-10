classdef tFlipOut < matlab.unittest.TestCase
% rheome.geom.draw.flipOut against the viewer's FlipOut (nxr-compute WASM, geometry-central): every pair of
% tests/goldens/flipout has the same points (1e-9 m), the same vertex/edge/face placement and the same length
% (1e-9 relative). fsaverage5 lh pial runs when NXR_FSAVERAGE5 is set.
%
% Author: Diellor Basha, 2026
    properties (TestParameter)
        mesh = {'bumpy642', 'bumpy2562', 'fsaverage5_lh_pial'};
    end
    methods (Test)
        function sameAsTheViewer(t, mesh)
            R = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'goldens', 'flipout', 'flipout-reference.json')));
            m = R.(mesh);
            if isfield(m, 'vertices')
                S = struct('Vertices', m.vertices, 'Faces', m.faces + 1);
            else
                root = getenv('NXR_FSAVERAGE5');
                t.assumeFalse(isempty(root), 'NXR_FSAVERAGE5 is not set: fsaverage5 skipped.');
                S = readSurf(fullfile(root, 'surf', 'lh.pial'));
            end
            cache = [];
            for j = 1:numel(m.paths)
                p = m.paths(j);  at = sprintf('%s %d -> %d', mesh, p.from, p.to);
                [P, info] = rheome.geom.draw.flipOut(S, p.from + 1, p.to + 1, cache);  cache = info.cache;
                t.assertEqual(size(P, 1), p.nPoints, [at ': number of points']);
                t.verifyEqual(P, reshape(p.positions, 3, []).', 'AbsTol', 1e-9, [at ': positions']);
                t.verifyEqual(info.vertices, reshape(p.vertices, 3, []).' + 1, [at ': vertices']);
                t.verifyEqual(info.weights, reshape(p.weights, 3, []).', 'AbsTol', 1e-9, [at ': weights']);
                t.verifyEqual(info.length, p.length, 'RelTol', 1e-9, [at ': length']);
            end
        end
        function theSameVertexAndBadInput(t)
            [V, F] = rheome.geom.icosphere(2);  S = struct('Vertices', V, 'Faces', F);
            [P, info] = rheome.geom.draw.flipOut(S, 5, 5);
            t.verifyEqual(P, V(5,:));  t.verifyEqual(info.length, 0);
            t.verifyError(@() rheome.geom.draw.flipOut(S, 0, 5), 'rheome:flipOut:vertex');
        end
    end
end

function S = readSurf(f)
% the FreeSurfer surface as the reference read it: float32 mm / 1000 (rheome.io.read.fssurf scales by 1e-3, a
% different last bit)
    fid = fopen(f, 'r', 'ieee-be');  c = onCleanup(@() fclose(fid));
    fread(fid, 3, 'uchar');  fgetl(fid);  fgetl(fid);
    n = fread(fid, 2, 'int32');
    V = reshape(double(fread(fid, n(1)*3, 'float32')), 3, []).' / 1000;
    F = reshape(fread(fid, n(2)*3, 'int32'), 3, []).' + 1;
    S = struct('Vertices', V, 'Faces', F);
end

% Author: Diellor Basha, 2026
