classdef tSurfaceLadder < matlab.unittest.TestCase
% The atlas's tiles in MATLAB: rheome.geom.surfaceladder is bit-exact to nsp's atlas ladder
% (nxr_convert/atlas/ladder.py), and rheome.geom.tree's Ladder / Codes routes carry its codes.
%
% Golden: ladderGolden.mat (ladderWriteMeshes.m -> ladderCaptureGolden.py, ladder.py under CPython
% 3.12) and pysortGolden.mat (pysortCaptureGolden.py -> pysortWriteGolden.m). On the ten VR1
% PREVENT-AD hemispheres the port is identical to ladder.py at all 12 levels (request d03bd647).
%
% Author: Diellor Basha, 2026

    properties
        G
    end

    methods (TestClassSetup)
        function golden(tc)
            tc.G = load(fullfile(fileparts(mfilename('fullpath')), 'ladderGolden.mat'));
        end
    end

    methods (Test)
        function pysortIsPythonsSort(tc)
            % Python's sorted() with an inconsistent (NaN) comparator: the steps decide the result
            g = load(fullfile(fileparts(mfilename('fullpath')), 'pysortGolden.mat'));
            for i = 1:numel(g.keys)
                tc.verifyEqual(rheome.utils.pysort(g.keys{i}), g.perm{i}, sprintf('case %d', i));
            end
        end

        function ladderIsBitExactToNsp(tc)
            for i = 1:numel(tc.G.meshes)
                m = tc.G.meshes{i};
                lab = rheome.geom.surfaceladder(m.V, m.F);
                tc.verifyEqual(double(lab), tc.G.labels{i}, sprintf('mesh %d', i));
            end
        end

        function ladderLevelsAreArithmetic(tc)
            m = tc.G.meshes{2};  [lab, mass] = rheome.geom.surfaceladder(m.V, m.F);
            L = size(lab, 2) - 1;
            for l = 0:L
                tc.verifyEqual(lab(:, l + 1), int32(floor(double(lab(:, end)) / 2^(L - l))));
            end
            ta = accumarray(double(lab(:, 5)) + 1, mass);          % equal mass at level 4 (160 vertices a tile)
            tc.verifyLessThan(std(ta) / mean(ta), 0.05);
        end

        function treeReadsTheLadder(tc)
            m = tc.G.meshes{2};  S = struct('Vertices', m.V, 'Faces', m.F);  D = 6;
            T = rheome.geom.tree(S, Method="Ladder", MaxDepth=D);
            tc.verifyEqual(height(T), 2^(D + 1) - 1);
            tc.verifyEqual(T.node_id', 1:2^(D + 1) - 1);                       % heap: children 2k, 2k+1
            tc.verifyEqual(T.parent_id, floor(T.node_id / 2));
            ref = tc.G.labels{2}(:, D + 1);
            leaves = T(T.is_leaf, :);
            for k = 1:height(leaves)                                           % code = node_id - 2^depth
                tc.verifyEqual(leaves.members{k}, find(ref == leaves.node_id(k) - 2^D)');
            end
            T2 = rheome.geom.tree(S, Codes=ref, MaxDepth=D);                   % read, not re-built
            tc.verifyEqual(T2, T);
            G = rheome.geom.tiles(T, S, D);
            tc.verifyTrue(G.heap);
        end

        function treeRejectsWhatTheLadderCannotDo(tc)
            m = tc.G.meshes{1};  S = struct('Vertices', m.V, 'Faces', m.F);
            tc.verifyError(@() rheome.geom.tree(S, Method="Ladder", MinArea=1e-3), 'rheome:geom:tree');
            tc.verifyError(@() rheome.geom.tree(S, Codes=zeros(size(m.V, 1), 1)), 'rheome:geom:tree');      % no MaxDepth
            tc.verifyError(@() rheome.geom.tree(S, Codes=zeros(size(m.V, 1), 1), MaxDepth=1), 'rheome:geom:tree');  % empty tile
        end
    end
end
