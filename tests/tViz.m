classdef tViz < matlab.unittest.TestCase
% TVIZ  rheome.viz (the visualization package's builders) and the bundle's labellings and point sets: tiles by level as a
% labelling, critical points as one set a type at their frames, the marks sidecar, and rheome.view handing them on.
%
% Author: Diellor Basha, 2026
    properties
        S = struct('Vertices', [0 0 0; 1 0 0; 0 1 0; 0 0 1], 'Faces', [1 2 3; 1 3 4; 1 4 2; 2 4 3]);
    end
    methods (Test)
        function tilesCoverEveryVertexAtALevel(tc)
            [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', 0.07 * V, 'Faces', F);
            T = rheome.geom.tree(S, MaxDepth=3);
            L = rheome.viz.tiles(T, 2);
            tc.verifyEqual(numel(L.entries), 4);
            tc.verifyEqual(numel(L.codes), size(V, 1));
            tc.verifyTrue(all(L.codes >= 0));
            tc.verifyEqual(unique(double(L.codes))', 0:3);
            % nested: a level-3 tile lies in one level-2 tile
            L3 = rheome.viz.tiles(T, 3);
            for c = unique(L3.codes)', tc.verifyNumElements(unique(L.codes(L3.codes == c)), 1); end
        end
        function tilesOfHemisphereTrees(tc)
            [V, F] = rheome.geom.icosphere(2);  S = struct('Vertices', V, 'Faces', F);
            T = rheome.geom.tree(S, MaxDepth=1);  n = size(V, 1);
            L = rheome.viz.tiles({T, T}, 1, Vertices={1:n, n + (1:n)}, NV=2 * n);
            tc.verifyEqual(numel(L.entries), 4);
            tc.verifyTrue(all(L.codes >= 0));
            tc.verifyError(@() rheome.viz.tiles({T, T}, 1), 'rheome:viz:tiles');
        end
        function pointsOneSetAType(tc)
            c1 = struct('pos', [0.9 0 0; 0 0 1], 'type', {{'vortex'; 'source'}}, 'charge', [1; 1]);
            c2 = struct('pos', [0 1 0], 'type', {{'vortex'}}, 'charge', -1);
            P = rheome.viz.points(struct('cp', [c1 c2], 'frames', [3 7]), tc.S, Tokens="none");
            tc.verifyEqual({P.name}, {'critical: vortex', 'critical: source'});
            tc.verifyEqual(P(1).vertices', [2 3]);
            tc.verifyEqual(P(1).frames', [3 7]);
            tc.verifyEqual(P(1).signs', [1 -1]);
            tc.verifyEqual(P(2).frames, 3);
            tc.verifyEqual(P(1).color, '');   % no tokens: the viewer's default
            Q = rheome.viz.points(c1, tc.S);  % one frame: every frame
            tc.verifyEqual([Q.frames], [0 0]);
        end
        function pointsTakeTheTokensColours(tc)
            f = [tempname '.json'];  c = onCleanup(@() delete(f));
            fid = fopen(f, 'w');  fwrite(fid, '{"colour":{"semantic":{"flow-source":{"$value":"#1b9e77"},"flow-vortex":{"$value":"#c2378a"}}}}');  fclose(fid);
            cp = struct('pos', [1 0 0; 0 1 0], 'type', {{'sink'; 'vortex'}}, 'charge', [1; 1]);
            P = rheome.viz.points(cp, tc.S, Tokens=f);
            tc.verifyEqual({P.color}, {'#1b9e77', '#c2378a'});
        end
        function bundleWritesLabellingsAndPoints(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));
            L = struct('name', 'tiles', 'codes', [0 0 1 -1]', 'entries', struct('code', {0, 1}, 'name', {'a', 'b'}));
            P = struct('name', 'cp', 'vertices', [1 4]', 'frames', [0 2]', 'signs', [1 -1]', 'color', '#1b9e77');
            m = rheome.io.write.bundle(d, tc.S, rand(4, 2), 'Labellings', L, 'Points', P);
            j = jsondecode(fileread(m));
            tc.verifyEqual(j.labellings.name, 'tiles');
            tc.verifyEqual(rd(d, j.labellings.codes, 'int32')', int32([0 0 1 -1]));
            tc.verifyEqual(numel(j.labellings.entries), 2);
            tc.verifyEqual(rd(d, j.points.vertices, 'int32')', int32([0 3]));   % 0-based
            tc.verifyEqual(rd(d, j.points.frames, 'int32')', int32([-1 1]));    % 0: every frame -> -1
            tc.verifyEqual(rd(d, j.points.signs, 'int32')', int32([1 -1]));
            tc.verifyEqual(j.points.color, '#1b9e77');
        end
        function bundleRefusesBadSets(tc)
            d = tempname;  c = onCleanup(@() rmdir(d, 's'));  mkdir(d);
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(4, 2), 'Labellings', struct('name', 'x', 'codes', [0 1])), 'rheome:bundle:labellings');
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(4, 2), 'Points', struct('name', 'x', 'vertices', 5)), 'rheome:bundle:points');
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(4, 2), 'Points', struct('name', 'x', 'vertices', 1, 'frames', 3)), 'rheome:bundle:points');
            tc.verifyError(@() rheome.io.write.bundle(d, tc.S, rand(4, 2), 'Points', struct('vertices', 1)), 'rheome:bundle:points');
        end
        function viewHandsThemOn(tc)
            d = tempname;  old = getenv('RHEOME_VIEWER');  setenv('RHEOME_VIEWER', d);
            c = onCleanup(@() cleanup(d, old));
            L = struct('name', 'tiles', 'codes', [0 0 1 1]');
            out = rheome.view(tc.S, rand(4, 2), 'Viewer', 'none', 'Labellings', L, 'Labelling', "tiles", ...
                'Points', struct('name', 'cp', 'vertices', 2));
            j = jsondecode(fileread(out.manifest));
            tc.verifyEqual(j.settings.labelling, 'tiles');
            tc.verifyEqual(j.labellings.name, 'tiles');
            tc.verifyEqual(j.points.name, 'cp');
        end
        function marksSidecar(tc)
            d = tempname;  mkdir(d);  c = onCleanup(@() rmdir(d, 's'));
            m = struct('format', 'nxr.marks@1', 'surface', struct('name', 'pial', 'folded', true), ...
                'overlays', {{struct('overlay', 'arrows', 'style', 'ink', 'key', 'arrows: J')}});
            f = rheome.viz.marks(m, fullfile(d, 'F2b.png'));
            tc.verifyEqual(f, fullfile(d, 'F2b.marks.json'));
            j = jsondecode(fileread(f));
            tc.verifyEqual(j.overlays.key, 'arrows: J');
            tc.verifyError(@() rheome.viz.marks(struct('format', 'x'), fullfile(d, 'a.png')), 'rheome:viz:marks');
        end
    end
end

function x = rd(d, name, cls)
    fid = fopen(fullfile(d, name), 'r', 'ieee-le');  x = fread(fid, Inf, ['*' cls]);  fclose(fid);
end

function cleanup(d, old)
    setenv('RHEOME_VIEWER', old);
    if isfolder(d), rmdir(d, 's'); end
end
