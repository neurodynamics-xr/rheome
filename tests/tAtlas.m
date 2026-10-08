classdef tAtlas < matlab.unittest.TestCase
% The ROI axis: rheome.io.read.atlas / rheome.import.atlas / rheome.load.atlas.
%
% Three traps are pinned here because each one was found in the real Brainstorm files
% and each fails SILENTLY rather than loudly:
%   * scout Vertices are stored as ROW vectors -- vertcat over scouts errors outright
%   * scouts carry Label, NOT Name (unlike almost every other Brainstorm struct)
%   * a parcellation does NOT cover the surface (Desikan-Killiany: 18743 of 20484).
%     Treating it as a partition quietly drops the medial wall from every total.
%
% Author: Diellor Basha, 2026

    properties
        SurfFile        % synthetic Brainstorm surface
        DsName = '_tmp_atlas_test'
        DsDir
        nV = 10
    end

    methods (TestClassSetup)
        function writeSurface(tc)
            d = tempname;  mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            tc.SurfFile = fullfile(d, 'tess_cortex_test.mat');

            Vertices = rand(tc.nV, 3);                     %#ok<NASGU>
            Faces    = [1 2 3; 4 5 6; 7 8 9];              %#ok<NASGU>

            empt = struct('Vertices', {}, 'Seed', {}, 'Color', {}, 'Label', {}, ...
                          'Function', {}, 'Region', {}, 'Handles', {});
            Atlas = struct('Name', {'User scouts'}, 'Scouts', {empt});

            % Vertices deliberately stored as ROW vectors, as Brainstorm does
            Atlas(2).Name   = 'Test-Parcellation';
            Atlas(2).Scouts = [ ...
                tc.scout('alpha L', 'LT', 1:3, 2), ...
                tc.scout('beta R',  'RF', 4:6, 5), ...
                tc.scout('gamma L', 'LO', 7:8, 7)];

            Atlas(3).Name   = 'Overlapping';
            Atlas(3).Scouts = [tc.scout('one L','LT',1:5,1), tc.scout('two L','LT',5:9,6)];

            Atlas(4).Name   = 'Structures';
            Atlas(4).Scouts = [tc.scout('Cortex L','LT',1:5,1), tc.scout('Cortex R','RT',6:10,6)];

            Comment = 'synthetic cortex';                  %#ok<NASGU>
            save(tc.SurfFile, 'Vertices', 'Faces', 'Atlas', 'Comment');

            % a cache entry pointing at it, for the import/load round trip
            tc.DsDir = fullfile(rheome.load.root(), tc.DsName);
            tc.assumeFalse(exist(tc.DsDir,'dir') == 7, 'scratch dataset already exists');
            mkdir(tc.DsDir);
            tc.addTeardown(@() rmdir(tc.DsDir, 's'));
            S = rheome.io.read.surface(tc.SurfFile);              %#ok<NASGU>
            cortexFile = tc.SurfFile;                      %#ok<NASGU>
            L = []; M = [];                                %#ok<NASGU>
            save(fullfile(tc.DsDir, 'surface.mat'), 'S', 'L', 'M', 'cortexFile', '-v7.3');
        end
    end

    methods
        function s = scout(~, label, region, verts, seed)
            s = struct('Vertices', verts(:).', 'Seed', seed, 'Color', [1 0 0], ...
                       'Label', label, 'Function', 'Mean', 'Region', region, 'Handles', []);
        end
    end

    methods (Test)

        function readsEveryAtlasIncludingEmptyOnes(tc)
            % An empty atlas must survive: a caller selecting by name should find it
            % rather than get a confusing "no such atlas".
            A = rheome.io.read.atlas(tc.SurfFile);
            tc.verifyEqual(numel(A), 4);
            tc.verifyEqual({A.Name}, {'User scouts','Test-Parcellation','Overlapping','Structures'});
            tc.verifyEqual(A(1).nScout, 0);
            tc.verifyEqual(A(2).nScout, 3);
        end

        function verticesAreNormalisedToColumns(tc)
            % Stored as rows; vertcat over scouts errors unless normalised.
            A = rheome.io.read.atlas(tc.SurfFile);
            v = A(2).Vertices;
            tc.verifySize(v{1}, [3 1]);
            tc.verifyEqual(v{1}, (1:3)');
            tc.verifyWarningFree(@() vertcat(v{:}));
        end

        function scoutLabelsComeFromLabelNotName(tc)
            A = rheome.io.read.atlas(tc.SurfFile);
            tc.verifyEqual(A(2).Label, {'alpha L','beta R','gamma L'});
        end

        function hemisphereIsDerivedFromTheRegionCode(tc)
            A = rheome.io.read.atlas(tc.SurfFile);
            tc.verifyEqual(A(2).Region, {'LT','RF','LO'});
            tc.verifyEqual(A(2).Hemi,   {'L','R','L'});
        end

        function membershipMatrixSelectsEachScoutsVertices(tc)
            % This is the ROI axis in operator form: [nScout x nV], so a reduction over
            % ROIs is one sparse multiply rather than a loop.
            A = rheome.io.read.atlas(tc.SurfFile);
            Mm = A(2).Membership;
            tc.verifySize(Mm, [3 tc.nV]);
            tc.verifyEqual(find(Mm(1,:)), 1:3);
            tc.verifyEqual(find(Mm(3,:)), 7:8);
            x = (1:tc.nV)';
            tc.verifyEqual(full(double(Mm) * x), [6; 15; 15], 'AbsTol', 1e-12);
        end

        function partialCoverageIsReportedNotHidden(tc)
            % 8 of 10 vertices -- exactly the medial-wall situation in the real files.
            A = rheome.io.read.atlas(tc.SurfFile);
            tc.verifyEqual(A(2).nCovered, 8);
            tc.verifyEqual(A(2).Coverage, 0.8, 'AbsTol', 1e-12);
            tc.verifyTrue(A(2).Disjoint);
        end

        function overlappingScoutsAreFlagged(tc)
            % Not an error -- but summing over ROIs double-counts, so the caller must know.
            A = rheome.io.read.atlas(tc.SurfFile);
            iOv = strcmp({A.Name}, 'Overlapping');
            tc.verifyFalse(A(iOv).Disjoint);
            tc.verifyEqual(A(iOv).nCovered, 9);
        end

        function errorsOnASurfaceWithoutAnAtlas(tc)
            f = fullfile(tempname); mkdir(f);
            tc.addTeardown(@() rmdir(f, 's'));
            g = fullfile(f, 'bare.mat');
            Vertices = rand(4,3); Faces = [1 2 3];   %#ok<NASGU>
            save(g, 'Vertices', 'Faces');
            tc.verifyError(@() rheome.io.read.atlas(g), 'io:read:atlas:noAtlas');
        end

        function importCachesAndLoadReturnsTheSameAtlases(tc)
            f = rheome.import.atlas(tc.DsName);
            tc.verifyEqual(exist(f, 'file'), 2);
            A = rheome.load.atlas(tc.DsName);
            tc.verifyEqual({A.Name}, {'User scouts','Test-Parcellation','Overlapping','Structures'});
            Aref = rheome.io.read.atlas(tc.SurfFile);
            tc.verifyEqual(A(2).Membership, Aref(2).Membership);
            tc.verifyEqual(A(2).Vertices,   Aref(2).Vertices);
        end

        function loadSelectsASingleAtlasByName(tc)
            rheome.import.atlas(tc.DsName);
            A = rheome.load.atlas(tc.DsName, 'Test-Parcellation');
            tc.verifyEqual(A.Name, 'Test-Parcellation');
            tc.verifyEqual(A.nScout, 3);
        end

        function loadErrorsClearlyOnAnUnknownAtlasName(tc)
            rheome.import.atlas(tc.DsName);
            tc.verifyError(@() rheome.load.atlas(tc.DsName, 'Nonexistent'), 'load:atlas:noSuchAtlas');
        end

    end
end

% Author: Diellor Basha, 2026
