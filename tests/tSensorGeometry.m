classdef tSensorGeometry < matlab.unittest.TestCase
% Author: Diellor Basha, 2026

    methods (Test)

        function linearProbeHasTheRequestedPitchAndSpan(tc)
            arr = rheome.sensors.linear('NumSites', 32, 'Pitch', 100e-6);
            tc.verifyEqual(arr.nCh, 32);
            tc.verifyEqual(arr.Dim, 1);
            tc.verifyEqual(arr.Pitch, 100e-6, 'RelTol', 1e-12);
            % 32 sites at 100 um span 31 gaps, not 32.
            tc.verifyEqual(arr.Aperture, 31*100e-6, 'RelTol', 1e-12);
        end

        function utahApertureIsTheDiagonalNotTheSide(tc)
            % Aperture is the MAXIMUM PAIRWISE distance, so a square grid's aperture is
            % its diagonal. Getting this wrong understates the window ceiling by 41%.
            arr = rheome.sensors.grid('utah');
            tc.verifyEqual(arr.nCh, 100);
            tc.verifyEqual(arr.Pitch, 400e-6, 'RelTol', 1e-12);
            tc.verifyEqual(arr.Aperture, 9*400e-6*sqrt(2), 'RelTol', 1e-12);
        end

        function ecogPresetIsEightByEightAtOneCentimetre(tc)
            arr = rheome.sensors.grid('ecog');
            tc.verifyEqual(arr.nCh, 64);
            tc.verifyEqual(arr.Pitch, 10e-3, 'RelTol', 1e-12);
            tc.verifyEqual(arr.Aperture, 7*10e-3*sqrt(2), 'RelTol', 1e-12);
        end

        function pitchIsMedianNearestNeighbourNotMeanSpacing(tc)
            % Every interior nearest-neighbour distance in a lattice is exactly the pitch,
            % so the median must be the pitch even though the MEAN pairwise distance is far
            % larger. This pins the definition the window depends on.
            arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 2e-3);
            tc.verifyEqual(arr.Pitch, 2e-3, 'RelTol', 1e-12);
        end

        function gridLabelsAreUniqueAndCountMatches(tc)
            arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 1e-3);
            tc.verifyEqual(arr.nCh, 12);
            tc.verifyNumElements(arr.Labels, 12);
            tc.verifyNumElements(unique(arr.Labels), 12);
        end

        function gridIsPlanarAndOrderedByNdgrid(tc)
            arr = rheome.sensors.grid('Size', [3 4], 'Pitch', 1e-3);
            tc.verifyEqual(arr.Pos(:,3), zeros(12,1), 'AbsTol', 1e-15);
            tc.verifyEqual(arr.Pos(1,:), [0 0 0], 'AbsTol', 1e-15);
            tc.verifyEqual(arr.Pos(2,:), [1e-3 0 0], 'AbsTol', 1e-15);   % rows vary fastest
            tc.verifyEqual(arr.Pos(4,:), [0 1e-3 0], 'AbsTol', 1e-15);
        end

        function anUnknownGridPresetErrors(tc)
            tc.verifyError(@() rheome.sensors.grid('nosucharray'), 'sensors:grid:preset');
        end

        function eegFromAChannelStructRejectsPrefixMatchedTypes(tc)
            % ⚠ THE SAME BUG THAT LET 26 REFERENCE GRADIOMETERS INTO THE MEG ARRAY.
            % strncmpi(Type,'EEG',3) matches 'EEG REF' too, so a prefix match would keep 3
            % channels here and put a reference electrode into the array geometry. Under the
            % old implementation nCh would read 3 and this test fails on the count alone.
            C = struct();
            C.Channel = struct('Name', {'Fz','Cz','REF1'}, ...
                               'Type', {'EEG','EEG','EEG REF'}, ...
                               'Loc',  {[0;0;0.09], [0.01;0;0.089], [0.05;0;0.07]});
            C.Type = {C.Channel.Type};
            C.Name = {C.Channel.Name};
            C.nCh  = numel(C.Channel);

            arr = rheome.sensors.eeg('Channel', C);
            tc.verifyEqual(arr.nCh, 2);
            tc.verifyEqual(arr.Kind, 'eeg');
            tc.verifyEqual(arr.Dim, 2);
            tc.verifyEqual(arr.Labels, {'Fz','Cz'});
            tc.verifyFalse(any(strcmp(arr.Labels, 'REF1')));
        end

        function eegWithNoPlainEegChannelsErrors(tc)
            % A struct of ONLY 'EEG REF' must error, not silently build a reference array.
            C = struct('Channel', struct('Name', 'REF1', 'Type', 'EEG REF', ...
                                         'Loc', [0.05;0;0.07]), ...
                       'Type', {{'EEG REF'}}, 'Name', {{'REF1'}}, 'nCh', 1);
            tc.verifyError(@() rheome.sensors.eeg('Channel', C), 'sensors:eeg:noChannels');
        end

        function eegCapHasTheRequestedCountAndSitsOnItsSphere(tc)
            arr = rheome.sensors.eeg('NumChannels', 64, 'Radius', 0.09, 'CapAngle', 120);
            tc.verifyEqual(arr.nCh, 64);
            tc.verifyEqual(arr.Dim, 2);
            tc.verifyEqual(vecnorm(arr.Pos, 2, 2), 0.09*ones(64,1), 'RelTol', 1e-12);
        end

        function eegCapApertureIsTheChordAcrossTheCap(tc)
            % A 120-degree cap extends past the equator (cap > 90 degrees), so it contains
            % points whose ANTIPODES also lie in the cap: a point at colatitude theta and one
            % at colatitude 180-theta are both within the cap whenever 180-cap <= theta <= cap,
            % a non-empty range exactly when cap >= 90 degrees. Such a pair is exactly
            % antipodal, so the widest baseline the cap can contain is the full sphere
            % DIAMETER 2*R -- not the rim-to-rim chord 2*R*sin(cap), which only bounds the
            % aperture when cap <= 90 degrees. Quasi-uniform sampling approaches the diameter
            % from below as N grows.
            R = 0.09;
            arr = rheome.sensors.eeg('NumChannels', 64, 'Radius', R, 'CapAngle', 120);
            tc.verifyLessThanOrEqual(arr.Aperture, 2*R + 1e-12);
            tc.verifyGreaterThan(arr.Aperture, 0.9 * 2*R);
        end

        function eegCapIsQuasiUniformSoPitchIsTight(tc)
            % Spherical-Fibonacci sampling: no sensor should have a nearest neighbour more
            % than 1.6x the median. A clustered layout would fail this and would make the
            % single 'pitch' number meaningless.
            arr = rheome.sensors.eeg('NumChannels', 64);
            D = tSensorGeometry.distanceMatrix(arr.Pos);
            D(1:size(D,1)+1:end) = Inf;
            nn = min(D, [], 2);
            tc.verifyLessThan(max(nn)/median(nn), 1.6);
        end

        function megFromAChannelStructKeepsOnlyMegChannels(tc)
            % Built by hand so the test needs no cached dataset: 4 MEG on a shell plus
            % 2 EEG, 1 STIM and 1 MEG REF that must all be dropped.
            %
            % ⚠ 'MEG REF' IS THE ONE THAT MATTERS. strncmpi(Type, 'MEG', 3) matches it too --
            % a CTF reference gradiometer sits in a rack behind the head, not over the
            % helmet, and including it corrupts every geometric quantity the array reports.
            % A fixture without this type cannot catch a three-character prefix match; it
            % must be here, not merely 'EEG'/'STIM', for this test to guard what it claims.
            ang = [0 0.4 0.8 1.2];
            C = struct();
            C.Channel = struct('Name', {}, 'Type', {}, 'Loc', {});
            for i = 1:4
                C.Channel(end+1) = struct('Name', sprintf('M%d', i), 'Type', 'MEG', ...
                    'Loc', [0.1*cos(ang(i)); 0.1*sin(ang(i)); 0.05] * [1 1]);
            end
            C.Channel(end+1) = struct('Name', 'E1',  'Type', 'EEG',  'Loc', [0;0;0.09]*[1 1]);
            C.Channel(end+1) = struct('Name', 'E2',  'Type', 'EEG',  'Loc', [0;0.01;0.09]*[1 1]);
            C.Channel(end+1) = struct('Name', 'TRG', 'Type', 'STIM', 'Loc', []);
            C.Channel(end+1) = struct('Name', 'BG1', 'Type', 'MEG REF', ...
                'Loc', [0.3; 0; 0.05] * [1 1]);
            C.Type = {C.Channel.Type};
            C.Name = {C.Channel.Name};
            C.nCh  = numel(C.Channel);

            arr = rheome.sensors.meg('Channel', C);
            tc.verifyEqual(arr.nCh, 4);
            tc.verifyEqual(arr.Kind, 'meg');
            tc.verifyEqual(arr.Dim, 2);
            tc.verifyEqual(arr.Labels, {'M1','M2','M3','M4'});
            tc.verifyFalse(any(strcmp(arr.Labels, 'BG1')));
        end

        function megUsesThePrimaryCoilNotTheCoilMean(tc)
            % A CTF gradiometer's Loc has TWO columns: the pickup coil and the reference
            % coil ~5 cm above it. Averaging them puts the sensor in mid-air, halfway to
            % the reference -- so the primary coil, column 1, is the position.
            C = struct();
            C.Channel = struct('Name', {'M1','M2'}, 'Type', {'MEG','MEG'}, ...
                'Loc', {[0 0; 0 0; 0.10 0.15], [0.01 0.01; 0 0; 0.10 0.15]});
            C.Type = {C.Channel.Type};  C.Name = {C.Channel.Name};  C.nCh = 2;
            arr = rheome.sensors.meg('Channel', C);
            tc.verifyEqual(arr.Pos(1,:), [0 0 0.10], 'AbsTol', 1e-15);
            tc.verifyEqual(arr.Pos(2,:), [0.01 0 0.10], 'AbsTol', 1e-15);
        end

        function megWithNoMegChannelsErrors(tc)
            C = struct('Channel', struct('Name', 'E1', 'Type', 'EEG', 'Loc', [0;0;0.09]), ...
                       'Type', {{'EEG'}}, 'Name', {{'E1'}}, 'nCh', 1);
            tc.verifyError(@() rheome.sensors.meg('Channel', C), 'sensors:meg:noChannels');
        end

        function aProbeHasNoFacesAndThatIsTheResultNotAGap(tc)
            % A winding number is summed around a face (a 2-cell), and a chain admits no
            % faces, so a winding number is UNDEFINED on it. Empty faces is how the machinery
            % states that, rather than special-casing it later.
            G = rheome.sensors.graph(rheome.sensors.linear('NumSites', 16));
            tc.verifyEmpty(G.Faces);
        end

        function aFlatGridTriangulatesToTwoFacesPerCell(tc)
            % An m x n lattice has (m-1)*(n-1) cells and a triangulation splits each in two.
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [5 6], 'Pitch', 1e-3));
            tc.verifyEqual(size(G.Faces, 2), 3);
            tc.verifyEqual(size(G.Faces, 1), 2*4*5);
        end

        function everyFaceIndexIsAValidSensor(tc)
            G = rheome.sensors.graph(rheome.sensors.grid('Size', [5 6], 'Pitch', 1e-3));
            tc.verifyGreaterThanOrEqual(min(G.Faces(:)), 1);
            tc.verifyLessThanOrEqual(max(G.Faces(:)), G.nV);
            tc.verifyEqual(numel(unique(G.Faces(:))), G.nV);   % no orphan sensors
        end

        function aCurvedCapTriangulatesWithoutClosingTheBottom(tc)
            % convhull of a CAP returns a closed hull, and the faces that close it span the
            % opening -- huge, spurious triangles. They must be dropped, or the readout
            % reports loops that do not exist on the array.
            arr = rheome.sensors.eeg('NumChannels', 64);
            G   = rheome.sensors.graph(arr);
            tc.verifyGreaterThan(size(G.Faces, 1), 90);
            e = tSensorGeometry.faceEdgeLengths(G.Vertices, G.Faces);
            tc.verifyLessThan(max(e)/median(e), 3.0 + 1e-9);
        end

        function faceCountObeysEulerForADisc(tc)
            % A triangulated disc: V - E + F = 1. A hull that still had its closing faces
            % would be a sphere (chi = 2) and fail this.
            G  = rheome.sensors.graph(rheome.sensors.eeg('NumChannels', 64));
            F  = G.Faces;
            ed = sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2);
            nE = size(unique(ed, 'rows'), 1);
            tc.verifyEqual(G.nV - nE + size(F,1), 1);
        end

        function aSensorAtTheSphereCentreErrorsInsteadOfGarbageFaces(tc)
            % A near-zero-norm radial vector normalises to an arbitrary-but-finite
            % direction that convhull would silently accept. Six sensors on an octahedron
            % (radius R, centred at the origin) pin the sphere fit to near-exact c = 0,
            % R = 0.09; a 7th sensor placed 1e-12 m from the origin then lands within
            % 1e-9*R of the fitted centre and must be refused, not triangulated.
            R = 0.09;
            oct = R * [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
            P = [oct; 1e-12 0 0];
            C = struct();
            C.Channel = struct('Name', {}, 'Type', {}, 'Loc', {});
            for i = 1:size(P, 1)
                C.Channel(end+1) = struct('Name', sprintf('M%d', i), 'Type', 'MEG', ...
                    'Loc', P(i,:).' * [1 1]);
            end
            C.Type = {C.Channel.Type};
            C.Name = {C.Channel.Name};
            C.nCh  = numel(C.Channel);
            arr = rheome.sensors.meg('Channel', C);

            tc.verifyError(@() rheome.sensors.graph(arr), 'sensors:faces:sensorAtSphereCentre');
        end

    end

    methods (Static)
        function D = distanceMatrix(P)
            G = P*P';
            D = sqrt(max(diag(G) + diag(G).' - 2*G, 0));
        end

        function e = faceEdgeLengths(V, F)
            e = [vecnorm(V(F(:,1),:) - V(F(:,2),:), 2, 2); ...
                 vecnorm(V(F(:,2),:) - V(F(:,3),:), 2, 2); ...
                 vecnorm(V(F(:,3),:) - V(F(:,1),:), 2, 2)];
        end
    end
end

% Author: Diellor Basha, 2026
