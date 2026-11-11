classdef tPeakTiles < matlab.unittest.TestCase
% Peaks on a mesh and the tiles that address them: rheome.detect.peaks and rheome.geom.tiles.
%
% A peak is a vertex that beats its one-ring; a tile is a node of rheome.geom.tree turned into a graph.
% What is asserted is that planted peaks come back at their vertex with sub-edge refinement, that
% a plateau yields one peak, that a tile's boundary maximum is not reported as a peak, and that the
% tile adjacency is exactly the cut-edge structure of the partition.
%
% Author: Diellor Basha, 2026

    properties
        S; L; M; basis; R = 0.07            % head-sized sphere, lengths in metres
    end

    methods (TestClassSetup)
        function mesh(tc)
            [V, F] = rheome.geom.icosphere(3);                         % 642 vertices
            tc.S = struct('Vertices', V * tc.R, 'Faces', F, 'nV', size(V,1), 'nF', size(F,1));
            [tc.L, tc.M] = rheome.operators.laplace_beltrami(tc.S.Vertices, tc.S.Faces);
            tc.basis = rheome.eigen.modes(tc.L, tc.M, 120);
        end
    end

    methods (Test)

        % ---------------- rheome.detect.peaks ----------------

        function twoPlantedBlobsComeBackAtTheirVertices(tc)
            x = tc.blob(10, 0.012) + 0.7 * tc.blob(400, 0.012);
            P = rheome.detect.peaks(x, tc.S, MinHeight=0.3);
            tc.verifyEqual(sort(P.table.vertex), sort([10; 400]));
            tc.verifyEqual(P.table.vertex(1), 10);               % highest first
            tc.verifyEqual(P.count, 2);
        end

        function aPlateauYieldsExactlyOnePeak(tc)
            x = zeros(tc.S.nV, 1);  x(tc.ring(5)) = 1;  x(5) = 1;
            P = rheome.detect.peaks(x, tc.S, MinHeight=0.5);
            tc.verifyEqual(P.count, 1);
        end

        function relativeThresholdIsPerFrame(tc)
            b = tc.blob(10, 0.012);
            P = rheome.detect.peaks([b 1e-3*b], tc.S, MinHeight=0.5, Relative=true);
            tc.verifyEqual(P.count, [1 1]);
            P = rheome.detect.peaks([b 1e-3*b], tc.S, MinHeight=0.5);
            tc.verifyEqual(P.count, [1 0]);
        end

        function refinedPositionFollowsABlendBetweenVertices(tc)
            % ⭐ A blend of two neighbouring atoms moves the refined position monotonically between
            % them while the vertex answer stays put: the sub-edge refinement is what makes a slow
            % peak's speed readable at all.
            v1 = 10;  nb = tc.ring(v1);  v2 = nb(1);
            a1 = tc.blob(v1, 0.012);  a2 = tc.blob(v2, 0.012);
            w = linspace(0, 0.45, 6);
            X = a1 .* (1 - w) + a2 .* w;
            P = rheome.detect.peaks(X, tc.S, MinHeight=0.5, Relative=true);
            d = vecnorm(P.table.pos - tc.S.Vertices(v1,:), 2, 2);
            tc.verifyTrue(all(P.table.vertex == v1));
            tc.verifyTrue(all(diff(d) > 0), 'refined position must move toward the second vertex');
        end

        function widthGrowsWithTheBlob(tc)
            P1 = rheome.detect.peaks(tc.blob(10, 0.008), tc.S, MinHeight=0.5);
            P2 = rheome.detect.peaks(tc.blob(10, 0.016), tc.S, MinHeight=0.5);
            tc.verifyGreaterThan(P2.table.widthMM, 1.5 * P1.table.widthMM);
        end

        function aTileBoundaryMaximumIsNotAPeak(tc)
            % ⚠ The trap the design avoids: one blob, many tiles. Per-tile maxima would report a
            % "peak" in every tile; true peaks number one.
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=3);  G = rheome.geom.tiles(T, tc.S, 3);
            x = tc.blob(10, 0.02);
            P = rheome.detect.peaks(x, tc.S, Tiles=G, MinHeight=0.5);
            tc.verifyEqual(P.count, 1);
            tc.verifyEqual(P.table.tile, G.tileOf(10));
            tileMax = arrayfun(@(k) max(x(G.tileOf == k)), 1:numel(G.node_id));
            tc.verifyGreaterThan(sum(tileMax > 0), 1);          % the naive rule would report more
        end

        function featuresFeedTheTracker(tc)
            X = [tc.blob(10, 0.012) tc.blob(10, 0.012) tc.blob(10, 0.012) tc.blob(10, 0.012)];
            P = rheome.detect.peaks(X, tc.S, MinHeight=0.5);
            out = rheome.detect.track([], tc.S, struct('frameFeatures', {P.features}, 'types', {{'peak'}}, ...
                'minStrength', 0, 'highStrength', 0));
            tc.verifyEqual(out.summary.nTracks, 1);
            tc.verifyEqual(numel(out.tracks(1).frames), 4);
        end

        % ---------------- rheome.geom.edgegraph ----------------

        function edgegraphReadsTheGreatCircleWithinThreePercent(tc)
            % ⭐ The ruler a speed is read with. Edges alone zig-zag (+7% median); the intrinsic
            % diagonals bring it to ~+1.4%, and it never saturates the way the heat method does.
            V = tc.S.Vertices / tc.R;  rng(1);
            a = randi(tc.S.nV, 200, 1);  b = randi(tc.S.nV, 200, 1);
            t = tc.R * acos(max(-1, min(1, sum(V(a,:) .* V(b,:), 2))));
            g0 = rheome.geom.edgegraph(tc.S, Diagonals=false);  g1 = rheome.geom.edgegraph(tc.S);
            d0 = arrayfun(@(k) distances(g0, a(k), b(k)), 1:200)';
            d1 = arrayfun(@(k) distances(g1, a(k), b(k)), 1:200)';
            m = t > 0.02;
            tc.verifyLessThan(median(d1(m) ./ t(m) - 1), 0.03);
            tc.verifyGreaterThan(min(d1(m) ./ t(m) - 1), -0.01);  % never shorter than the truth
            tc.verifyLessThan(median(d1(m) ./ t(m)), median(d0(m) ./ t(m)));
            far = t > 0.8 * pi * tc.R;                              % where the heat method saturates
            tc.verifyLessThan(max(abs(d1(far) ./ t(far) - 1)), 0.05);
        end

        function edgegraphIsSymmetric(tc)
            g = rheome.geom.edgegraph(tc.S);  D = distances(g, 1:50, 1:50);
            tc.verifyEqual(D, D', 'AbsTol', 1e-15);
        end

        % ---------------- rheome.geom.tiles ----------------

        function tilesPartitionAndAdjacencyIsTheCutEdges(tc)
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=3);  G = rheome.geom.tiles(T, tc.S, 3);
            tc.verifyEqual(numel(G.node_id), 8);
            tc.verifyEqual(full(sum(G.P, 2)), ones(tc.S.nV, 1));
            tc.verifyEqual(G.A, G.A');
            tc.verifyEqual(full(diag(G.A)), zeros(8, 1));
            F = double(tc.S.Faces);
            E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
            tc.verifyEqual(full(sum(G.A(:))) / 2, sum(G.tileOf(E(:,1)) ~= G.tileOf(E(:,2))));
            tc.verifyTrue(all(sum(G.A > 0, 2) >= 1));             % a sphere has no isolated tile
        end

        function heapArithmeticHoldsForAFullTree(tc)
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=3);  G = rheome.geom.tiles(T, tc.S, 3);
            tc.verifyTrue(G.heap);
            tc.verifyEqual(G.node_id', 8:15);
        end

        function heapArithmeticIsReportedBrokenWhenABranchStopsEarly(tc)
            % ⚠ One early stop and floor(k/2) no longer names the parent.
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=4, MinArea=0.9 * sum(tc.M(:)) / 8);
            T2 = T;  T2.node_id(end) = T2.node_id(end) + 100;
            G = rheome.geom.tiles(T2, tc.S, 2);
            tc.verifyFalse(G.heap);
        end
    end

    methods
        function x = blob(tc, v, sigma)
            d = zeros(tc.S.nV, 1);  d(v) = 1;
            b = tc.basis;
            x = b.Phi * (exp(-b.Lambda(:) * sigma^2 / 2) .* (b.Phi' * (b.Mass * d)));
            x = x / max(x);
        end
        function nb = ring(tc, v)
            F = tc.S.Faces;  f = any(F == v, 2);
            nb = setdiff(unique(F(f, :)), v);
        end
    end
end

% Author: Diellor Basha, 2026
