classdef tTileTrack < matlab.unittest.TestCase
% Tracking at tile resolution: rheome.flow.movingpeak (the plant), rheome.geom.tilemean (the observation),
% rheome.detect.tilepeaks (the peak as a tile) and rheome.detect.tiletrack (the trajectory as a walk on the tile
% graph). What is asserted: the plant moves at the speed asked; the tile mean is exact and rolls up
% the tree; a tile peak is a tile beating its neighbours; a clean mover is followed as ONE track whose
% steps are all to adjacent tiles; and the margin hysteresis stops a boundary peak from flickering.
%
% Author: Diellor Basha, 2026

    properties
        S; L; M; basis; R = 0.07; T; G; a; ge
    end

    methods (TestClassSetup)
        function mesh(tc)
            [V, F] = rheome.geom.icosphere(4);                         % 2562 vertices
            tc.S = struct('Vertices', V * tc.R, 'Faces', F, 'nV', size(V,1), 'nF', size(F,1));
            [tc.L, tc.M] = rheome.operators.laplace_beltrami(tc.S.Vertices, tc.S.Faces);
            tc.basis = rheome.eigen.modes(tc.L, tc.M, 300);
            tc.T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=6);
            tc.G = rheome.geom.tiles(tc.T, tc.S, 6);
            tc.a = full(sum(tc.M, 2));
            tc.ge = rheome.geom.edgegraph(tc.S);
        end
    end

    methods (Test)

        function movingpeakCoversThePathAtTheAskedSpeed(tc)
            [pv, cu] = tc.path(1, 0.06);
            P = rheome.flow.movingpeak(tc.atoms(pv, 0.008), cu, SpeedMS=0.1, SampleRate=200);
            tc.verifyEqual(P.s(end), cu(end), 'AbsTol', 0.1 / 200);
            tc.verifyEqual(diff(P.s), 0.1 / 200 * ones(numel(P.s) - 1, 1), 'AbsTol', 1e-12);
            tc.verifyEqual(P.k(1), 1);  tc.verifyEqual(P.k(end), numel(pv));
        end

        function movingpeakAtZeroSpeedHolds(tc)
            [pv, cu] = tc.path(1, 0.04);
            P = rheome.flow.movingpeak(tc.atoms(pv, 0.008), cu, SpeedMS=0, SampleRate=100, StillSec=0.5);
            tc.verifySize(P.X, [tc.S.nV 50]);
            tc.verifyEqual(P.X(:,1), P.X(:,end));
        end

        function tilemeanIsTheAreaWeightedMeanAndRollsUp(tc)
            X = randn(tc.S.nV, 3);
            Y = rheome.geom.tilemean(tc.G, X, tc.a);
            k = 5;  m = tc.G.tileOf == k;
            tc.verifyEqual(Y(k,:), sum(tc.a(m) .* X(m,:)) / sum(tc.a(m)), 'AbsTol', 1e-12);
            G5 = rheome.geom.tiles(tc.T, tc.S, 5);  Y5 = rheome.geom.tilemean(G5, X, tc.a);
            % a parent is the area-weighted mean of its two children
            aT = full(double(tc.G.P)' * tc.a);
            for i = 1:3
                kids = find(floor(tc.G.node_id / 2) == G5.node_id(i));
                tc.verifyEqual(Y5(i,:), sum(aT(kids) .* Y(kids,:)) / sum(aT(kids)), 'AbsTol', 1e-12);
            end
            tc.verifyTrue(tc.G.heap);
        end

        function tilepeaksBeatEveryNeighbour(tc)
            Y = randn(numel(tc.G.node_id), 20);
            P = rheome.detect.tilepeaks(Y, tc.G);
            [I, J] = find(tc.G.A);
            for t = 1:20
                for k = P.tiles{t}'
                    nb = J(I == k);
                    tc.verifyTrue(all(Y(k,t) > Y(nb,t) | (Y(k,t) == Y(nb,t) & k < nb)));
                end
            end
            tc.verifyGreaterThan(sum(P.count), 20);            % a random field has several
        end

        function aCleanMoverIsOneTrackOfAdjacentSteps(tc)
            [pv, cu] = tc.path(7, 0.09);
            P = rheome.flow.movingpeak(tc.atoms(pv, 0.012), cu, SpeedMS=0.1, SampleRate=100);
            Y = rheome.geom.tilemean(tc.G, P.X, tc.a);
            Pk = rheome.detect.tilepeaks(Y, tc.G, MinHeight=0.5, Relative=true);
            out = rheome.detect.tiletrack(Pk, tc.G, Y, LostFrames=10, SampleRate=100);
            tc.verifyEqual(out.nTracks, 1);
            x = out.tracks(1);
            tc.verifyGreaterThanOrEqual(x.hops, 2, 'the peak should cross tiles');
            st = [x.tiles(1:end-1) x.tiles(2:end)];  st = st(st(:,1) ~= st(:,2), :);
            tc.verifyTrue(all(full(tc.G.A(sub2ind(size(tc.G.A), st(:,1), st(:,2)))) > 0), ...
                'every step must be to an adjacent tile');
        end

        function marginStopsBoundaryFlicker(tc)
            % two adjacent tiles trading the lead by a hair every frame: with no margin the track hops
            % every frame; with a 10% margin it holds.
            [I, J] = find(tc.G.A);  k1 = I(1);  k2 = J(1);
            nT = 40;  Y = zeros(numel(tc.G.node_id), nT);
            Y(k1,:) = 1 + 0.01 * (-1).^(1:nT);  Y(k2,:) = 1 - 0.01 * (-1).^(1:nT);
            Pk = rheome.detect.tilepeaks(Y, tc.G, MinHeight=0.5, Relative=true);
            o0 = rheome.detect.tiletrack(Pk, tc.G, Y, Margin=0);
            o1 = rheome.detect.tiletrack(Pk, tc.G, Y, Margin=0.1);
            tc.verifyGreaterThan(max([o0.tracks.hops]), 10);
            tc.verifyEqual(o1.nTracks, 1);
            tc.verifyEqual(o1.tracks(1).hops, 0);
        end

        function aJumpToANonAdjacentTileStartsANewTrack(tc)
            nTile = numel(tc.G.node_id);  far = find(~tc.G.A(1,:) & (1:nTile) ~= 1, 1);
            Y = zeros(nTile, 20);  Y(1, 1:10) = 1;  Y(far, 11:20) = 1;
            Pk = rheome.detect.tilepeaks(Y, tc.G, MinHeight=0.5, Relative=true);
            out = rheome.detect.tiletrack(Pk, tc.G, Y, LostFrames=0);
            tc.verifyEqual(out.nTracks, 2);
        end
    end

    methods
        function [pv, cu] = path(tc, s0, len)
            d = distances(tc.ge, s0);  [~, tgt] = min(abs(d - len));
            [pv, ~, eg] = shortestpath(tc.ge, s0, tgt);
            cu = [0; cumsum(tc.ge.Edges.Weight(eg))];
        end
        function A = atoms(tc, pv, sigma)
            b = tc.basis;  A = zeros(tc.S.nV, numel(pv));
            g = rheome.filters.mexhat(b.Lambda(:), sigma^2 / 2);
            for k = 1:numel(pv)
                d = zeros(tc.S.nV, 1);  d(pv(k)) = 1;
                x = b.Phi * (g .* (b.Phi' * (b.Mass * d)));  A(:,k) = x / max(x);
            end
        end
    end
end

% Author: Diellor Basha, 2026
