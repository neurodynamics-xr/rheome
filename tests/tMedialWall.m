classdef tMedialWall < matlab.unittest.TestCase
% The medial wall (decision 3882c4ce, 2026-10-14): excluded from measurements, harmonically filled for operators.
% rheome.geom.cortexmask, rheome.operators.harmonicfill, rheome.geom.tiles/.tilemean (cortex area),
% rheome.group.spin (Cortex=), rheome.detect.criticalPoints (.wall).
%
% Author: Diellor Basha, 2026

    properties
        S       % icosphere, radius 0.1 m, with a polar cap marked as wall
        wall    % [nV x 1] logical
    end

    methods (TestClassSetup)
        function sphere(tc)
            [V, F] = rheome.geom.icosphere(4);  V = 0.1 * V ./ vecnorm(V, 2, 2);
            tc.wall = V(:, 3) > 0.08;                                % a cap of ~10 % of the area
            tc.S = struct('Vertices', V, 'Faces', F, 'nV', size(V, 1), 'Cortex', ~tc.wall, ...
                          'VertNormals', V ./ vecnorm(V, 2, 2));
        end
    end

    methods (Test)
        function maskByPrefix(tc)
            % 'unknown' and 'corpuscallosum' are the wall; 'medialorbitofrontal' is cortex
            A = struct('Label', {{'unknown L', 'medialorbitofrontal L', 'precentral L', 'corpuscallosum L'}}, ...
                       'Membership', sparse([1 2 3 4], [1 2 3 4], true, 4, 6));
            c = rheome.geom.cortexmask(A, 6);
            tc.verifyEqual(c, [false; true; true; false; false; false]);
            tc.verifyEqual(rheome.geom.cortexmask(rmfield(tc.S, 'Cortex')), true(tc.S.nV, 1));
        end

        function fillIsHarmonic(tc)
            S = tc.S;  rng(1);  X = randn(S.nV, 4);
            Y = rheome.operators.harmonicfill(S, X, S.Cortex);
            tc.verifyEqual(Y(S.Cortex, :), X(S.Cortex, :));          % cortex untouched
            L = rheome.operators.laplace_beltrami(S.Vertices, S.Faces);
            r = L(tc.wall, :) * Y;                                   % Laplace's equation on the wall
            tc.verifyLessThan(max(abs(r), [], 'all') / max(abs(L * X), [], 'all'), 1e-10);
            k = rheome.operators.harmonicfill(S, ones(S.nV, 1), S.Cortex);
            tc.verifyEqual(k, ones(S.nV, 1), 'AbsTol', 1e-12);       % a constant stays constant (max principle)
            X3 = randn(3 * S.nV, 2);  Y3 = rheome.operators.harmonicfill(S, X3, S.Cortex);
            tc.verifyEqual(Y3(2:3:end, :), rheome.operators.harmonicfill(S, X3(2:3:end, :), S.Cortex), 'AbsTol', 1e-12);
        end

        function islandErrors(tc)
            tc.verifyError(@() rheome.operators.harmonicfill(tc.S, ones(tc.S.nV, 1), false(tc.S.nV, 1)), ...
                'operators:harmonicfill:island');
        end

        function tileCortexArea(tc)
            S = tc.S;  T = rheome.geom.tree(S, MaxDepth=3);  G = rheome.geom.tiles(T, S, 3);
            tc.verifyEqual(sum(G.cortex_area), sum(G.area(S.Cortex)), 'RelTol', 1e-12);
            tc.verifyTrue(all(G.cortex_fraction >= 0 & G.cortex_fraction <= 1));
            tc.verifyLessThan(min(G.cortex_fraction), 1);           % some tile holds wall
            X = randn(S.nV, 1);  Xw = X;  Xw(tc.wall) = 1e6;         % wall values never reach a tile mean
            tc.verifyEqual(rheome.geom.tilemean(G, Xw, G.area), rheome.geom.tilemean(G, X, G.area), 'AbsTol', 1e-9);
            m = rheome.geom.tilemean(G, Xw, [], 'median');  tc.verifyLessThan(max(m), 1e6);
            G2 = rheome.geom.tiles(T, S, 2);  [~, par] = ismember(floor(G.node_id / 2), G2.node_id);
            tc.verifyEqual(accumarray(par, G.cortex_area), G2.cortex_area, 'RelTol', 1e-12);   % roll-up exact
        end

        function spinDropsWallTiles(tc)
            rng(2);  m = 12;  U = randn(m, 3);  U = U ./ vecnorm(U, 2, 2);  h = repmat("L", m, 1);
            M = randn(m, 1);  Ref = M + 0.1 * randn(m, 1);  k = true(m, 1);  k([3 7]) = false;
            R = rheome.group.spin(M, Ref, U, h, NSpin=50, Cortex=k);
            R0 = rheome.group.spin(M(k), Ref(k), U(k, :), h(k), NSpin=50);
            tc.verifyEqual(R.r, R0.r);  tc.verifyEqual(R.p, R0.p);  tc.verifyEqual(R.cortex, k);
        end

        function wallCoresFlagged(tc)
            % a rotation about z: index +1 at both poles; the north pole sits in the wall cap
            S = tc.S;  n = S.VertNormals;  J3 = cross(repmat([0 0 1], S.nV, 1), S.Vertices, 2);
            J3 = J3 - sum(J3 .* n, 2) .* n;  J = reshape(J3', [], 1);
            cp = rheome.detect.criticalPoints(J, S);
            tc.verifyEqual(sum(cp.chi), 2);                          % Poincare-Hopf still over every face
            tc.verifyTrue(all(cp.pos(:, 3) < 0.08));                 % nothing reported on the wall
            tc.verifyTrue(any(cp.wall.pos(:, 3) > 0.08));            % the north core is flagged
        end
    end
end
