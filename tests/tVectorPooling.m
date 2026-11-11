classdef tVectorPooling < matlab.unittest.TestCase
    % Tangent vectors summed over a tile or at a place on a FOLDED surface are carried to the anchor first, never added in 3-D (nxr 1.15.0,
    % tests/test_vector_pooling.py). A sheet corrugated across x (2.5 mm ridges, 8 mm apart) carries a field whose vectors are parallel along the
    % surface; in 3-D they point up one wall of a ridge and down the other, so the raw sum shortens; carried to the centre it is 1.
    %
    % Author: Diellor Basha, 2026

    properties
        V, F, U, Nrm, Kit
    end

    methods (TestClassSetup)
        function sheet(t)
            n = 51;  h = 0.5;  a = 2.5;  wl = 8;
            [i, j] = ndgrid(0:n-1, 0:n-1);  x = i(:) * h;  y = j(:) * h;  w = 2*pi/wl;
            H = a * sin(w * x);  dH = a * w * cos(w * x);
            q = reshape(i(1:end-1, 1:end-1) * n + j(1:end-1, 1:end-1), [], 1);
            t.F = [[q q+n q+1]; [q+1 q+n q+n+1]] + 1;
            t.V = [x y H];
            u = [ones(size(x)) zeros(size(x)) dH];  t.U = u ./ vecnorm(u, 2, 2);
            nr = cross(t.U, repmat([0 1 0], numel(x), 1), 2);  nr = nr .* sign(nr(:, 3));  t.Nrm = nr;
            t.Kit = rheome.geom.GeometryKit(t.V, t.F);
        end
    end

    methods (Test)
        function liftAndEmbedAreInverse(t)
            v = (1:97:size(t.V, 1))';
            w = t.Kit.lift(t.U(v, :), v);
            X = t.Kit.embedAtVertices(v, w);
            t.verifyLessThan(max(vecnorm(X - t.U(v, :), 2, 2)), 0.05);          % the unit tangent, through the face-by-face gauge
            t.verifyEqual(abs(w), ones(size(w)), 'AbsTol', 1e-9);
        end

        function parallelFieldAddsToOneWhenCarried(t)
            centre = find(vecnorm(t.V(:, 1:2) - max(t.V(:, 1:2), [], 1) / 2, 2, 2) == min(vecnorm(t.V(:, 1:2) - max(t.V(:, 1:2), [], 1) / 2, 2, 2)), 1);
            tile = find(vecnorm(t.V(:, 1:2) - t.V(centre, 1:2), 2, 2) < 9);
            a = zeros(size(t.V, 1), 1);  a(tile) = centre;
            pool = rheome.cores.Pooling(t.Kit, a);
            raw = norm(sum(t.U(tile, :), 1)) / numel(tile);
            w = sum(pool.carry(t.U(tile, :), tile)) / numel(tile);
            c = pool.embed(w, centre);
            t.verifyLessThan(raw, 0.8);  t.verifyEqual(abs(w), 1, 'AbsTol', 0.03);
            t.verifyLessThan(acosd(min(1, c * t.U(centre, :)' / norm(c))), 3);
            t.verifyLessThan(abs(c * t.Nrm(centre, :)'), 0.05);
            fprintf('parallel field over a ridge: raw 3-D sum %.3f, carried %.3f\n', raw, abs(w));
        end

        function tileCrossingsAndTheirRAreOneWhereTheRawSumsWereNot(t)
            nv = size(t.V, 1);  rng(0);
            lab = double(t.V(:, 1) > max(t.V(:, 1)) / 2) * 2 + double(t.V(:, 2) > max(t.V(:, 2)) / 2);
            sup = (1:nv)';
            pool = rheome.cores.tilepool(t.Kit, t.V, lab, sup);
            vs = randi(nv, 400, 1);  S = zeros(800, 9);
            for k = 1:400
                S(2*k-1, :) = [1 t.V(vs(k), :) 3 * t.U(vs(k), :) k 0];  S(2*k, :) = [1 t.V(vs(k), :) 3 * t.U(vs(k), :) k 1];
            end
            [s, n, key, vec] = rheome.cores.crossings(S, 0:1, lab, t.V, t.Nrm, sup, pool);
            R = rheome.cores.pathR(s, n);
            Rc = rheome.cores.pathR(accumarray(key, real(vec), [8 1]) + 1i * accumarray(key, imag(vec), [8 1]), n);
            % the old computation: raw 3-D unit vectors of a track added per tile
            near = knnsearch(t.V, S(:, 2:4));  per = containers.Map('KeyType', 'double', 'ValueType', 'any');
            [~, ~, tr] = unique(S(:, 8));  acc = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for r = 1:size(S, 1)
                kk = sprintf('%d_%d', lab(near(r)), tr(r));
                if isKey(acc, kk), acc(kk) = acc(kk) + t.U(near(r), :); else, acc(kk) = t.U(near(r), :); end
            end
            sumsT = zeros(4, 3);  nT = zeros(4, 1);
            for kk = keys(acc)
                p = sscanf(kk{1}, '%d_%d');  m = acc(kk{1});  sumsT(p(1) + 1, :) = sumsT(p(1) + 1, :) + m / norm(m);  nT(p(1) + 1) = nT(p(1) + 1) + 1;
            end
            old = sum(vecnorm(sumsT, 2, 2)) / sum(nT);
            t.verifyLessThan(old, 0.8);  t.verifyGreaterThan(R, 0.97);
            t.verifyEqual(Rc, R, 'AbsTol', 1e-6);                              % the complex sums give the same R as the 3-D ones
            fprintf('tile crossings R: raw 3-D %.3f, carried %.3f\n', old, R);
            for q = 1:4, t.verifyLessThan(abs(s(2*q-1, :) * t.Nrm(pool.Centre(q), :)'), 1e-3 * max(1, norm(s(2*q-1, :)))); end   % tangent vectors at the centres
        end

        function directionConsistencyCarriesToThePlace(t)
            nv = size(t.V, 1);  rng(1);  sup = (1:nv)';
            pl = zeros(4, 1);  c = [6 6; 6 19; 19 6; 19 19];
            for q = 1:4, [~, pl(q)] = min(vecnorm(t.V(:, 1:2) - c(q, :), 2, 2)); end
            vs = randi(nv, 1500, 1);  S = [ones(1500, 1) t.V(vs, :) 3 * t.U(vs, :) zeros(1500, 2)];
            pool = rheome.cores.placepool(t.Kit, t.V, sup, pl);
            R = rheome.cores.placeconsistency(S, pool, sup, t.V, t.Nrm, pl, 5, 10);
            [~, of] = min(sum((permute(t.V, [1 3 2]) - permute(t.V(pl, :), [3 1 2])).^2, 3), [], 2);
            [~, nn] = min(sum((permute(S(:, 2:4), [1 3 2]) - permute(t.V, [3 1 2])).^2, 3), [], 2);
            old = 0;  po = of(nn);
            for p = 1:4, old = old + norm(sum(S(po == p, 5:7) / 3, 1)); end
            old = old / size(S, 1);
            t.verifyLessThan(old, 0.8);  t.verifyGreaterThan(R, 0.95);
            fprintf('direction consistency: raw 3-D %.3f, carried %.3f\n', old, R);
        end
    end
end
