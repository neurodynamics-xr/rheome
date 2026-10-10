classdef tConnectomeGroup < matlab.unittest.TestCase
% TCONNECTOMEGROUP  Group connectome tools for the PREVENT-AD cards (R-T1): the Betzel consensus keeps
% each distance bin's mean edge count, the spin null permutes within a hemisphere, matched controls
% obey the length and decile rules, and the parcel geodesic matches the great-circle distance on a sphere.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function consensusOfIdenticalSubjectsIsTheGraph(tc)
            rng(1);  n = 20;  W = rand(n) .* (rand(n) < 0.3);  W = triu(W, 1);  W = W + W';
            A = repmat(W, 1, 1, 5);
            D = squareform(pdist(rand(n, 3)));  hemi = (1:n)' > n/2;
            [G, Gw] = rheome.connectome.consensus(A, D, hemi, 5);
            tc.verifyEqual(G, W > 0);
            tc.verifyEqual(Gw, W, 'AbsTol', 1e-12);
        end

        function consensusKeepsTheMeanCountPerBin(tc)
            rng(2);  n = 30;  S = 12;  nBins = 4;
            D = squareform(pdist(rand(n, 3)));  hemi = (1:n)' > n/2;
            A = zeros(n, n, S);
            for s = 1:S, M = triu(rand(n) < 0.25, 1);  A(:, :, s) = M + M'; end
            G = rheome.connectome.consensus(A, D, hemi, nBins);
            tc.verifyEqual(G, G');
            tc.verifyFalse(any(diag(G)));
            upper = triu(true(n), 1);  same = hemi == hemi';
            for cls = [true false]
                idx = find(upper & (same == cls));
                d = D(idx);  bin = discretize(d, linspace(min(d), max(d) + eps(max(d)), nBins + 1));
                for b = 1:nBins
                    inBin = idx(bin == b);  B = reshape(A > 0, n*n, S);
                    tc.verifyEqual(nnz(G(inBin)), round(mean(sum(B(inBin, :), 1))));
                end
            end
        end

        function spinPermutesWithinHemisphere(tc)
            [V, ~] = rheome.geom.icosphere(1);  V = V ./ vecnorm(V, 2, 2);
            n = size(V, 1);  C = [V; V];  hemi = [false(n, 1); true(n, 1)];
            P = rheome.connectome.spin(C, hemi, 50, 7);
            tc.verifyEqual(sort(P(1:n, :)), repmat((1:n)', 1, 50));
            tc.verifyEqual(sort(P(n+1:end, :)), repmat((n+1:2*n)', 1, 50));
            tc.verifyEqual(P, rheome.connectome.spin(C, hemi, 50, 7));           % seeded
            tc.verifyLessThan(mean(P(:) == repmat((1:2*n)', 50, 1)), 0.5);   % it actually rotates
        end

        function controlsObeyLengthAndDecile(tc)
            rng(3);  n = 40;  N = round(100 * rand(n));  N = triu(N, 1);  N = N + N';
            L = 20 + 100 * rand(n);  L = triu(L, 1);  L = L + L';
            [i, j] = find(triu(true(n), 1));  pairs = [i j];
            targets = pairs(1:10, :);  pool = pairs(11:end, :);
            [sets, nc] = rheome.connectome.matchcontrols(L, N, targets, pool, 100);
            q = quantile(N(triu(true(n), 1) & N > 0), 0.1:0.1:0.9);  dec = @(x) 1 + sum(x(:) > q, 2);
            at = @(M, P) M(sub2ind(size(M), P(:, 1), P(:, 2)));
            for t = find(~isnan(sets(:, 1)))'
                c = pool(sets(t, :), :);
                tc.verifyLessThanOrEqual(abs(at(L, c) - at(L, targets(t, :))), 0.1 * at(L, targets(t, :)) + 1e-12);
                tc.verifyEqual(dec(at(N, c)), repmat(dec(at(N, targets(t, :))), 100, 1));
            end
            tc.verifyEqual(isnan(sets(:, 1)), nc == 0 | at(N, targets) == 0);
        end

        function geodesicIsGreatCircleOnASphere(tc)
            [V, F] = rheome.geom.icosphere(4);  V = V ./ vecnorm(V, 2, 2);  n = size(V, 1);
            S = struct('Vertices', [V; V + [5 0 0]], 'Faces', [F; F + n]);
            north = V(:, 3) > 0.95;  south = V(:, 3) < -0.95;
            member = sparse(logical([north' false(1, n); south' false(1, n); false(1, n) true(1, n)]));
            d = rheome.connectome.geodesicfrom(S, member, 1);
            nv = find(north);  [~, k] = min(sum((V(nv, :) - mean(V(nv, :), 1)).^2, 2));
            exact = mean(acos(min(1, V(south, :) * V(nv(k), :)')));      % great-circle mean over the south cap
            tc.verifyEqual(d(2), exact, 'RelTol', 0.02);
            tc.verifyLessThan(d(1), 0.2);
            tc.verifyTrue(isnan(d(3)));                                    % the other component
        end
    end
end

% Author: Diellor Basha, 2026
