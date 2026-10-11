function [C, lab, info] = cluster(X, K, opts)
% MODES.CLUSTER  Phase-invariant complex k-means: complex spatial modes from amplitude-normalised snapshots.
%
%   [C, lab] = rheome.modes.cluster(X, K)               % X [d x n] complex snapshots (columns)
%   [C, lab, info] = rheome.modes.cluster(X, K, Replicates=5, MaxIter=200, Seed=1)
%
% Each snapshot is normalised to unit length and assigned to the centroid with the largest |<x, c_k>|^2, which
% does not depend on the snapshot's global phase. The centroid that maximises sum |<x, c>|^2 over its members
% is the principal eigenvector of their scatter sum x x' -- the complex analogue of the polarity-invariant
% ("modified") k-means of EEG microstates. |c_k| is the mode's amplitude topography and arg c_k its phase
% layout, defined up to one global phase, fixed here so the largest-amplitude element is real and positive.
% k-means++ seeding on 1 - |<x, c>|^2; the best of Replicates runs (largest explained share) is returned.
%
%   C     [d x K] unit-norm centroids, ordered by decreasing membership
%   lab   [1 x n] labels
%   info  .gev  sum_n max_k |<x_n, c_k>|^2 / n (explained share, 1 = every snapshot equals a centroid)
%         .iter  iterations of the returned run   .converged
%
% Author: Diellor Basha, 2026

    arguments
        X {mustBeNumeric}
        K (1,1) double {mustBeInteger, mustBePositive}
        opts.Replicates (1,1) double {mustBeInteger, mustBePositive} = 5
        opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.Seed (1,1) double = 1
    end
    X = double(X);  X = X ./ vecnorm(X);  X(:, ~all(isfinite(X), 1)) = [];
    [d, n] = size(X);
    if n < K, error('modes:cluster:few', '%d snapshots for K = %d', n, K); end
    s = RandStream('mt19937ar', 'Seed', opts.Seed);
    best = -Inf;
    for rep = 1:opts.Replicates
        % k-means++ on 1 - |<x, c>|^2
        C = zeros(d, K);  C(:, 1) = X(:, randi(s, n));  D = 1 - abs(C(:, 1)' * X).^2;
        for k = 2:K
            p = max(D, 0);  p = cumsum(p) / sum(p);
            C(:, k) = X(:, find(rand(s) <= p, 1));  D = min(D, 1 - abs(C(:, k)' * X).^2);
        end
        lab = zeros(1, n);  conv = false;
        for it = 1:opts.MaxIter
            [~, nl] = max(abs(C' * X).^2, [], 1);
            if isequal(nl, lab), conv = true;  break; end
            lab = nl;
            for k = 1:K
                Xk = X(:, lab == k);
                if isempty(Xk), [~, j] = min(max(abs(C' * X).^2, [], 1));  C(:, k) = X(:, j);  continue; end
                [V, L] = eig(Xk * Xk', 'vector');  [~, j] = max(real(L));  C(:, k) = V(:, j);
            end
        end
        g = mean(max(abs(C' * X).^2, [], 1));
        if g > best, best = g;  Cb = C;  lb = lab;  itb = it;  cb = conv; end
    end
    cnt = accumarray(lb(:), 1, [K 1]);  [~, o] = sort(cnt, 'descend');
    C = Cb(:, o);  [~, inv] = sort(o);  lab = inv(lb).';
    lab = reshape(lab, 1, []);
    for k = 1:K, [~, j] = max(abs(C(:, k)));  C(:, k) = C(:, k) * exp(-1i * angle(C(j, k))); end
    info = struct('gev', best, 'iter', itb, 'converged', cb);
end

% Author: Diellor Basha, 2026
