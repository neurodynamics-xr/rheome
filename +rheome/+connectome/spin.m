function perm = spin(C, hemi, nPerm, seed)
% CONNECTOME.SPIN  Parcel-level spin permutations (Alexander-Bloch 2018; Vasa 2018, Hungarian assignment).
%
%   perm = rheome.connectome.spin(C, hemi, nPerm)
%   perm = rheome.connectome.spin(C, hemi, nPerm, seed)
%
% Rotates the parcel centroids on the sphere at random and reassigns each parcel to a rotated
% one, so a null map x(perm(:,k)) keeps the spatial autocorrelation of x. The left hemisphere
% gets a random rotation R, the right the mirrored F*R*F (F flips x), as in Alexander-Bloch.
% Each hemisphere is assigned one-to-one by the optimal (Hungarian) matching of original to rotated
% centroids (MATLAB matchpairs), so every null map is a permutation within each hemisphere.
%
% INPUTS:
%   C      [n x 3] parcel centroids on the SPHERE surface (e.g. ?h.sphere), each hemisphere centred at 0
%   hemi   [n x 1] logical, true = right hemisphere
%   nPerm  number of rotations (the cards ask 10 000)
%   seed   rng seed (default 0) -- the null is reproducible
%
% OUTPUT:
%   perm [n x nPerm] parcel indices: the k-th null map of x is x(perm(:,k))
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(seed), seed = 0; end
    hemi = logical(hemi(:));
    n = size(C, 1);
    assert(size(C, 2) == 3 && numel(hemi) == n, 'rheome:spin:size', 'C must be n x 3 and hemi n x 1.');
    stream = RandStream('mt19937ar', 'Seed', seed);
    F = diag([-1 1 1]);
    perm = zeros(n, nPerm);
    for k = 1:nPerm
        [Q, R] = qr(randn(stream, 3));
        Q = Q * diag(sign(diag(R)));                    % Haar-uniform orthogonal
        if det(Q) < 0, Q(:, 1) = -Q(:, 1); end          % a rotation, not a reflection
        for right = [false true]
            idx = find(hemi == right);
            if isempty(idx), continue; end
            Rh = Q;  if right, Rh = F * Q * F; end
            rotated = C(idx, :) * Rh';
            cost = sqrt(sum((permute(C(idx, :), [1 3 2]) - permute(rotated, [3 1 2])).^2, 3));
            M = matchpairs(cost, 1e3 * max(cost(:)) + 1);  % unmatched cost high: a full assignment
            perm(idx(M(:, 1)), k) = idx(M(:, 2));
        end
    end
end

% Author: Diellor Basha, 2026
