function [idx, w] = shepard_weights(X, Y, neighbours, power)
% REGISTER.SHEPARD_WEIGHTS  Brainstorm's Shepard weights of the points Y from the points X (unit-sphere coordinates).
%
%   [idx, w] = rheome.register.shepard_weights(X, Y)                  % 8 neighbours, power 2
%   [idx, w] = rheome.register.shepard_weights(X, Y, neighbours, power)
%
% X [n x 3] the source points, Y [m x 3] the points to place.
% idx [m x p] the source index of each weight (1-based), w [m x p] the weights, p = min(neighbours, n) - 1 (1 for a single
%     neighbour); each row non-negative and summing to 1.
%
% As nxr.register.sphere.shepard_weights, a port of Brainstorm's bst_shepards (excludeParam = 0, as tess_interp_tess2tess
% calls it): the k source points nearest each point (Euclidean); d2 their squared distances (0 becomes eps); R2 the k-th one's;
% w_i = ((R2 - d2_i) / (R2 d2_i))^power for the first k - 1 neighbours only (the k-th has weight 0 by construction), normalised
% to sum 1; a row summing to 0 (the point overlaps a vertex) takes weight 1 on its nearest; k = 1 is the nearest vertex.
%
% See also: rheome.register.registration
%
% Author: Diellor Basha, 2026
arguments
    X (:,3) double, Y (:,3) double
    neighbours (1,1) double {mustBeInteger, mustBePositive} = 8
    power (1,1) double = 2
end
k = min(neighbours, size(X, 1));
[idx, dist] = knnsearch(X, Y, 'K', k);
if k == 1, w = ones(size(Y, 1), 1); return, end
d2 = dist .^ 2; d2(d2 == 0) = eps;
W = ((d2(:, k) - d2) ./ (d2(:, k) .* d2)) .^ power;
W = W(:, 1:k-1); idx = idx(:, 1:k-1);
s = sum(W, 2); zero = s == 0;
W(zero, :) = 0; W(zero, 1) = 1; s(zero) = 1;
w = W ./ s;
end
