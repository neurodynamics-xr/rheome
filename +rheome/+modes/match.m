function [r, perm, R] = match(A, B)
% MODES.MATCH  Hungarian matching of two sets of complex modes by phase-invariant |r|.
%
%   [r, perm] = rheome.modes.match(A, B)     % A [d x K], B [d x K2], K <= K2
%
% R = |A' B| after unit-normalising every column (global phase does not matter). perm(k) is the column of B
% matched to A(:, k) (MATLAB matchpairs on 1 - R, every A column matched), r(k) = R(k, perm(k)).
%
% Author: Diellor Basha, 2026

    A = A ./ vecnorm(A);  B = B ./ vecnorm(B);
    R = abs(A' * B);
    M = matchpairs(1 - R, 1e6);                    % cost of leaving unmatched never chosen
    perm = zeros(size(A, 2), 1);  perm(M(:, 1)) = M(:, 2);
    r = R(sub2ind(size(R), (1:size(A, 2))', perm));
end

% Author: Diellor Basha, 2026
