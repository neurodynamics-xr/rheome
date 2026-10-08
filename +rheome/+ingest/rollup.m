function T = rollup(T0, g)
% INGEST.ROLLUP  Every coarser level from level 0, by pairing rows. No signal needed.
%
%   T = rheome.ingest.rollup(T0, g)
%
% Level L+1 row k is rows 2k-1 and 2k of level L: sums and counts add, maxima take the
% max, minima the min. A missing second child (odd row count) is the neutral element --
% 0 for sums, counts, absMax and envMax; -Inf for max; +Inf for min -- which is exactly
% what an empty frame already holds, so partial and absent children are the same case.
%
% ⚠ This is the property the whole store rests on: a parent equals its children's merge
% EXACTLY, so a query can prune a subtree from a bound and a rebuild at any level equals a
% direct computation at that level (tests/tIngestRollup.m pins this to round-off for sums
% and exactly for extrema). Nothing non-mergeable may be added to T0 (design §2).
%
% INPUTS:
%   T0  level-0 statistics from rheome.ingest.reduce
%   g   grid from rheome.ingest.grid
% OUTPUT:
%   T   {1 x Lmax+1} cell; T{L+1} has the fields of T0 at level L with K(L+1) rows
%
% See also: rheome.ingest.reduce, rheome.ingest.grid
%
% Author: Diellor Basha, 2026

    arguments
        T0 (1,1) struct
        g  (1,1) struct
    end

    T = cell(1, g.Lmax + 1);
    T{1} = T0;
    for L = 1:g.Lmax
        P = T{L};  K = g.K(L+1);
        Q = struct();
        Q.n      = i_pair(P.n,      K, @plus, 0);
        Q.sumX   = i_pair(P.sumX,   K, @plus, 0);
        Q.sumX2  = i_pair(P.sumX2,  K, @plus, 0);
        for f = ["sumAbs","sumSqrt","sumX3","sumX4"]            % the shape moments: sums too
            if isfield(P, f), Q.(f) = i_pair(P.(f), K, @plus, 0); end
        end
        Q.absMax = i_pair(P.absMax, K, @max,  0);
        Q.min    = i_pair(P.min,    K, @min,  +Inf);
        Q.max    = i_pair(P.max,    K, @max,  -Inf);
        Q.energy = i_pair(P.energy, K, @plus, 0);
        Q.envMax = i_pair(P.envMax, K, @max,  0);
        Q.nCoi   = i_pair(P.nCoi,   K, @plus, 0);
        if isfield(P, 'senergy')                                   % the wavelength axis, when built
            Q.senergy = i_pair(P.senergy, K, @plus, 0);
            Q.senvMax = i_pair(P.senvMax, K, @max,  0);
        end
        Q.level  = L;
        Q.K      = K;
        T{L+1} = Q;
    end
end

function B = i_pair(A, K, op, neutral)
% Combine rows (2k-1, 2k) of A along dim 1 into row k; pad to 2K rows with the neutral.
    Kc = size(A, 1);
    if 2*K > Kc
        A(Kc+1:2*K, :, :) = cast(neutral, 'like', A);
    end
    B = op(A(1:2:end, :, :), A(2:2:end, :, :));
end
% Author: Diellor Basha, 2026
