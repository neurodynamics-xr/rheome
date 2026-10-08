function y = bound(x, dir)
% INGEST.BOUND  Cast an extremum to single ROUNDING OUTWARD, so it stays a bound.
%
%   y = rheome.ingest.bound(x, 'down')     y <= x, in single
%   y = rheome.ingest.bound(x, 'up')       y >= x, in single
%
% ⚠ WHY ROUND-TO-NEAREST IS WRONG HERE. single(x) can land above a minimum or below a
% maximum by one ulp, and a stored extremum that is not a bound is not an extremum: the
% min/max envelope drawn from it excludes samples it covers, and a pruning descent that
% trusts "parent >= child" can drop a subtree that should have survived. Measured on the
% test fixture before this existed: every tile at every level failed to contain its samples,
% and the raw trace fell outside the drawn band.
%
% The slip is one ulp, 6e-8 relative, so no number visibly changes -- but the GUARANTEE
% does, and the guarantee is what makes the store readable at a glance. Monotone rounding
% also keeps the roll-up exact: max over children of up(x) equals up of max over children.
%
% Infinities (a fully masked tile carries +/-Inf) are left alone: they are already exact,
% and eps(Inf) is Inf, so nudging would give NaN.
%
% See also: rheome.ingest.reduce, rheome.ingest.rollup
%
% Author: Diellor Basha, 2026

    y = single(x);
    d = double(y);
    switch dir
        case 'down'
            i = d > x & isfinite(y);
            y(i) = y(i) - eps(y(i));
        case 'up'
            i = d < x & isfinite(y);
            y(i) = y(i) + eps(y(i));
        otherwise
            error('ingest:bound:dir', 'dir must be ''down'' or ''up''.');
    end
end
% Author: Diellor Basha, 2026
