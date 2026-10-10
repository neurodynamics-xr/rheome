function O = occupancy(b, dt, opts)
% DETECT.OCCUPANCY  A dyadic pyramid of on/off occupancy: how long a state lasts, at every scale.
%
%   O = rheome.detect.occupancy(on, dt)                 % on [nCh x n] logical at the finest tile length dt (s)
%   O = rheome.detect.occupancy(on, dt, Levels=10, Theta=0.75)
%
% Level 0 holds the binary state of each finest time tile. Level L holds, for each tile of 2^L finest
% tiles, its OCCUPANCY: the fraction of them that are on. Occupancy is a mean, so it rolls up EXACTLY
% (a parent is the mean of its two children), and a tile is ON at its level when occupancy >= Theta.
%
% ⭐ WHY A PYRAMID AND NOT ONE DETECTOR. An episode lasting tens of seconds with brief dips is, to a
% single-resolution detector, a string of sub-second bursts -- which is how the alpha envelope first
% measured 0.3-1.3 s. In the pyramid the dips stay visible at level 0 and are absorbed at the levels
% where they no longer matter, so duration becomes a function of SCALE rather than of one threshold:
%   .runs{L}   maximal runs of consecutive ON tiles at level L, in seconds -- episode durations at
%              that resolution
%   .depth     for every finest tile that is on, the LONGEST enclosing dyadic tile that is still on
%              (in seconds): how long the state lasts around that moment, tolerant to short gaps
%
% ⚠ Tiles are ALIGNED to the dyadic grid, so an episode straddling a boundary is held by two tiles of
% the level below. .depth therefore understates a straddling episode by up to one level; .runs, which
% concatenates adjacent on-tiles, does not.
%
% INPUTS
%   b       [nCh x n] logical, finest-level state (e.g. rheome.detect.onstate(...).on per channel)
%   dt      finest tile length (s)
%   Levels  number of coarser levels (default: as many as fit)
%   Theta   occupancy for a tile to count as on (0.75)
%
% OUTPUT (struct O)
%   .occ {1 x Levels+1}  [nCh x n/2^L] occupancy per level      .tileS [1 x Levels+1] tile lengths (s)
%   .runs {1 x Levels+1} {nCh x 1} run durations (s)            .depth [nCh x n] (s; NaN where off)
%
% See also: rheome.detect.onstate, rheome.ingest.rollup
%
% Author: Diellor Basha, 2026

    arguments
        b logical
        dt (1,1) double {mustBePositive}
        opts.Levels = []
        opts.Theta (1,1) double {mustBeInRange(opts.Theta, 0, 1)} = 0.75
    end
    [nCh, n] = size(b);
    Lmax = floor(log2(n));  if ~isempty(opts.Levels), Lmax = min(Lmax, opts.Levels); end
    O.occ = cell(1, Lmax + 1);  O.runs = cell(1, Lmax + 1);  O.tileS = dt * 2.^(0:Lmax);
    O.depth = nan(nCh, n);
    x = double(b);
    for L = 0:Lmax
        if L > 0
            m = floor(size(x, 2) / 2);  x = (x(:, 1:2:2*m) + x(:, 2:2:2*m)) / 2;   % exact roll-up
        end
        O.occ{L+1} = x;  on = x >= opts.Theta;  if L == 0, on = x > 0.5; end
        r = cell(nCh, 1);
        for c = 1:nCh
            d = diff([false on(c, :) false]);  s0 = find(d == 1);  s1 = find(d == -1) - 1;
            r{c} = (s1 - s0 + 1)' * O.tileS(L+1);
            % depth: the finest tiles under an on-tile of this level get this level's length
            k = find(on(c, :));
            if ~isempty(k)
                idx = (k - 1) * 2^L + (1:2^L)';  idx = idx(idx <= n);
                keep = b(c, idx);  O.depth(c, idx(keep)) = O.tileS(L+1);
            end
        end
        O.runs{L+1} = r;
    end
end

% Author: Diellor Basha, 2026
