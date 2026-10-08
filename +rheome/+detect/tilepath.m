function out = tilepath(Y, G, opts)
% DETECT.TILEPATH  Track-before-detect on the tile graph: the best-scoring walk through a window.
%
%   out = rheome.detect.tilepath(Y, G, Baseline=b)
%   out = rheome.detect.tilepath(Y, G, Baseline=b, StepCost=0, NumPaths=3, SampleRate=300)
%
% Y is a tile field over ONE observation window ([nTile x nF], rheome.geom.tilemean then a dyadic block);
% G is rheome.geom.tiles. A path is a walk on the tile graph -- each frame it stays or steps to an adjacent
% tile -- and its score is the sum of (Y - Baseline) along it, less StepCost per step. The best path
% is found exactly by dynamic programming (Viterbi), with a local-alignment reset so it may START AND
% END ANYWHERE in the window:
%
%     S(k,t) = Y(k,t) - b + max(0, max_{j in {k} u adj(k)} S(j,t-1) - StepCost*[j ~= k])
%
% ⭐⭐ WHY NOT DETECT, THEN LINK. plant_select_omega.m: at 10 dB the peak falls under a noise blob every
% few frames, so a chain of per-frame peaks (rheome.detect.tiletrack) breaks, and each fragment of a 120 mm
% mover covered 36-88 mm -- no further than noise tracks of the same length. The information was lost
% when each frame committed to its own peaks. Here no frame commits: a frame where noise wins costs
% the path a little score instead of breaking it, and evidence accumulates along the motion. It is
% the tile-graph form of velocity-matched integration, i.e. of the joint speed filter -- the window
% (tile size, rate) comes from rheome.geom.jointcell, the path is the measurement made inside it.
%
% ⭐ THE BASELINE MAKES NOISE LOSE. Without it every value adds, and the best path is simply the
% longest. With b above typical noise, a walk through noise scores negative on average -- it can
% follow a noise blob for its lifetime, then pays for the gap -- while a walk along a real peak gains.
% Set b from noise alone (a high quantile of noise-only tile values), not by hand.
%
% ⭐ The adjacency-only step is the motion model and it is exact for the window: the ladder's rate is
% chosen so a trackable feature moves at most one tile per frame (speedTrackMS), so no real path needs
% a longer step.
%
% INPUTS
%   Y           [nTile x nF] tile values        G  rheome.geom.tiles (.A)
%   Baseline    subtracted from every value (0)
%   StepCost    charged per tile change (0)
%   NumPaths    extract this many disjoint paths, best first (1); later paths may not use a
%               (tile, frame) an earlier one used
%   SampleRate  frames per second, for .duration (1)
%   Span        "local" (default): the path may start and end anywhere; "full": the path must run
%               from the first frame to the last -- for a window that is ALREADY the event, e.g. an
%               alpha episode from its fitted onset to its offset, where dropping frames would drop
%               the phases being compared
%
% OUTPUT (struct out)  -- the shape of rheome.detect.tiletrack's output, so rheome.detect.trackstats applies
%   .tracks  struct array: frames tiles values held hops firstFrame lastFrame duration score
%   .nTracks
%
% See also: rheome.detect.tiletrack, rheome.detect.trackstats, rheome.detect.tracknull, rheome.geom.jointcell
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        G struct
        opts.Baseline (1,1) double = 0
        opts.StepCost (1,1) double {mustBeNonnegative} = 0
        opts.NumPaths (1,1) double {mustBeInteger, mustBePositive} = 1
        opts.SampleRate (1,1) double {mustBePositive} = 1
        opts.Span (1,1) string {mustBeMember(opts.Span, ["local","full"])} = "local"
    end
    spanFull = opts.Span == "full";
    [nTile, nF] = size(Y);
    A = G.A > 0;  A(1:nTile+1:end) = false;
    deg = full(sum(A, 1));  NB = ones(nTile, max(deg));  pad = true(nTile, max(deg));
    for k = 1:nTile, nb = find(A(:, k));  NB(k, 1:numel(nb)) = nb;  pad(k, 1:numel(nb)) = false; end
    Z = Y - opts.Baseline;
    used = false(nTile, nF);
    out.tracks = struct('frames', {}, 'tiles', {}, 'values', {}, 'held', {}, 'hops', {}, ...
        'firstFrame', {}, 'lastFrame', {}, 'duration', {}, 'score', {});
    for n = 1:opts.NumPaths
        Zn = Z;  Zn(used) = -Inf;
        S = -Inf(nTile, nF);  back = zeros(nTile, nF);          % 0 = the path starts here
        S(:, 1) = Zn(:, 1);
        for t = 2:nF
            prev = S(:, t-1);
            pn = prev(NB) - opts.StepCost;  pn(pad) = -Inf;       % [nTile x maxDeg] neighbours
            [bn, ib] = max(pn, [], 2);
            step = bn > prev;
            best = prev;  best(step) = bn(step);
            from = (1:nTile)';  from(step) = NB(sub2ind(size(NB), find(step), ib(step)));
            go = best > 0 | spanFull;                               % else the path starts here
            S(:, t) = Zn(:, t);  S(go, t) = S(go, t) + best(go);   % (0 * -Inf would be NaN)
            back(go, t) = from(go);
        end
        if spanFull
            [sc, k] = max(S(:, nF));  t = nF;
            if ~isfinite(sc), break; end
        else
            [sc, idx] = max(S(:));
            if ~isfinite(sc) || sc <= 0, break; end
            [k, t] = ind2sub(size(S), idx);
        end
        tiles = k;  frames = t;
        while t > 1 && back(k, t) ~= 0
            k = back(k, t);  t = t - 1;
            tiles(end+1) = k;  frames(end+1) = t; %#ok<AGROW>
        end
        tiles = flip(tiles(:));  frames = flip(frames(:));
        used(sub2ind(size(used), tiles, frames)) = true;
        v = Y(sub2ind(size(Y), tiles, frames));
        out.tracks(end+1) = struct('frames', frames, 'tiles', tiles, 'values', v, ...
            'held', false(size(frames)), 'hops', sum(diff(tiles) ~= 0), 'firstFrame', frames(1), ...
            'lastFrame', frames(end), 'duration', (frames(end) - frames(1)) / opts.SampleRate, 'score', sc);
    end
    out.nTracks = numel(out.tracks);
end

% Author: Diellor Basha, 2026
