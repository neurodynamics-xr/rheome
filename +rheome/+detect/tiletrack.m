function out = tiletrack(P, G, Y, opts)
% DETECT.TILETRACK  Link tile peaks across frames into trajectories: walks on the tile graph.
%
%   out = rheome.detect.tiletrack(P, G, Y)
%   out = rheome.detect.tiletrack(P, G, Y, LostFrames=30, Margin=0.1, MinFrames=10, SampleRate=300)
%
% P is rheome.detect.tilepeaks of Y. A track sitting on tile c may continue at frame t only to a peak on c
% itself or on a tile adjacent to c (G.A > 0). That is the whole motion model: at a rate and tile size
% chosen so a feature cannot cross a whole tile in one frame, adjacency lists every legal destination,
% and a jump to a non-adjacent tile is not motion -- it starts a new track.
%
% Assignment each frame is global (matchpairs), with cost 0 to stay, 1 to step to a neighbour, less a
% small reward for a higher peak, so two tracks cannot claim the same peak and a track prefers staying.
%
% ⭐ `Margin` IS A HYSTERESIS ON THE TILE GRAPH. A peak sitting on a boundary makes two tiles nearly
% equal and the argmax flickers between them. With Margin = h, a track on c steps to an adjacent peak
% d only when Y(d) >= (1 + h) * Y(c); otherwise it stays on c, even though c is momentarily the flank
% rather than the peak. It is expressed in the field's own units, so it does not depend on how fast the
% feature moves -- unlike a hold-for-n-frames debounce, which erases the hops of anything that crosses
% a tile in fewer than n frames.
%
% INPUTS
%   P           rheome.detect.tilepeaks output          G  rheome.geom.tiles struct (.A)
%   Y           [nTile x nT] the tile field P was detected on (needed for Margin)
%   LostFrames  frames a track survives with no peak before it ends (30 = one alpha cycle at 300 Hz)
%   Margin      tile-step hysteresis, fraction (0 = off)
%   MinFrames   drop tracks observed on fewer frames (3)
%   SampleRate  frames per second, for .duration (1)
%
% OUTPUT (struct out)
%   .tracks  struct array: frames tiles values held hops firstFrame lastFrame duration
%            `held` marks frames where Margin kept the track on a non-peak tile
%   .nTracks
%
% See also: rheome.detect.tilepeaks, rheome.geom.tilemean, rheome.geom.tiles, rheome.detect.track
%
% Author: Diellor Basha, 2026

    arguments
        P struct
        G struct
        Y double
        opts.LostFrames (1,1) double {mustBeNonnegative} = 30
        opts.Margin (1,1) double {mustBeNonnegative} = 0
        opts.MinFrames (1,1) double {mustBePositive} = 3
        opts.SampleRate (1,1) double {mustBePositive} = 1
    end
    nT = numel(P.tiles);  A = G.A > 0;
    tr = struct('frames', {}, 'tiles', {}, 'values', {}, 'held', {}, 'last', {}, 'open', {});
    for t = 1:nT
        pk = P.tiles{t}(:);  pv = P.values{t}(:);
        open = find([tr.open]);
        claimed = false(numel(pk), 1);  moved = false(1, numel(open));
        if ~isempty(open) && ~isempty(pk)
            C = inf(numel(open), numel(pk));
            for a = 1:numel(open)
                c = tr(open(a)).tiles(end);
                for b = 1:numel(pk)
                    if pk(b) == c
                        C(a, b) = 0;
                    elseif A(c, pk(b)) && Y(pk(b), t) >= (1 + opts.Margin) * Y(c, t)
                        C(a, b) = 1 - 0.1 * pv(b) / max(abs(pv(1)), realmin);
                    end
                end
            end
            C(~isfinite(C)) = 1e6;                     % matchpairs needs finite costs
            M = matchpairs(C, 2);
            M = M(C(sub2ind(size(C), M(:,1), M(:,2))) < 1e6, :);
            for r = 1:size(M, 1)
                k = open(M(r,1));  b = M(r,2);
                tr(k) = i_add(tr(k), t, pk(b), pv(b), false);
                claimed(b) = true;  moved(M(r,1)) = true;
            end
        end
        % hysteresis: an unmatched track whose own tile is still (1+Margin)-dominant over every
        % adjacent peak stays put rather than being dropped
        for a = find(~moved)
            k = open(a);  c = tr(k).tiles(end);
            if opts.Margin > 0 && t - tr(k).last == 1
                nb = pk(A(c, pk));
                if ~isempty(nb) && all(Y(nb, t) < (1 + opts.Margin) * Y(c, t))
                    tr(k) = i_add(tr(k), t, c, Y(c, t), true);
                    claimed(ismember(pk, nb)) = true;          % the flank is this track, not a birth
                    continue
                end
            end
            if t - tr(k).last > opts.LostFrames, tr(k).open = false; end
        end
        for b = find(~claimed)'
            tr(end+1) = struct('frames', t, 'tiles', pk(b), 'values', pv(b), 'held', false, ...
                'last', t, 'open', true); %#ok<AGROW>
        end
    end
    keep = arrayfun(@(x) numel(x.frames) >= opts.MinFrames, tr);
    tr = tr(keep);
    out.tracks = struct('frames', {}, 'tiles', {}, 'values', {}, 'held', {}, 'hops', {}, ...
        'firstFrame', {}, 'lastFrame', {}, 'duration', {});
    for k = 1:numel(tr)
        x = tr(k);
        out.tracks(k) = struct('frames', x.frames(:), 'tiles', x.tiles(:), 'values', x.values(:), ...
            'held', x.held(:), 'hops', sum(diff(x.tiles) ~= 0), 'firstFrame', x.frames(1), ...
            'lastFrame', x.frames(end), 'duration', (x.frames(end) - x.frames(1)) / opts.SampleRate);
    end
    out.nTracks = numel(out.tracks);
end

function x = i_add(x, t, tile, val, held)
    x.frames(end+1) = t;  x.tiles(end+1) = tile;  x.values(end+1) = val;  x.held(end+1) = held;
    x.last = t;
end

% Author: Diellor Basha, 2026
