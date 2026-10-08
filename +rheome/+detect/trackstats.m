function S = trackstats(out, G, rate)
% DETECT.TRACKSTATS  Per-track features of a tile trajectory: how long, how far, how straight.
%
%   S = rheome.detect.trackstats(out, G, rate)      % out = rheome.detect.tiletrack(...), G = rheome.geom.tiles(..., Ruler=)
%
% The features a trajectory is judged by. A still object and a moving one are both PERSISTENT; only
% the mover has DISPLACEMENT; noise produces tracks that are short, or long but going nowhere. Every
% distance is between tile centres on G's ruler (G.D), so it is read on the same ruler as the truth.
%
% OUTPUT (table S, one row per track)
%   track       index into out.tracks
%   nFrames     frames observed        durationS  (last - first) / rate
%   netMM       distance from the first tile's centre to the last's
%   pathMM      sum of the distances of each tile change
%   straight    netMM / pathMM (1 for a straight walk; NaN with no hops)
%   hops        tile changes           meanValue  mean tile value along the track
%   speedMS     netMM / durationS -- a NET speed: jitter cancels, so it is the one to quote
%
% ⚠ `netMM` is quantised to the tiling: a track that never leaves its tile has netMM = 0 however its
% object moves inside it. Displacement below one tile diameter is not measured here, by design.
%
% See also: rheome.detect.tiletrack, rheome.detect.tracknull, rheome.detect.trackselect, rheome.geom.tiles
%
% Author: Diellor Basha, 2026

    if isempty(G.D)
        error('detect:trackstats:ruler', 'G has no .D: build it with rheome.geom.tiles(..., Ruler=rheome.geom.edgegraph(S)).');
    end
    n = numel(out.tracks);
    v = zeros(n, 9);
    for q = 1:n
        x = out.tracks(q);  t = x.tiles(:);
        dur = (x.frames(end) - x.frames(1)) / rate;
        net = G.D(t(1), t(end));
        ch = find(diff(t) ~= 0);
        pth = sum(G.D(sub2ind(size(G.D), t(ch), t(ch + 1))));
        st = NaN;  if pth > 0, st = net / pth; end
        v(q,:) = [q numel(x.frames) dur 1e3*net 1e3*pth st numel(ch) mean(x.values) net/max(dur, eps)];
    end
    S = array2table(v, 'VariableNames', {'track','nFrames','durationS','netMM','pathMM','straight', ...
        'hops','meanValue','speedMS'});
end

% Author: Diellor Basha, 2026
