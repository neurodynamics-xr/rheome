function [keep, kind] = trackselect(stats, thr)
% DETECT.TRACKSELECT  Keep the trajectories that beat the noise-only null.
%
%   [keep, kind] = rheome.detect.trackselect(stats, thr)     % stats = rheome.detect.trackstats, thr = rheome.detect.tracknull
%
% A track is kept if it is PERSISTENT (durationS > thr.durationS) or MOVING (netMM > thr.netMM).
% `kind` says which, so a still object is reported as present-and-still and a brief fast one as moving,
% rather than both being ranked on one score.
%
% OUTPUT
%   keep  [nTracks x 1] logical
%   kind  [nTracks x 1] string: "moving", "persistent", "both" or ""
%
% See also: rheome.detect.tracknull, rheome.detect.trackstats, rheome.detect.tiletrack
%
% Author: Diellor Basha, 2026

    p = stats.durationS > thr.durationS;
    m = stats.netMM > thr.netMM;
    keep = p | m;
    kind = repmat("", height(stats), 1);
    kind(p & ~m) = "persistent";  kind(m & ~p) = "moving";  kind(p & m) = "both";
end

% Author: Diellor Basha, 2026
