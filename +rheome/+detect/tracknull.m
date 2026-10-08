function thr = tracknull(nullStats, nullSeconds, opts)
% DETECT.TRACKNULL  Selection thresholds from tracks found in NOISE ALONE.
%
%   thr = rheome.detect.tracknull(nullStats, nullSeconds)              % Alpha = 0.1 per second
%   thr = rheome.detect.tracknull(nullStats, nullSeconds, Alpha=0.05)
%
% nullStats is rheome.detect.trackstats pooled over noise-only records totalling nullSeconds. A track is
% worth reporting if it is PERSISTENT (lasts longer than noise tracks do) or MOVING (gets further than
% noise tracks do). Each rule gets half of Alpha, a false-trajectory RATE in tracks per second of
% record: the threshold is the smallest value that fewer than Alpha/2 * nullSeconds noise tracks
% exceed.
%
% ⭐ WHY TWO RULES AND A UNION. A still object is persistent and goes nowhere; a fast one over a short
% path is brief and goes far -- plant_tiletrack_omega.m's 1 m/s mover lived 0.12 s, no longer than a
% noise bump, so a duration-only rule (the "longest track") picked noise. Neither rule alone keeps
% both; the union does, at a false rate of at most Alpha.
%
% ⭐ A RATE, NOT A PER-TRACK p-VALUE. Noise produces tens of tracks per second at 10 dB, so a per-track
% 5% would pass several false trajectories every second. What a reader of a trajectory map needs is
% how many of the reported trajectories per second are noise, and that is what Alpha sets.
%
% ⚠ If the null holds too few tracks to reach Alpha/2 (a short null, or a clean one), the threshold is
% the null's maximum and `.saturated` says so: the false rate is then bounded only by 1/nullSeconds.
%
% OUTPUT (struct thr)
%   .durationS  .netMM       the two thresholds (a track passes on value > threshold)
%   .alpha  .nullSeconds  .nNull   .saturated  [1 x 2] (duration, net)
%
% See also: rheome.detect.trackselect, rheome.detect.trackstats
%
% Author: Diellor Basha, 2026

    arguments
        nullStats table
        nullSeconds (1,1) double {mustBePositive}
        opts.Alpha (1,1) double {mustBePositive} = 0.1
    end
    allow = floor(opts.Alpha / 2 * nullSeconds);        % null tracks allowed above each threshold
    [thr.durationS, s1] = i_thr(nullStats.durationS, allow);
    [thr.netMM, s2] = i_thr(nullStats.netMM, allow);
    thr.alpha = opts.Alpha;  thr.nullSeconds = nullSeconds;  thr.nNull = height(nullStats);
    thr.saturated = [s1 s2];
end

% The smallest threshold with at most `allow` values strictly above it.
function [t, sat] = i_thr(x, allow)
    x = sort(x(:), 'descend');
    if isempty(x), t = 0; sat = true; return; end
    if allow >= numel(x), t = 0; sat = false; return; end
    t = x(allow + 1);  sat = allow == 0;
end

% Author: Diellor Basha, 2026
