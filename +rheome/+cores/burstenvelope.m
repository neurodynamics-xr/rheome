function [c, th] = burstenvelope(a, good, per, opts)
% CORES.BURSTENVELOPE  Candidate bursts of an amplitude envelope, durations in cycles (nxr spindle_envelope).
%
%   c = rheome.cores.burstenvelope(a, good, per)
%   [c, th] = rheome.cores.burstenvelope(a, good, per, HiSD=1.5, LoSD=1, MinCycles=3, MergeCycles=1, SmoothCycles=2)
%
% a [1 x m] the sensors' band amplitude (after a 1-cycle smoothing in burststages); good [1 x m] logical; per frames per cycle.
% The envelope is smoothed again (rheome.cores.boxsmooth, SmoothCycles*per frames, zero-padded); a robust global baseline over
% the good frames: median and 1.4826 x MAD (MAD about that median) -- never a percentile, and NOT guarded against a zero MAD
% (rheome.detect.spindle skips the channel then). hi = median + HiSD*SD, lo = median + LoSD*SD. Runs of (e >= lo) & good are kept
% when any frame is >= hi; a run closer than MergeCycles*per frames to the previous (gap = start - previous stop, in frames, less
% than the limit) is merged into it, before the length gate; kept when its length is >= MinCycles*per frames. No maximum length.
%
% OUTPUT c [n x 2] first and last frame of each burst, 1-BASED INCLUSIVE (nxr returns 0-based [start, stop): start+1 and stop).
% Second output th = [median, robust SD, hi, lo] of the smoothed envelope over the good frames.
%
% Differs from rheome.detect.spindle (which is unchanged): thresholds 1.5 / 1.0 SD (spindle: 3.0 / 2.5), durations in cycles (not
% seconds), envelope branch only, one channel, baseline always the good frames, no maximum duration. Ports nxr-analysis
% src/nxr/vortex.py spindle_envelope (line 829, commit 30b9622).
%
% See also: rheome.cores.burststages, rheome.cores.burstfit, rheome.detect.spindle
%
% Author: Diellor Basha, 2026
    arguments
        a (1,:) double
        good (1,:) logical
        per (1,1) double {mustBeInteger, mustBePositive}
        opts.HiSD (1,1) double = 1.5
        opts.LoSD (1,1) double = 1.0
        opts.MinCycles (1,1) double = 3
        opts.MergeCycles (1,1) double = 1
        opts.SmoothCycles (1,1) double = 2
    end
    k = max(1, evenround(opts.SmoothCycles * per));
    e = rheome.cores.boxsmooth(a, k);
    med = median(e(good));
    sd = 1.4826 * median(abs(e(good) - med));
    hi = med + opts.HiSD * sd;  lo = med + opts.LoSD * sd;
    th = [med sd hi lo];
    r = runsof((e >= lo) & good);
    merged = zeros(0, 2);
    for q = 1:size(r, 1)
        i = r(q, 1);  j = r(q, 2);
        if ~any(e(i+1:j) >= hi), continue; end
        if ~isempty(merged) && i - merged(end, 2) < opts.MergeCycles * per
            merged(end, 2) = j;
        else
            merged(end+1, :) = [i j]; %#ok<AGROW>
        end
    end
    keep = (merged(:, 2) - merged(:, 1)) >= opts.MinCycles * per;
    merged = merged(keep, :);
    c = [merged(:, 1) + 1, merged(:, 2)];
end
