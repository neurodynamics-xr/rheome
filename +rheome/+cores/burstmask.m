function [keep, info] = burstmask(z, good, per, opts)
% CORES.BURSTMASK  The frames inside bursts by a percentile threshold (nxr's older burst_mask).
%
%   [keep, info] = rheome.cores.burstmask(z, good, per)
%
% The band amplitude sqrt(mean over channels |z|^2) smoothed over one cycle (rheome.cores.boxsmooth, zero-padded) at or above its
% Pct-th percentile (75) over the good frames (numpy's linear rule, rheome.cores.linpercentile; inclusive), AND good, in runs of
% at least MinCycles (3) cycles. keep [1 x m] logical. info: .nBursts, .burstCycles (median run length in cycles, NaN when none),
% .burstShare (kept frames / good frames).
%
% Ports nxr-analysis src/nxr/vortex.py burst_mask (line 784, commit 30b9622).
%
% See also: rheome.cores.burststages, rheome.cores.linpercentile
%
% Author: Diellor Basha, 2026
    arguments
        z (:,:) {mustBeNumeric}
        good (1,:) logical
        per (1,1) double {mustBeInteger, mustBePositive}
        opts.Pct (1,1) double = 75
        opts.MinCycles (1,1) double = 3
    end
    a = sqrt(mean(abs(double(z)).^2, 1));
    a = rheome.cores.boxsmooth(a, per);
    on = (a >= rheome.cores.linpercentile(a(good), opts.Pct)) & good;
    keep = false(size(on));  lens = [];
    r = runsof(on);
    for q = 1:size(r, 1)
        if r(q, 2) - r(q, 1) >= opts.MinCycles * per
            keep(r(q,1)+1:r(q,2)) = true;  lens(end+1) = (r(q,2) - r(q,1)) / per; %#ok<AGROW>
        end
    end
    if isempty(lens), mc = NaN; else, mc = median(lens); end
    info = struct('nBursts', numel(lens), 'burstCycles', mc, 'burstShare', nnz(keep) / max(nnz(good), 1));
end
