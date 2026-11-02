function [thr, hist] = threshold(C, opts)
% CORES.THRESHOLD  The |curl| threshold of the core detector: a percentile read from a log histogram (nxr _core_threshold).
%
%   thr = rheome.cores.threshold(C)                  % C [V x T] single curl movie, 98th percentile
%   [thr, hist] = rheome.cores.threshold(C, Pct=98, Chunk=3000)
%
% 4096 bins: bin 0 for |c| < lo, then geometric up to hi = 4*max|c| of the FIRST chunk of Chunk frames (lo = 1e-6*hi); the
% index is floor(log(|c|/lo)*4095/log(hi/lo))+1, capped at 4095; values clipped to hi. The threshold is the point of
% geomspace(lo,hi,4096) at the first bin whose cumulative fraction reaches Pct/100 (a grid point, not a bin edge).
% The binning is done in single, as nxr does. Cores are then C >= thr and C <= -thr (inclusive).
%
% Ports nxr-analysis src/nxr/vortex.py _core_threshold (commit 903fdc2).
%
% See also: rheome.cores.patches, rheome.cores.stats
%
% Author: Diellor Basha, 2026
    arguments
        C
        opts.Pct (1,1) double = 98
        opts.Chunk (1,1) double = 3000
    end
    C = single(C); nb = 4096; nT = size(C, 2); hist = zeros(nb, 1);
    hi = []; 
    for i = 1:opts.Chunk:nT
        a = abs(C(:, i:min(i + opts.Chunk - 1, nT)));  a = a(:);
        if isempty(hi)
            hi = double(max(a)) * 4;  lo = hi * 1e-6;  sc = (nb - 1) / log(hi / lo);
        end
        a = min(a, single(hi));
        b = zeros(numel(a), 1);
        m = a >= single(lo);
        b(m) = min(floor(log(a(m) / single(lo)) * single(sc)) + 1, nb - 1);
        hist = hist + accumarray(b + 1, 1, [nb 1]);
    end
    grid = geomspace(lo, hi, nb);
    k = find(cumsum(hist) / sum(hist) >= opts.Pct / 100, 1);
    thr = grid(k);
end

function x = geomspace(lo, hi, n)
% numpy geomspace: 10.^linspace(log10(lo), log10(hi)), end points exact
    x = 10 .^ linspace(log10(lo), log10(hi), n);  x(1) = lo;  x(end) = hi;
end
