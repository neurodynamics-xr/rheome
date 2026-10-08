function C = cleanspan(F, fs, opts)
% SCALE.CLEANSPAN  Where a recording is stationary enough to tile: artefact blocks and trimmed edges.
%
%   C = rheome.scale.cleanspan(F, fs)
%   C = rheome.scale.cleanspan(F, fs, BlockS=0.5, Factor=3, EdgeS=10, Dilate=1, LowHz=4)
%
% Cuts the recording into BlockS blocks, takes the MEDIAN over channels of each block's RMS about
% the channel's recording median, and calls a block bad when that exceeds Factor x its median over blocks --
% measured twice, broadband and below LowHz, and a block is bad if either level says so.
% Bad blocks are dilated by Dilate blocks on each side. A bad block inside the first or last
% EdgeS seconds moves the usable start (end) past it, so tiles are never placed on a recording
% edge; a bad block further in is left in .mask for the caller to reject the tile that touches it.
%
% ⚠ WHY THIS EXISTS. About 3% of the subjects of a resting cohort failed the periodicflow measure
% with spectral:decompose:scale on their LAST tile, whose window ends on the last
% sample. On one of them the final ~3.5 s carry a common-mode artefact (first PC 96% of the variance,
% against 18% mid-recording) at 150x-340000x the mid-recording power from 300 Hz down to 1 Hz,
% and the first 0.5 s sits at 10x. It is data, not units: on clean tiles the data-derived |C|^2/PSD
% scale equals the analytic fs*nT/2 to 0.88-1.13. The same edges bias the subjects that passed: tile 12
% read > 60% background-above-observed in ~5% of them (<1% for any other tile) and tile 1's median
% was 13.6% against 27% for every other tile.
%
% ⚠ THE BROADBAND LEVEL ALONE SEES AN ARTEFACT 2.5 s LATE. Broadband std is dominated by the high
% frequencies; that subject's end artefact is red, and its sub-1 Hz level rises at 120.5 s while the
% broadband level first crosses 3x at 123 s. Below the fit band is not harmless: periodicflow's
% unwindowed tile FFT leaks it into 1-45 Hz, and a tile ending at 122.5 s read 3.5x the recording's
% power. Hence the second, low-passed level. On the 600 Hz cache the start shows the same thing:
% 1-4 Hz at 204x / 41x / 8x in the first three blocks, broadband 7x / 1.3x / 1.0x.
%
% INPUTS:
%   F   [nCh x N] recording, channels in rows
%   fs  sampling rate (Hz)
% OUTPUT (struct C):
%   .first .last   usable sample range [first, last] after edge trimming (1-based, inclusive)
%   .mask          [1 x N] logical, true on samples of a (dilated) bad block
%   .ratio         [1 x nB] block level / median block level, the larger of broadband and < LowHz
%   .bad           [1 x nB] logical, the dilated bad blocks
%   .blockS .typical .nTrimStart .nTrimEnd   (trims in samples)
%
% See also: rheome.scale.measure_periodicflow
%
% Author: Diellor Basha, 2026

    arguments
        F double
        fs (1,1) double {mustBePositive}
        opts.BlockS (1,1) double {mustBePositive} = 0.5
        opts.Factor (1,1) double {mustBePositive} = 3
        opts.EdgeS (1,1) double {mustBeNonnegative} = 10
        opts.Dilate (1,1) double {mustBeInteger, mustBeNonnegative} = 1
        opts.LowHz (1,1) double {mustBeNonnegative} = 4      % 0 = broadband level only
    end
    N = size(F, 2);  w = max(1, round(opts.BlockS * fs));  nB = ceil(N / w);
    lev = i_level(F, w, nB, N);  typ = median(lev);  ratio = lev / max(typ, realmin);
    if opts.LowHz > 0 && opts.LowHz < fs/2 && N > 30
        [b, a] = butter(4, opts.LowHz / (fs/2), 'low');
        ratio = max(ratio, i_ratio(filtfilt(b, a, F.').', w, nB, N));
    end
    bad0 = ratio > opts.Factor;
    bad = bad0;
    for d = 1:opts.Dilate
        bad = bad | [bad0(1+d:end) false(1,d)] | [false(1,d) bad0(1:end-d)];
    end
    mask = repelem(bad, w);  mask = mask(1:N);

    nE = min(nB, ceil(opts.EdgeS * fs / w));
    first = 1;  last = N;
    kS = find(bad(1:nE), 1, 'last');
    if ~isempty(kS), first = min(N, kS*w + 1); end
    kE = find(bad(nB-nE+1:nB), 1, 'first');
    if ~isempty(kE), last = max(1, (nB - nE + kE - 1)*w); end

    C = struct('first', first, 'last', last, 'mask', mask, 'ratio', ratio, 'bad', bad, ...
               'blockS', w / fs, 'typical', typ, 'nTrimStart', first - 1, 'nTrimEnd', N - last);
end

function lev = i_level(F, w, nB, N)
% median over channels of each block's RMS about the channel's RECORDING median
% ⚠ not the block std: a swell slower than the block is a shift of the block mean, which std removes
    F = F - median(F, 2);
    lev = zeros(1, nB);
    for k = 1:nB
        lev(k) = median(sqrt(mean(F(:, (k-1)*w + 1 : min(k*w, N)).^2, 2)));
    end
end

function r = i_ratio(F, w, nB, N)
    lev = i_level(F, w, nB, N);  r = lev / max(median(lev), realmin);
end

% Author: Diellor Basha, 2026
