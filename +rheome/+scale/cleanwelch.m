function [P, f, info] = cleanwelch(F, fs, span, mask, nf)
% SCALE.CLEANWELCH  Welch PSD of a span, averaging only the segments that touch no artefact block.
%
%   [P, f] = rheome.scale.cleanwelch(F, fs, first:last, mask, nf)
%   [P, f, info] = rheome.scale.cleanwelch(...)      % info.nSeg, info.nUsed, info.used
%
% pwelch(F(:, span).', hann(nf), nf/2, nf, fs) with the segments that overlap a true sample of
% mask left out of the average. The segments are pwelch's own: start 1 + (k-1)*nf/2 within the
% span, k = 1 .. fix((L - nf/2) / (nf/2)).
%
% ⭐ WHEN NO SEGMENT TOUCHES THE MASK THIS IS THE pwelch CALL ITSELF, bit for bit -- the masked
% path is only taken when it has something to exclude.
%
% ⚠ WHY THIS EXISTS. periodicflow screens each tile by power_ratio = (tile's fit-band |C|^2 /
% this PSD) / (fs*nT/2), and fits the aperiodic background on this PSD. rheome.scale.cleanspan masks
% interior artefact blocks so no TILE touches them, but before 2026-10-07 the Welch average still
% included them. On one subject, one 8 s segment (an artefact at 20-25 s, 50x the median block) read
% 806x the median segment's fit-band power, the mean PSD 20x a typical stretch, and EVERY clean
% tile then read power_ratio 0.03-0.07 and was rejected. On another (artefacts at 74, 160 and
% 245-274 s): mean 4.4x, every clean tile 0.15-0.33. The artefact is red (first subject: tile/PSD
% 0.01 at 1-2 Hz, 0.2 at 8-45 Hz), so it bends the aperiodic fit as well as its level.
%
% INPUTS:
%   F     [nCh x N] recording, channels in rows
%   fs    sampling rate (Hz)
%   span  [1 x L] contiguous sample indices (rheome.scale.cleanspan's first:last)
%   mask  [1 x N] logical, true on artefact samples (rheome.scale.cleanspan's .mask)
%   nf    segment length = nfft (samples); overlap nf/2, Hann window
% OUTPUTS:
%   P     [nFreq x nCh] one-sided PSD,  f [nFreq x 1] Hz -- pwelch's shapes
%   info  .nSeg (pwelch's segment count) .nUsed (clean segments averaged) .used [1 x nSeg] logical
%
% See also: rheome.scale.cleanspan, rheome.scale.measure_periodicflow, pwelch
%
% Author: Diellor Basha, 2026

    L = numel(span);  hop = nf / 2;
    if L < nf, error('scale:cleanwelch:short', 'span of %d samples is shorter than one %d-sample segment', L, nf); end
    nSeg = fix((L - hop) / hop);  st = 1 + (0:nSeg-1) * hop;
    used = true(1, nSeg);
    m = mask(span);
    for k = 1:nSeg, used(k) = ~any(m(st(k) : st(k) + nf - 1)); end
    info = struct('nSeg', nSeg, 'nUsed', sum(used), 'used', used);
    if all(used)
        [P, f] = pwelch(F(:, span).', hann(nf), hop, nf, fs);
        return
    end
    if ~any(used)
        error('scale:cleanwelch:nosegments', 'every one of the %d Welch segments touches an artefact block', nSeg);
    end
    P = 0;  w = hann(nf);
    for k = find(used)
        [Pk, f] = pwelch(F(:, span(st(k) : st(k) + nf - 1)).', w, 0, nf, fs);
        P = P + Pk;
    end
    P = P / info.nUsed;
end

% Author: Diellor Basha, 2026
