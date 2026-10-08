function T = support(W, opts)
% INGEST.SUPPORT  Which bands a window of a given length can resolve.
%
%   T = rheome.ingest.support(2)                       % a 2 s window
%   T = rheome.ingest.support([2 16 64], Voices=4)     % several, at the bank's voices
%   [~] = rheome.ingest.support(2, Anchor=1)
%
% ⭐ A WINDOW LENGTH IS A STATEMENT ABOUT FREQUENCY. Constant Q ties a band's time support to
% its centre frequency: support = k / f, with k fixed by the bank's voices (measured on
% @timefilterbank: support*fc = 4.8, 8.0, 11.5, 15.0, 22.4, 37.1 at 1, 2, 3, 4, 6, 10
% voices). So a window of W seconds resolves everything above k / W and nothing below it.
% Choosing a stretch in a preview therefore already chooses the analysis that fits in it,
% which is why the explorer prints this next to the window.
%
% Returns one row per window: the lowest frequency it supports, the lowest OCTAVE band fully
% inside it (anchored the way the store's bands are), and how many cycles of that band the
% window holds.
%
% ⚠ THIS IS THE BANK'S SUPPORT, NOT A DISPLAY LIMIT. A shorter window still SHOWS low
% frequencies; it cannot separate them, because the kernel that would do it is longer than
% the window. The store's own tile for a band (rheome.select.ladder) is the first dyadic level at
% least this long.
%
% See also: rheome.select.ladder, rheome.ingest.preview, rheome.timefilterbank
%
% Author: Diellor Basha, 2026

    arguments
        W double {mustBePositive}
        opts.Voices (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.Anchor (1,1) double {mustBePositive} = 1
    end
    k = i_timesfc(opts.Voices);
    W = double(W(:));
    fLow = k ./ W;                                                   % the lowest f with support <= W
    oct = floor(log2(fLow / opts.Anchor)) + 1;                       % the lowest octave fully inside
    bLo = opts.Anchor * 2.^oct;
    T = table(W, repmat(k, numel(W), 1), fLow, bLo, 2*bLo, W .* sqrt(bLo .* 2 .* bLo), ...
              'VariableNames', {'window_s','support_times_fc','f_lowest','band_lo','band_hi','cycles_in_window'});
end

% The measured time-bandwidth product of the frame bank's kernels, by voices (99.9 % of the
% energy). Interpolated in log for a voice count not measured; the four-voice default is the
% store's.
function k = i_timesfc(V)
    v = [1 2 3 4 6 10];  s = [4.8 8.0 11.5 15.0 22.4 37.1];
    if any(v == V), k = s(v == V); return; end
    k = exp(interp1(log(v), log(s), log(V), 'linear', 'extrap'));
end
% Author: Diellor Basha, 2026
