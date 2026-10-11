function [z, rate, valid] = analytic(Y, fs, edges, mask)
% MODES.ANALYTIC  Band-pass, Hilbert, decimate: the complex band signal of tile series.
%
%   [z, rate, valid] = rheome.modes.analytic(Y, fs, [8 13])
%   [z, rate, valid] = rheome.modes.analytic(Y, fs, [8 13], mask)
%
% Order-8 Butterworth band-pass applied forward and backward (zero phase, as the MS1 measures), analytic
% signal by Hilbert, decimated to the largest rate fs/m >= 3 x the upper edge (the analytic signal is
% band-limited, so no anti-alias filter is needed). valid is false on mask (decimated, any sample of the
% block) and within 3 cycles of the lower edge of either end (filter transients).
%
%   z [nTiles x T] single complex   rate (Hz)   valid [1 x T] logical
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(mask), mask = false(1, size(Y, 2)); end
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', edges(1), ...
                    'HalfPowerFrequency2', edges(2), 'SampleRate', fs);
    Za = hilbert(filtfilt(bp, Y.')).';
    m = max(1, floor(fs / (3 * edges(2))));  rate = fs / m;
    z = single(Za(:, 1:m:end));  T = size(z, 2);
    bad = any(reshape([logical(mask(:)'), false(1, m*T - numel(mask))], m, T), 1);
    e = ceil(3 / edges(1) * rate);  bad([1:min(e, T), max(1, T-e+1):T]) = true;
    valid = ~bad;
end

% Author: Diellor Basha, 2026
