function [Y, fsOut, info] = bandlimit(X, fs, band, opts)
% FILTERS.BANDLIMIT  Temporal band-pass + rate matching, with a DECLARED transfer function.
%
%   [Y, fsOut, info] = rheome.filters.bandlimit(X, fs, [f1 f2])
%   [Y, fsOut, info] = rheome.filters.bandlimit(X, fs, [f1 f2], opts)
%
% WHY THIS EXISTS. Feature tracking on a band-limited field must be sampled at a rate matched to
% the BAND, not to the acquisition. Sampling faster adds interpolation, not information: the
% frame-to-frame displacement of a detected feature is then position-estimation error d rather
% than motion, and the apparent speed it produces is
%
%       v_floor ~ d * f_track                         (grows LINEARLY with tracking rate)
%
% so tracking a 5 Hz-wide band at 2400 Hz manufactures metres per second out of sub-millimetre
% jitter. This function band-limits and then decimates to a rate set by the band, and RETURNS the
% gain it applied, so the scale of the measurement is reported rather than implied.
%
% The filter is a zero-phase raised-cosine band-pass applied in the frequency domain: an explicit
% multiplication by a gain vector g(f), invertible wherever g > 0, with no phase distortion and no
% filter-order/edge heuristics. Decimation after band-limiting is plain subsampling and is
% alias-free because everything above f2+w has already been removed.
%
% INPUTS:
%   X      [nRow x nTime] signal (rows = channels or source components)
%   fs     input sampling rate (Hz)
%   band   [f1 f2] pass band in Hz (e.g. [8 13] for alpha)
%   opts   .oversample  target rate = oversample * f2   (default 8; >=2 by Nyquist, 8 gives
%                       ~8 samples per cycle at the top of the band)
%          .transition  raised-cosine roll-off width as a fraction of the band width
%                       (default 0.25; for [8 13] that is 1.25 Hz on each side)
%          .decimate    explicit integer decimation factor, overriding .oversample
%          .detrend     remove the per-row mean first (default true)
%
% OUTPUTS:
%   Y      [nRow x nOut] band-limited, decimated signal
%   fsOut  output sampling rate (Hz) = fs / info.decimation
%   info   .band .oversample .transition .decimation .fsIn .fsOut
%          .gain  [1 x nTime] the applied one-sided gain (the declared transfer function)
%          .freq  [1 x nTime] the frequency axis the gain is defined on
%          .passbandFraction  fraction of input variance retained (a lossiness receipt)
%
% EXAMPLE (alpha, from 2400 Hz):
%   [Fa, fsA, info] = rheome.filters.bandlimit(F(good,:), 2400, [8 13]);
%   % -> fsA ~ 104 Hz, decimation 23; track with struct('dt', 1/fsA)
%
% See also: rheome.filters.bandpass (spectral gain on the joint (lambda,omega) axis -- a different
%           object: that one weights eigenbasis coefficients, this one filters a time series),
%           rheome.detect.track
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct; end
    if ~isfield(opts,'oversample') || isempty(opts.oversample), opts.oversample = 8;    end
    if ~isfield(opts,'transition') || isempty(opts.transition), opts.transition = 0.25; end
    if ~isfield(opts,'detrend')    || isempty(opts.detrend),    opts.detrend    = true; end

    band = double(band(:)');
    if numel(band) ~= 2 || band(1) <= 0 || band(2) <= band(1)
        error('filters:bandlimit:band', 'band must be [f1 f2] with 0 < f1 < f2.');
    end
    f1 = band(1);  f2 = band(2);
    if f2 >= fs/2
        error('filters:bandlimit:nyquist', 'f2 (%.3g Hz) must be below Nyquist (%.3g Hz).', f2, fs/2);
    end

    X = double(X);
    [nRow, nTime] = size(X); %#ok<ASGLU>
    rowMean = mean(X, 2);
    if opts.detrend, X = X - rowMean; end

    % ---- 1. declared gain on the one-sided |f| axis ----
    freq = (0:nTime-1) / nTime * fs;
    fAbs = min(freq, fs - freq);                    % fold to |f|, so the gain is even => real output
    w    = max(opts.transition * (f2 - f1), eps);   % roll-off half-width (Hz)
    gain = zeros(1, nTime);
    gain(fAbs >= f1 & fAbs <= f2) = 1;
    rise = fAbs > (f1 - w) & fAbs < f1;
    fall = fAbs > f2       & fAbs < (f2 + w);
    gain(rise) = 0.5 * (1 - cos(pi * (fAbs(rise) - (f1 - w)) / w));
    gain(fall) = 0.5 * (1 + cos(pi * (fAbs(fall) -  f2     ) / w));

    % ---- 2. apply (zero-phase: a real, even gain multiplies the spectrum) ----
    Xf = real(ifft(fft(X, [], 2) .* gain, [], 2));

    % ---- 3. decimate to a band-matched rate. Alias-free: content above f2+w is gone. ----
    if isfield(opts,'decimate') && ~isempty(opts.decimate)
        dec = max(1, round(opts.decimate));
    else
        fsTarget = opts.oversample * f2;
        dec      = max(1, floor(fs / fsTarget));
    end
    Y     = Xf(:, 1:dec:end);
    fsOut = fs / dec;

    % ---- 4. receipts ----
    vIn  = sum(X(:).^2);
    vOut = sum(Xf(:).^2);
    info = struct('band', band, 'oversample', opts.oversample, 'transition', opts.transition, ...
                  'decimation', dec, 'fsIn', fs, 'fsOut', fsOut, ...
                  'gain', gain, 'freq', freq, ...
                  'passbandFraction', vOut / max(vIn, eps));
end

% Author: Diellor Basha, 2026
