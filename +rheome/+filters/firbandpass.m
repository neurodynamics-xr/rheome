function [y, spec] = firbandpass(x, Fs, band, opts)
% FILTERS.FIRBANDPASS  Linear-phase FIR bandpass with a DECLARED transition band.
%
%   [y, spec] = rheome.filters.firbandpass(x, Fs, band)
%   [y, spec] = rheome.filters.firbandpass(x, Fs, band, opts)
%   [~, spec] = rheome.filters.firbandpass([], Fs, band)      % design only, no data
%
% Follows the design used by Brainstorm's bst_bandpass_hfilter: order and window from kaiserord /
% kaiser, coefficients from fir1, applied with fftfilt (or filtfilt for the zero-phase variant).
% Ripple, stopband attenuation and transition band are PARAMETERS, so the filter is declared
% rather than implied.
%
% WHY NOT JUST TRUNCATE THE FFT. Keeping a set of bins and discarding the rest IS a bandpass -- with
% a RECTANGULAR transfer function. Two things follow, and both are avoidable:
%   * a rectangular band is a sinc in time, decaying as 1/t, so its sidelobes are long; and
%   * a whole-record FFT assumes the signal is periodic, so x(T) ~= x(1) puts a discontinuity at
%     the seam whose Gibbs oscillation the retained band picks up.
% An FIR has FINITE support, so edge effects are confined to the filter length instead of spread
% across the record, and mirroring or filtfilt handles what remains.
%
% USE IT BEFORE rheome.flow.joint. Filter properly in time first; the joint-domain truncation is then
% nearly lossless, because everything outside the band has already been pushed into the stopband,
% and it still buys the full cost saving (the GEMM runs over the retained bins either way).
%
% INPUTS:
%   x     [nChan x nTime] signal, or [] to design only
%   Fs    sampling rate (Hz)
%   band  [f1 f2] passband (Hz); f1 = 0 for low-pass, f2 = 0 or Inf for high-pass
%   opts  .tranBand  transition width (Hz). Default: max(0.1*f1, 0.5) at the low edge
%         .ripple    passband ripple (dB), default 0.05
%         .atten     stopband attenuation (dB), default 60
%         .mode      'fftfilt' (linear phase, group delay removed) | 'filtfilt' (zero phase,
%                    forward-backward -- squares the magnitude response, so the effective
%                    attenuation doubles and the effective passband narrows slightly)
%         .mirror    mirror-pad by the filter length before filtering (default true)
%
% OUTPUT:
%   y     filtered signal, same size as x
%   spec  .b coefficients  .order  .tranBand  .ripple  .atten  .mode  .edgeSamples
%         .edgeSamples is the number of samples at each end that the filter's own support touches
%         -- data within it is not trustworthy and should be guarded or discarded.
%
% See also: rheome.flow.joint, rheome.filters.bandlimit
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'ripple') || isempty(opts.ripple), opts.ripple = 0.05; end
    if ~isfield(opts,'atten')  || isempty(opts.atten),  opts.atten  = 60;   end
    if ~isfield(opts,'mode')   || isempty(opts.mode),   opts.mode   = 'fftfilt'; end
    if ~isfield(opts,'mirror') || isempty(opts.mirror), opts.mirror = true; end

    f1 = band(1);  f2 = band(2);
    if ~isfield(opts,'tranBand') || isempty(opts.tranBand)
        opts.tranBand = max(0.1*max(f1, eps), 0.5);
    end
    tb = opts.tranBand;
    nyq = Fs/2;

    % ---- design: kaiserord gives the order and beta for the requested ripple/attenuation ----
    dev = [10^(-opts.atten/20), (10^(opts.ripple/20)-1)/(10^(opts.ripple/20)+1), 10^(-opts.atten/20)];
    fcuts = [f1-tb, f1, f2, f2+tb];
    if fcuts(1) <= 0 || fcuts(4) >= nyq
        error('filters:firbandpass:band', ...
            'transition band %.2f Hz does not fit: [%g %g] Hz at Fs = %g.', tb, f1, f2, Fs);
    end
    [n, Wn, beta, ftype] = kaiserord(fcuts, [0 1 0], dev, Fs);
    n = n + rem(n,2);                                  % even order -> odd length, integer delay
    b = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');

    spec.b = b;  spec.order = n;  spec.tranBand = tb;
    spec.ripple = opts.ripple;  spec.atten = opts.atten;  spec.mode = opts.mode;
    spec.edgeSamples = ceil(n/2);
    if strcmpi(opts.mode,'filtfilt'), spec.edgeSamples = n; end

    if isempty(x), y = []; return; end

    nPad = 0;
    if opts.mirror
        nPad = min(n, size(x,2)-1);
        x = [x(:, nPad+1:-1:2), x, x(:, end-1:-1:end-nPad)];
    end

    switch lower(opts.mode)
        case 'fftfilt'
            y = fftfilt(b, x.').';                     % linear phase
            y = [y(:, n/2+1:end), zeros(size(y,1), n/2)];   % remove the group delay
        case 'filtfilt'
            y = filtfilt(b, 1, x.').';                 % zero phase, forward + backward
        otherwise
            error('filters:firbandpass:mode', 'unknown mode ''%s''.', opts.mode);
    end
    if nPad > 0
        y = y(:, nPad+1:end-nPad);
    end
end

% Author: Diellor Basha, 2026
