function g = bandpass(f, f0, bw, half)
% FILTERS.BANDPASS  Gaussian temporal band-pass gain on the joint-frequency axis.
%
%   g = rheome.filters.bandpass(f, f0, bw)
%   g = rheome.filters.bandpass(f, f0, bw, half)
%
% A separable (frequency-only) factor for building joint filters g(lambda,f) =
% g_space(lambda) .* g_bandpass(f). Centred on +/- f0 Hz with Gaussian half-width bw:
%       g(f) = exp( -(|f_signed| - f0)^2 / (2 bw^2) ),
% all in HERTZ.
% Even in frequency, so the resulting time-domain filter is real.
%
% INPUTS:
%   f      [1 x N] frequency axis (Hz) from rheome.filters.jspectrum (full FFT grid unless half = 'positive')
%   f0     centre frequency (Hz), e.g. 10.5 for alpha
%   bw     Gaussian half-width (Hz), e.g. 1.5 for a ~9-12 Hz band
%   half   'full' (default) | 'positive': see rheome.filters.signedfrequencies
%
% OUTPUT:
%   g      [1 x N] gain per frequency bin
%
% See also: rheome.filters.jspectrum, rheome.filters.travwave, rheome.filters.stmatern
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(half), half = 'full'; end
    ws = rheome.filters.signedfrequencies(f, half);
    g  = exp( -((abs(ws) - f0).^2) / (2 * max(bw, eps)^2) );
end


% Author: Diellor Basha, 2026
