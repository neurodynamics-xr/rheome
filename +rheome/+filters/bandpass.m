function g = bandpass(omega, f0, bw)
% FILTERS.BANDPASS  Gaussian temporal band-pass gain on the joint-frequency axis.
%
%   g = rheome.filters.bandpass(omega, f0, bw)
%
% A separable (frequency-only) factor for building joint filters g(lambda,omega) =
% g_space(lambda) .* g_bandpass(omega). Centred on +/- f0 Hz with Gaussian half-width bw:
%       g(omega) = exp( -(|omega_signed| - f0)^2 / (2 bw^2) ).
% Even in frequency, so the resulting time-domain filter is real.
%
% INPUTS:
%   omega  [1 x N] one-sided frequency axis (Hz) from rheome.filters.jspectrum
%   f0     centre frequency (Hz), e.g. 10.5 for alpha
%   bw     Gaussian half-width (Hz), e.g. 1.5 for a ~9-12 Hz band
%
% OUTPUT:
%   g      [1 x N] gain per frequency bin
%
% See also: rheome.filters.jspectrum, rheome.filters.travwave, rheome.filters.stmatern
%
% Author: Diellor Basha, 2026

    ws = i_signed(omega);
    g  = exp( -((abs(ws) - f0).^2) / (2 * max(bw, eps)^2) );
end

function ws = i_signed(w)
% Map a one-sided [0,Fs) Hz axis to signed frequencies (wrap w >= Fs/2 to negative),
% matching Brainstorm's js-kernel convention (bst_eigfilter_design_*).
    w = double(w(:)');  Fs = numel(w) * (w(2) - w(1));  ws = w - Fs .* (w >= Fs/2);
end

% Author: Diellor Basha, 2026
