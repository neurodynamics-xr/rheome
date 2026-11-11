function G = resonator(lambda, f, f0, Q, half)
% FILTERS.RESONATOR  Damped harmonic resonator (eigenmode x frequency, 'js' domain).
%
%   G = rheome.filters.resonator(lambda, f, f0, Q)
%   G = rheome.filters.resonator(lambda, f, f0, Q, half)
%
% g(f) = f0^2 / (f0^2 - f^2 + i*f*f0/Q), f the SIGNED frequency in HERTZ. The transfer function of a damped
% harmonic oscillator resonant at f0 Hz with quality Q -- lambda-independent (same at
% every spatial scale). Its impulse response is a decaying f0-Hz oscillation. COMPLEX
% (Hermitian in signed frequency). Port of bst_eigfilter_design_resonator.
%
% INPUTS:
%   lambda [K x 1] (used only for the output row count)
%   f      [1 x N] frequency axis (Hz) from rheome.filters.jspectrum (full FFT grid unless half = 'positive')
%   f0     resonant frequency (Hz, default 10)   Q quality factor (default 6)
%   half   'full' (default) | 'positive': see rheome.filters.signedfrequencies
% OUTPUT:
%   G      [K x N] complex gains (constant across rows)
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(f0), f0 = 10; end
    if nargin < 4 || isempty(Q),  Q  = 6;  end
    f0 = max(f0, eps);  Q = max(Q, eps);
    if nargin < 5 || isempty(half), half = 'full'; end
    ws = rheome.filters.signedfrequencies(f, half);
    H  = (f0.^2) ./ (f0.^2 - ws.^2 + 1i .* ws .* (f0/Q));   % [1 x N]
    G  = ones(numel(lambda), 1) * H;
end


% Author: Diellor Basha, 2026
