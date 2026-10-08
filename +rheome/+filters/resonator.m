function G = resonator(lambda, omega, f0, Q)
% FILTERS.RESONATOR  Damped harmonic resonator (eigenmode x frequency, 'js' domain).
%
%   G = rheome.filters.resonator(lambda, omega, f0, Q)
%
% g(omega) = f0^2 / (f0^2 - omega^2 + i*omega*f0/Q). The transfer function of a damped
% harmonic oscillator resonant at f0 Hz with quality Q -- lambda-independent (same at
% every spatial scale). Its impulse response is a decaying f0-Hz oscillation. COMPLEX
% (Hermitian in signed frequency). Port of bst_eigfilter_design_resonator.
%
% INPUTS:
%   lambda [K x 1] (used only for the output row count)
%   omega  [1 x N] one-sided frequency axis (Hz) from rheome.filters.jspectrum
%   f0     resonant frequency (Hz, default 10)   Q quality factor (default 6)
% OUTPUT:
%   G      [K x N] complex gains (constant across rows)
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(f0), f0 = 10; end
    if nargin < 4 || isempty(Q),  Q  = 6;  end
    f0 = max(f0, eps);  Q = max(Q, eps);
    ws = i_signed(omega);
    H  = (f0.^2) ./ (f0.^2 - ws.^2 + 1i .* ws .* (f0/Q));   % [1 x N]
    G  = ones(numel(lambda), 1) * H;
end

function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w) * (w(2) - w(1));  ws = w - Fs .* (w >= Fs/2);
end

% Author: Diellor Basha, 2026
