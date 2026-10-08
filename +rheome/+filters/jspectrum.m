function [chat, f] = jspectrum(c, sfreq, NFFT)
% FILTERS.JSPECTRUM  Joint eigenmode-frequency spectrum of a mode-coefficient series.
%
%   [chat, f] = rheome.filters.jspectrum(c, sfreq)
%   [chat, f] = rheome.filters.jspectrum(c, sfreq, NFFT)
%
% The mode coefficients c(t) = Phi'*M*U are the cortex (manifold) Fourier transform of
% a time-vertex field. Their unitary temporal FFT completes the JOINT time-vertex
% transform, giving coefficients on the (eigenvalue lambda, temporal-frequency omega)
% grid:  chat = fft(c, NFFT, 2) / sqrt(NFFT). Same normalization as manifold_jft, so
% joint energy equals the cortex (M-weighted) energy.
%
% INPUTS:
%   c      [nModes x nT] mode-coefficient time series (e.g. rheome.source.dirac out.c)
%   sfreq  sampling frequency (Hz)
%   NFFT   temporal bins (default nT)
%
% OUTPUTS:
%   chat   [nModes x NFFT] joint spectral coefficients
%   f      [1 x NFFT] one-sided temporal-frequency axis (Hz), spacing sfreq/NFFT.
%          This is the axis the js kernels (rheome.filters.travwave / stmatern / bandpass)
%          expect; they wrap f >= sfreq/2 to negative internally.
%
% See also: rheome.filters.ijspectrum, rheome.filters.travwave, rheome.filters.stmatern, rheome.show.spectrum2d
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(NFFT), NFFT = size(c, 2); end
    chat = fft(c, NFFT, 2) / sqrt(NFFT);
    f    = (0:NFFT-1) * (sfreq / NFFT);
end

% Author: Diellor Basha, 2026
