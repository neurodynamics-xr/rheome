function c = ijspectrum(chat, nKeep)
% FILTERS.IJSPECTRUM  Inverse joint transform: joint coefficients -> mode series.
%
%   c = rheome.filters.ijspectrum(chat)
%   c = rheome.filters.ijspectrum(chat, nKeep)
%
% Inverts rheome.filters.jspectrum along time:  c = ifft(chat, [], 2) * sqrt(NFFT), then keeps
% the first nKeep samples (to un-pad). Returns the real part -- for a real field filtered
% by an even-in-frequency kernel the imaginary part is numerical noise. Same convention
% as manifold_ijft (temporal half).
%
% INPUTS:
%   chat   [nModes x NFFT] (filtered) joint spectral coefficients
%   nKeep  keep the first nKeep time samples (default NFFT)
%
% OUTPUT:
%   c      [nModes x nKeep] mode-coefficient time series
%
% See also: rheome.filters.jspectrum, rheome.forward.reconstruct
%
% Author: Diellor Basha, 2026

    NFFT = size(chat, 2);
    if nargin < 2 || isempty(nKeep), nKeep = NFFT; end
    C = ifft(chat, [], 2) * sqrt(NFFT);
    c = real(C(:, 1:min(nKeep, NFFT)));
end

% Author: Diellor Basha, 2026
