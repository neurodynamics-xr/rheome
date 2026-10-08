function Gjs = jtv(G, domain, NFFT)
% FILTERS.JTV  Evaluate a kernel into the joint-spectral (eigenmode x frequency) grid.
%
%   Gjs = rheome.filters.jtv(G, domain [,NFFT])
%
% Puts any joint kernel on the common (lambda, omega) axes so its shape can be shown or
% multiplied against rheome.filters.jspectrum output. An eigen-TIME kernel g(lambda,t) ('ts')
% is FFT'd along time into the joint spectrum; an eigen-FREQUENCY kernel g(lambda,omega)
% ('js') is already there. Port of bst_eigfilter_jtv_evaluate.
%
% INPUTS:
%   G       [K x N] the kernel already evaluated on its native axis:
%             'ts' -> G = g(lambda, t)      (time-lag axis)
%             'js' -> G = g(lambda, omega)  (frequency axis)
%   domain  'ts' | 'js'
%   NFFT    temporal bins (default size(G,2)); the 1/sqrt(NFFT) matches manifold_jft
%
% OUTPUT:
%   Gjs     [K x NFFT] joint-spectral gains (multiply rheome.filters.jspectrum output elementwise)
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(NFFT), NFFT = size(G, 2); end
    switch lower(domain)
        case 'ts', Gjs = fft(G, NFFT, 2) / sqrt(NFFT);   % eigen-time -> joint-spectral
        case 'js', Gjs = G;                              % already eigen-frequency
        otherwise, error('filters:jtv:domain', 'domain must be ''ts'' or ''js''.');
    end
end

% Author: Diellor Basha, 2026
