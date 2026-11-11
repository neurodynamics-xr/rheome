function G = stmatern(lambda, f, kappa, nu, v, half)
% FILTERS.STMATERN  Spatiotemporal Whittle-Matern joint filter (the 1/f background).
%
%   G = rheome.filters.stmatern(lambda, f, kappa, nu, v)
%   G = rheome.filters.stmatern(lambda, f, kappa, nu, v, half)
%
% A non-separable, low-pass joint prior -- the aperiodic 1/f-like spatiotemporal
% background spectral density:
%       g(lambda,f) = ( kappa^2 + lambda + (f/v)^2 )^(-nu),
% f the SIGNED frequency in HERTZ.
% Faithful port of bst_eigfilter_design_stmatern (prior-admissible; even in freq).
%
% INPUTS:
%   lambda [K x 1] eigenvalues (dbasis.Lambda)
%   f      [1 x N] frequency axis (Hz) from rheome.filters.jspectrum (full FFT grid unless half = 'positive')
%   kappa  inverse spatial correlation length (default 126)
%   nu     smoothness (default 1.5)
%   v      space-time coupling speed (length*Hz, default 1)
%   half   'full' (default) | 'positive': see rheome.filters.signedfrequencies
%
% OUTPUT:
%   G      [K x N] joint gains
%
% See also: rheome.filters.jspectrum, rheome.filters.travwave
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(kappa), kappa = 126; end
    if nargin < 4 || isempty(nu),    nu    = 1.5; end
    if nargin < 5 || isempty(v),     v     = 1;   end
    if nargin < 6 || isempty(half), half = 'full'; end
    ws   = rheome.filters.signedfrequencies(f, half);
    base = (kappa.^2 + double(lambda(:))) * ones(1, numel(ws)) + ones(numel(lambda),1) * (ws./v).^2;
    G    = base .^ (-nu);
end


% Author: Diellor Basha, 2026
