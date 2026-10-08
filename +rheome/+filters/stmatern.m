function G = stmatern(lambda, omega, kappa, nu, v)
% FILTERS.STMATERN  Spatiotemporal Whittle-Matern joint filter (the 1/f background).
%
%   G = rheome.filters.stmatern(lambda, omega, kappa, nu, v)
%
% A non-separable, low-pass joint prior -- the aperiodic 1/f-like spatiotemporal
% background spectral density:
%       g(lambda,omega) = ( kappa^2 + lambda + (omega/v)^2 )^(-nu).
% Faithful port of bst_eigfilter_design_stmatern (prior-admissible; even in freq).
%
% INPUTS:
%   lambda [K x 1] eigenvalues (dbasis.Lambda)
%   omega  [1 x N] one-sided frequency axis (Hz) from rheome.filters.jspectrum
%   kappa  inverse spatial correlation length (default 126)
%   nu     smoothness (default 1.5)
%   v      space-time coupling speed m/s (default 1)
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
    ws   = i_signed(omega);
    base = (kappa.^2 + double(lambda(:))) * ones(1, numel(ws)) + ones(numel(lambda),1) * (ws./v).^2;
    G    = base .^ (-nu);
end

function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w) * (w(2) - w(1));  ws = w - Fs .* (w >= Fs/2);
end

% Author: Diellor Basha, 2026
