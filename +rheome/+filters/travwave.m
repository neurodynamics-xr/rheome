function G = travwave(lambda, omega, c, width)
% FILTERS.TRAVWAVE  Joint traveling-wave filter: energy on the dispersion ridge.
%
%   G = rheome.filters.travwave(lambda, omega, c, width)
%
% A non-separable joint filter g(lambda,omega) selecting the (lambda,omega) pairs that
% obey a linear dispersion relation -- i.e. waves propagating on the cortex at phase
% speed c. Energy concentrates on the ridge  f = c*sqrt(lambda)/(2*pi) Hz, with Gaussian
% width. Faithful port of bst_eigfilter_design_travwave (band-pass, even in freq -> real).
%
% INPUTS:
%   lambda [K x 1] eigenvalues (dbasis.Lambda)
%   omega  [1 x N] one-sided frequency axis (Hz) from rheome.filters.jspectrum
%   c      phase speed (m/s), default 1
%   width  ridge width (Hz), default 2
%
% OUTPUT:
%   G      [K x N] joint gains (multiply the joint spectrum elementwise)
%
% See also: rheome.filters.jspectrum, rheome.filters.stmatern, rheome.filters.bandpass
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(c),     c = 1;     end
    if nargin < 4 || isempty(width), width = 2; end
    ws = i_signed(omega);
    fr = c .* sqrt(max(double(lambda(:)), 0)) ./ (2*pi);   % [K x 1] ridge frequency
    G  = exp( -((abs(ws) - fr).^2) / (2 * max(width,eps)^2) );
end

function ws = i_signed(w)
    w = double(w(:)');  Fs = numel(w) * (w(2) - w(1));  ws = w - Fs .* (w >= Fs/2);
end

% Author: Diellor Basha, 2026
