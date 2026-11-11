function G = dampedwave(lambda, t, c, beta, varargin)
% FILTERS.DAMPEDWAVE  Damped wave kernel (eigenmode x time lag, 'ts' domain).
%
%   G = rheome.filters.dampedwave(lambda, t, c, beta)                       % continuum (default)
%   G = rheome.filters.dampedwave(lambda, t, c, beta, 'Variant','discrete', 'LambdaMax', lmax)
%
% CONTINUUM (the definition):  g(lambda,t) = exp(-beta*|t|) * cos( t * c * sqrt(lambda) ).
% A wave that propagates at speed c (omega_k = c*sqrt(lambda_k), rad/s) while every mode decays
% at the rate beta (1/s): an expanding wavefront that rings down. beta is a given amplitude
% decay rate, not derived from the equation (the exact solution of u_tt + 2*beta*u_t = -c^2 L u
% has a shifted frequency sqrt(c^2*lambda - beta^2)). In PHYSICAL UNITS, not normalised by lambda_max.
%
% DISCRETE (named variant, never the default): the Chebyshev-angle form
% exp(-beta*|t|) * cos( c * t * acos(1 - 2*lambda/lmax) ), as in rheome.filters.wave.
% Port of bst_eigfilter_design_dampedwave.
%
% A (lambda,t) kernel acts on a time LAG; it cannot multiply a joint (lambda,f) spectrum.
%
% INPUTS:
%   lambda [K x 1] (1/length^2)   t [1 x N] time-lag axis (s)
%   c      wave speed (length/s, default 1)   beta damping (1/s, default 0.1)
% NAME-VALUE:
%   'Variant'    'continuum' (default) | 'discrete'
%   'LambdaMax'  normaliser of the discrete variant (default max(lambda))
% OUTPUT:
%   G      [K x N]
%
% See also: rheome.filters.wave, rheome.filters.kleingordon
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(c),    c    = 1;   end
    if nargin < 4 || isempty(beta), beta = 0.1; end
    tt = double(t(:)');
    W  = rheome.filters.wave(lambda, tt, c, varargin{:});
    G  = W .* exp(-beta*abs(tt));
end

% Author: Diellor Basha, 2026
