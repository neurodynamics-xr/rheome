function G = dampedwave(lambda, t, alpha, beta, lmax)
% FILTERS.DAMPEDWAVE  Damped wave kernel (eigenmode x time, 'ts' domain).
%
%   G = rheome.filters.dampedwave(lambda, t, alpha, beta, lmax)
%
% g(lambda,t) = exp(-beta*|t|) * cos( alpha * t * acos(1 - 2*lambda/lmax) ). A wave that
% propagates (speed alpha) while decaying (rate beta) -- an expanding wavefront that
% rings down. Port of bst_eigfilter_design_dampedwave.
%
% INPUTS:
%   lambda [K x 1]   t [1 x N] time-lag axis
%   alpha  wave speed (default 1)   beta damping (default 0.1)   lmax (default max(lambda))
% OUTPUT:
%   G      [K x N]
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(alpha), alpha = 1;           end
    if nargin < 4 || isempty(beta),  beta  = 0.1;         end
    if nargin < 5 || isempty(lmax),  lmax  = max(lambda); end
    sc = 2 / max(lmax, eps);
    th = acos( max(min(1 - sc*double(lambda(:)), 1), -1) );
    tt = double(t(:)');
    G  = (ones(numel(lambda),1) * exp(-beta*abs(tt))) .* cos( alpha * (th * tt) );
end

% Author: Diellor Basha, 2026
