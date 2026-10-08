function G = wave(lambda, t, alpha, lmax)
% FILTERS.WAVE  Wave propagation kernel (eigenmode x time, 'ts' domain).
%
%   G = rheome.filters.wave(lambda, t, alpha, lmax)
%
% g(lambda,t) = cos( alpha * t * acos(1 - 2*lambda/lmax) ). The d'Alembert wave
% propagator: each mode oscillates at frequency omega_k = alpha*acos(1-2*lambda_k/lmax),
% so a delta launches an expanding wavefront at speed ~alpha. Port of
% bst_eigfilter_design_wave (non-separable, the dispersion is cortex geometry).
%
% INPUTS:
%   lambda [K x 1]   t [1 x N] time-lag axis
%   alpha  wave speed (default 1)   lmax spectral normalizer (default max(lambda))
% OUTPUT:
%   G      [K x N]
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(alpha), alpha = 1;           end
    if nargin < 4 || isempty(lmax),  lmax  = max(lambda); end
    sc = 2 / max(lmax, eps);
    th = acos( max(min(1 - sc*double(lambda(:)), 1), -1) );   % [K x 1] per-mode angular freq
    G  = cos( alpha * (th * double(t(:)')) );
end

% Author: Diellor Basha, 2026
