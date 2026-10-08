function G = kleingordon(lambda, t, alpha, mu, lmax)
% FILTERS.KLEINGORDON  Massive (Klein-Gordon) wave kernel (eigenmode x time, 'ts').
%
%   G = rheome.filters.kleingordon(lambda, t, alpha, mu, lmax)
%
% g(lambda,t) = cos( alpha * t * acos(1 - 2*(lambda+mu)/lmax) ). A wave with a mass term
% mu: even the DC mode oscillates (rest frequency ~ sqrt(mu)), and the dispersion is
% shifted, so wavefronts spread with a frequency floor. Port of
% bst_eigfilter_design_kleingordon.
%
% INPUTS:
%   lambda [K x 1]   t [1 x N] time-lag axis
%   alpha  speed (default 1)   mu mass (default 0.1)   lmax (default max(lambda))
% OUTPUT:
%   G      [K x N]
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(alpha), alpha = 1;           end
    if nargin < 4 || isempty(mu),    mu    = 0.1;         end
    if nargin < 5 || isempty(lmax),  lmax  = max(lambda); end
    sc = 2 / max(lmax, eps);
    th = acos( max(min(1 - sc*(double(lambda(:)) + mu), 1), -1) );
    G  = cos( alpha * (th * double(t(:)')) );
end

% Author: Diellor Basha, 2026
