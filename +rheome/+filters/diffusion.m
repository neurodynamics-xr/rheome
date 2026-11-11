function G = diffusion(lambda, t, D)
% FILTERS.DIFFUSION  Physical-time heat kernel (eigenmode x time lag, 'ts' domain).
%
%   G = rheome.filters.diffusion(lambda, t, D)
%
% g(lambda,t) = exp(-D * lambda * |t|). The Green's function of diffusion u_t = -D L u as a
% function of elapsed time t: applied to a delta it spreads into a widening positive bump.
% PHYSICAL UNITS: lambda in 1/length^2, D (the diffusivity) in length^2/s, t in seconds; the
% exponent is not normalised by lambda_max. (The earlier form exp(-tau*(lambda/lmax)*|t|) had
% a time unit relative to the largest eigenvalue; it is gone.) Port of
% bst_eigfilter_design_diffusion (separable, low-pass).
%
% A (lambda,t) kernel acts on a time LAG; it cannot multiply a joint (lambda,f) spectrum.
%
% INPUTS:
%   lambda [K x 1] eigenvalues (1/length^2)   t [1 x N] time-lag axis (s)
%   D      diffusivity (length^2/s, default 1)
% OUTPUT:
%   G      [K x N] gains g(lambda_k, t_n)
%
% See also: rheome.filters.heat, rheome.filters.schrodinger, rheome.filters.wave
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(D), D = 1; end
    G = exp( -D * double(lambda(:)) .* abs(double(t(:)')) );
end

% Author: Diellor Basha, 2026
