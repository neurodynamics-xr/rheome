function G = diffusion(lambda, t, tau, lmax)
% FILTERS.DIFFUSION  Physical-time heat kernel (eigenmode x time, 'ts' domain).
%
%   G = rheome.filters.diffusion(lambda, t, tau, lmax)
%
% g(lambda,t) = exp(-tau * (lambda/lmax) * |t|). The Green's function of diffusion as a
% function of elapsed time t: applied to a delta it spreads into a widening positive
% bump. Port of bst_eigfilter_design_diffusion (separable, low-pass).
%
% INPUTS:
%   lambda [K x 1] eigenvalues       t [1 x N] time-lag axis
%   tau    diffusion rate (default 1)   lmax spectral normalizer (default max(lambda))
% OUTPUT:
%   G      [K x N] gains g(lambda_k, t_n)
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(tau),  tau  = 1;            end
    if nargin < 4 || isempty(lmax), lmax = max(lambda);  end
    sc = 1 / max(lmax, eps);
    G  = exp( -tau * sc * double(lambda(:)) .* abs(double(t(:)')) );
end

% Author: Diellor Basha, 2026
