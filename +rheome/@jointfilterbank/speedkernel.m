function K = speedkernel(c, d)
% JOINTFILTERBANK.SPEEDKERNEL  A speed-selective joint kernel.
%
%   K = rheome.jointfilterbank.speedkernel(c)        % d = 0.1 * c, a 10% band
%   K = rheome.jointfilterbank.speedkernel(c, d)
%
%   K(lambda, omega) = exp( -(omega - c*sqrt(lambda)).^2 / (2*d^2) )
%
% A structure travelling at speed c satisfies omega = c*sqrt(lambda), so this selects a
% DIAGONAL BAND in the (sqrt(lambda), omega) plane. This is the canonical non-separable
% kernel -- it cannot be written as psi_g(lambda)*psi_t(omega), which is exactly why the
% bank stores factors and evaluates lazily rather than assuming separability.
%
% d is the half-width in rad/s of the ridge. Use dispersion() to FIT c from data rather
% than assuming it.
%
% See also: dispersion, spectralResponse
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(d), d = 0.1 * abs(c); end
    if ~(d > 0)
        error('jointfilterbank:speedkernel', 'width d must be positive, got %g.', d);
    end
    K = @(l, w) exp(-(w - c*sqrt(double(l))).^2 / (2*d^2));
end

% Author: Diellor Basha, 2026
