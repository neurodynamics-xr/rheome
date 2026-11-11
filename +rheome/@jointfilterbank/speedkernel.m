function K = speedkernel(c, d)
% JOINTFILTERBANK.SPEEDKERNEL  A speed-selective joint kernel.
%
%   K = rheome.jointfilterbank.speedkernel(c)        % d = 0.1 * c, a 10% band
%   K = rheome.jointfilterbank.speedkernel(c, d)
%
%   K(lambda, f) = exp( -(2*pi*f - c*sqrt(lambda)).^2 / (2*d^2) )
%
% A structure travelling at speed c satisfies omega = 2*pi*f = c*sqrt(lambda), so this selects a
% DIAGONAL BAND in the (sqrt(lambda), f) plane. This is the canonical non-separable
% kernel -- it cannot be written as psi_g(lambda)*psi_t(omega), which is exactly why the
% bank stores factors and evaluates lazily rather than assuming separability.
%
% UNITS. The kernel's argument f is the SIGNED frequency in HERTZ, like every (lambda, f) kernel; the
% angular frequency 2*pi*f is formed here, so c stays a speed in length/s (lambda in 1/length^2). d is
% the half-width of the ridge in rad/s (an ANGULAR width, the one angular quantity of the interface).
% The ridge is at +f: a signed ridge (one sense of travel) is a kernel of f's sign, see rheome.filters.travwave.
% Use dispersion() to FIT c from data rather than assuming it.
%
% See also: dispersion, spectralResponse
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(d), d = 0.1 * abs(c); end
    if ~(d > 0)
        error('jointfilterbank:speedkernel', 'width d must be positive, got %g.', d);
    end
    K = @(l, f) exp(-(2*pi*f - c*sqrt(double(l))).^2 / (2*d^2));
end

% Author: Diellor Basha, 2026
