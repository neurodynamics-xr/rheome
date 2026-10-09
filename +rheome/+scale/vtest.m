function [u, p] = vtest(th, mu)
% SCALE.VTEST  The V-test (modified Rayleigh test) for circular uniformity against a known mean direction.
%
%   [u, p] = rheome.scale.vtest(theta)        % mu = 0
%   [u, p] = rheome.scale.vtest(theta, mu)    % angles in radians
%
% V = sum cos(theta - mu), u = V sqrt(2/n); under uniformity u is standard normal for large n, and p is the
% one-sided upper tail (Zar, Biostatistical Analysis, 27.5). NaN angles are dropped; n = 0 gives NaN.
%
% See also: rheome.scale.measure_slowosc, rheome.scale.groupslowosc
%
% Author: Diellor Basha, 2026

    if nargin < 2, mu = 0; end
    th = th(isfinite(th));  n = numel(th);
    if n == 0, u = NaN; p = NaN; return, end
    u = sum(cos(th - mu)) * sqrt(2 / n);
    p = 0.5 * erfc(u / sqrt(2));
end

% Author: Diellor Basha, 2026
