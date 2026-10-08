function u = fbd_sectoral(S, l, phase)
% FBD_SECTORAL  A sectoral spherical harmonic, sin^l(theta) * cos(l*phi + phase).
%
%   u = fbd_sectoral(S, l, phase)
%
% Sectoral (m = l) harmonics are used because they are exactly writable in closed form
% and rotate rigidly: advancing phase by -l*Omega*dt rotates the pattern at Omega.
% A degree-l harmonic sits at lambda = l(l+1)/R^2, which is what makes the ridge in
% rheome.demos.filterbank_joint analytic.
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(phase), phase = 0; end
    x = S.V(:,1);  y = S.V(:,2);  z = S.V(:,3);
    th = acos(min(1, max(-1, z / S.R)));
    ph = atan2(y, x);
    u  = (sin(th).^l) .* cos(l*ph + phase);
end

% Author: Diellor Basha, 2026
