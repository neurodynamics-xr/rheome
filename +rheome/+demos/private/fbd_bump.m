function u = fbd_bump(S, seedIdx, sigma)
% FBD_BUMP  A Gaussian bump of known GEODESIC width on the sphere.
%
%   u = fbd_bump(S, seedIdx, sigma)      sigma in metres
%
%   u(v) = exp( -d(v, seed)^2 / (2*sigma^2) ),   d = R * acos(v.s / R^2)
%
% ⭐ WHY THIS IS THE RIGHT TEST FEATURE. Its spectrum is exp(-lambda*sigma^2/2), so
% against a mexhat member t*lambda*exp(-t*lambda) the response is t/(t + sigma^2/2)^2,
% maximised at t = sigma^2/2 -- i.e. the matched member has width sqrt(2t) = sigma
% EXACTLY. The demo can therefore assert a recovered width against the planted one with
% no fitted calibration.
%
% ⚠ NOT the same calibration as rheome.filters.frame's Gamma = 1/sqrt(2). That applies to a
% VORTEX, whose vorticity spectrum carries an extra factor of lambda, moving the maximum
% to t = sigma^2/4.
%
% Author: Diellor Basha, 2026

    s = S.V(seedIdx, :).';
    ct = min(1, max(-1, (S.V * s) / S.R^2));
    d  = S.R * acos(ct);
    u  = exp(-d.^2 / (2*sigma^2));
end

% Author: Diellor Basha, 2026
