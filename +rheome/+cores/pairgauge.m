function [fa, fb, u, theta] = pairgauge(fa, fb, u)
% CORES.PAIRGAUGE  The canonical gauge of a pair: rotate (fa, fb) in their plane so that fb = 0 and fa > 0 on the strongest sensor.
%
%   [fa, fb, u, theta] = rheome.cores.pairgauge(fa, fb, u)
%
% j = first argmax(fa^2 + fb^2); theta = atan2(fb(j), fa(j)); fa' = c fa + s fb, fb' = -s fa + c fb, u' = u e^{-i theta}.
% The orientation (sign of fb, q after p) is kept. This is the comparison gauge of nxr's vortex-pair goldens (the golden script
% applies it; nxr's pair does not), so two ports of pair can be compared. u is optional ([] to skip).
%
% See also: rheome.cores.pair
%
% Author: Diellor Basha, 2026
    arguments
        fa (1,:) double
        fb (1,:) double
        u {mustBeNumeric} = []
    end
    [~, j] = max(fa.^2 + fb.^2);
    theta = atan2(fb(j), fa(j)); c = cos(theta); s = sin(theta);
    fa0 = fa;
    fa = c * fa0 + s * fb;
    fb = -s * fa0 + c * fb;
    u = u * exp(-1i * theta);
end
