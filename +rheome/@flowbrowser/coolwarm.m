function c = coolwarm(n)
% COOLWARM  A perceptually balanced diverging colormap (Moreland cool-to-warm).
%
%   c = rheome.flowbrowser.coolwarm(n)
%
% ⚠ WHY NOT A NAIVE BLUE-WHITE-RED. Ramping straight from [0 0 1] to white to [1 0 0] gives
% the two arms UNEQUAL luminance: pure blue is much darker than pure red, so counter-clockwise
% and clockwise rotation of the same magnitude do not look the same magnitude, and the eye
% reads a handedness bias that is not in the data. This map's endpoints are chosen to have
% matched lightness either side of a neutral centre, so the two handednesses are visually
% symmetric.
%
% Zero maps to the neutral centre, which is what keeps the CCW/CW boundary readable.
%
% Endpoints after Moreland, "Diverging Color Maps for Scientific Visualization" (2009).
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(n), n = 256; end
    lo  = [0.230 0.299 0.754];      % cool
    mid = [0.865 0.865 0.865];      % neutral, NOT pure white -- keeps both arms symmetric
    hi  = [0.706 0.016 0.150];      % warm
    m   = floor(n/2);
    t1  = linspace(0, 1, m).';
    t2  = linspace(0, 1, n - m).';
    c   = [ lo  + t1.*(mid - lo) ; mid + t2.*(hi - mid) ];
    c   = min(max(c, 0), 1);
end

% Author: Diellor Basha, 2026
