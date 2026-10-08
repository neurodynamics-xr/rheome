function w = wavelengths(obj)
% WAVELENGTHS  2*pi/k per member. The centerPeriods analogue.
%   w = wavelengths(gfb)     [1 x M], in the operator's length units
%
% See also: centerWavenumbers, widths
%
% Author: Diellor Basha, 2026
    w = 2*pi ./ max(centerWavenumbers(obj), eps);
end

% Author: Diellor Basha, 2026
