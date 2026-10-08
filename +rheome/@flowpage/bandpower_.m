function P = bandpower_(obj)
% BANDPOWER_  Total in-band vorticity energy per frame.
%
%   P = bandpower_(fp)               -> [1 x nT]
%
% The independent variable the scale descriptors are read against: does the spatial scale of
% the flow change with how much band power there is. Trailing underscore because `bandpower`
% is a Signal Processing Toolbox function.
%
% Author: Diellor Basha, 2026

    P = sum(abs(double(obj.Coefficients)).^2, 1);
end

% Author: Diellor Basha, 2026
