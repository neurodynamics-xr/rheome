function U = sample(arr, fieldFcn, t)
% SENSORS.SAMPLE  Evaluate an analytic field at the sensor coordinates.
%
%   U = rheome.sensors.sample(arr, @(P) cos(k*P(:,1)))
%   U = rheome.sensors.sample(arr, @(P, tj) cos(k*P(:,1) - w*tj), tvec)
%
% The array does not host the field -- it SAMPLES one defined in the surrounding space.
% Keeping that distinction in the code is what makes the aperture and alias limits mean
% something: they are properties of the sampling, and the field is indifferent to them.
%
% INPUTS:
%   arr       a rheome.sensors.* array struct
%   fieldFcn  @(P) -> [N x 1], or @(P, tj) -> [N x 1] when t is given
%   t         [1 x nT] times, seconds (optional)
% OUTPUT:
%   U  [N x nT] (or [N x 1] with no t)
%
% See also: rheome.sensors.calibrate, rheome.sensors.grid
%
% Author: Diellor Basha, 2026

    P = arr.Pos;
    if nargin < 3 || isempty(t)
        U = fieldFcn(P);
        U = U(:);
        if numel(U) ~= arr.nCh
            error('sensors:sample:size', ...
                'fieldFcn returned %d values for %d sensors.', numel(U), arr.nCh);
        end
        return;
    end

    t  = double(t(:)).';
    nT = numel(t);
    U  = zeros(arr.nCh, nT);
    for j = 1:nT
        uj = fieldFcn(P, t(j));
        U(:, j) = uj(:);
    end
end

% Author: Diellor Basha, 2026
