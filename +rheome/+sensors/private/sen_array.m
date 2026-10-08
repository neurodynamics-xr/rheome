function arr = sen_array(name, kind, P, dim, labels)
% SEN_ARRAY  Assemble and validate the sensor-array struct.
%
%   arr = sen_array(name, kind, P, dim, labels)
%
% ONE definition of the array contract, so all four constructors agree on what an array is
% and measure pitch and aperture the same way.
%
% INPUTS:
%   name    char, display name
%   kind    'linear' | 'grid' | 'meg' | 'eeg'
%   P       [N x 3] coordinates, metres
%   dim     intrinsic dimension, 1 or 2
%   labels  {1 x N} cellstr channel names
% OUTPUT:
%   arr  .Name .Kind .Pos .Dim .Pitch .Aperture .Labels .nCh
%
% See also: rheome.sensors.linear, rheome.sensors.grid, rheome.sensors.meg, rheome.sensors.eeg
%
% Author: Diellor Basha, 2026

    if size(P, 2) ~= 3
        error('sensors:array:pos', 'Pos must be [N x 3], got %s.', mat2str(size(P)));
    end
    n = size(P, 1);
    if n < 2
        error('sensors:array:count', 'An array needs at least 2 sensors, got %d.', n);
    end
    if ~any(dim == [1 2])
        error('sensors:array:dim', 'Dim must be 1 or 2, got %g.', dim);
    end
    labels = labels(:).';
    if numel(labels) ~= n
        error('sensors:array:labels', '%d labels for %d sensors.', numel(labels), n);
    end

    [pitch, aperture] = sen_metrics(P);
    arr = struct('Name', char(name), 'Kind', char(kind), 'Pos', P, 'Dim', dim, ...
                 'Pitch', pitch, 'Aperture', aperture, 'Labels', {labels}, 'nCh', n);
end

% Author: Diellor Basha, 2026
