function [pitch, aperture, D] = sen_metrics(P)
% SEN_METRICS  Pitch and aperture of a sensor array, from its coordinates alone.
%
%   [pitch, aperture, D] = sen_metrics(P)
%
% pitch is the MEDIAN NEAREST-NEIGHBOUR distance and aperture is the MAXIMUM PAIRWISE
% distance. The two bound the array's measurement window: wavelengths below 2*pitch alias,
% and above the aperture the array measures a phase gradient rather than a wave.
%
% ⚠ APERTURE IS THE DIAGONAL, NOT THE SIDE. For a square grid the largest available
% baseline runs corner to corner; taking the side understates the ceiling by 41%.
%
% No toolbox: the arrays here are at most a few hundred sensors, so the full distance
% matrix is cheaper than any neighbour-search machinery.
%
% INPUTS:
%   P  [N x 3] sensor coordinates, metres
% OUTPUT:
%   pitch     median nearest-neighbour distance, metres
%   aperture  maximum pairwise distance, metres
%   D         [N x N] distance matrix, zero diagonal
%
% See also: rheome.sensors.window
%
% Author: Diellor Basha, 2026

    n  = size(P, 1);
    G  = P * P';
    d2 = max(diag(G) + diag(G).' - 2*G, 0);
    D  = sqrt(d2);
    D(1:n+1:end) = 0;                       % exact zero diagonal

    if n < 2
        pitch = 0;  aperture = 0;  return;
    end
    Dnn = D;
    Dnn(1:n+1:end) = Inf;                   % exclude self for the nearest-neighbour pass
    pitch    = median(min(Dnn, [], 2));
    aperture = max(D(:));
end

% Author: Diellor Basha, 2026
