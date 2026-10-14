function A = atlastime(f)
% DETECT.ATLASTIME  The dynamics atlas's time tile for a frequency: the tower band that holds it and its tile.
%
%   A = rheome.detect.atlastime(10)            % alpha: level -19, P = 0.1648 s, band 6.07-12.14 Hz
%   A = rheome.detect.atlastime(sqrt(8*16))    % MS1's 8-16 Hz band, by its geometric centre: the same tile
%
% The atlas divides time with a tower of cycles: band L covers [1/P, 2/P) Hz and its tiles
% last P = 86 400 s * 2^L, one cycle of the band's lowest frequency; a tile at L is exactly two at L-1. One
% wavelet per octave sits on each band (@timefilterbank at one voice per octave, centred on sqrt(2)/P,
% half-power edges 1/P and 2/P). This returns that tile for the band holding f.
%
% OUTPUT (struct A)
%   .level   L                                       .P       tile length, s
%   .band    [1/P 2/P], Hz (half-power edges)        .fc      sqrt(2)/P, the band's wavelet centre, Hz
%   .rate    1/P, tiles per second
%
% See also: rheome.detect.timetiles, rheome.filters.atlasband, rheome.scale.measure_atlas
%
% Author: Diellor Basha, 2026

    arguments
        f (1,1) double {mustBePositive}
    end
    L = ceil(log2(1/(f*86400)) - 1e-12);              % 1/P <= f < 2/P
    P = 86400 * 2^L;
    A = struct('level', L, 'P', P, 'band', [1 2]/P, 'fc', sqrt(2)/P, 'rate', 1/P);
end

% Author: Diellor Basha, 2026
