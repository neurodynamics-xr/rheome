function [Yt, info] = timetiles(Y, rate, P, opts)
% DETECT.TIMETILES  Means of a frame sequence over the atlas's time tiles (edges at whole multiples of P).
%
%   [Yt, info] = rheome.detect.timetiles(Y, rate, P)             % frame 1 at t = 0
%   [Yt, info] = rheome.detect.timetiles(Y, rate, P, T0=12.3)    % frame 1 at t = T0 s
%
% Y is [n x nFrames] at `rate` frames/s. Tile k covers [kP, (k+1)P) of the recording's time, so tiles of every
% analysis line up with the atlas's (rheome.detect.atlastime); only complete tiles are returned (a tile short
% of less than half a frame counts as complete). Each frame belongs to the tile its start falls in.
% ⚠ P is generally not a whole number of frames (alpha: 49.4 frames at 300/s): tiles hold 49 or 50 frames.
%
% OUTPUT  Yt [n x nTiles] tile means (double);  info.k absolute tile numbers, info.mid the frame nearest each
%         tile's centre, info.rate = 1/P
%
% See also: rheome.detect.atlastime, rheome.detect.tilepath
%
% Author: Diellor Basha, 2026

    arguments
        Y {mustBeNumericOrLogical}
        rate (1,1) double {mustBePositive}
        P (1,1) double {mustBePositive}
        opts.T0 (1,1) double = 0
    end
    nF = size(Y, 2);  t = opts.T0 + (0:nF-1)/rate;
    k = floor(t/P + 1e-9);
    kk = ceil(opts.T0/P - 0.5/(rate*P)) : floor((t(end) + 1.5/rate)/P) - 1;
    if isempty(kk), error('detect:timetiles:short', 'No complete tile of %.4g s in %.4g s.', P, nF/rate); end
    in = find(k >= kk(1) & k <= kk(end));
    A = sparse(in, k(in) - kk(1) + 1, 1, nF, numel(kk));
    Yt = (double(Y) * A) ./ full(sum(A, 1));
    info = struct('k', kk, 'mid', min(max(round(((kk + 0.5)*P - opts.T0)*rate) + 1, 1), nF), 'rate', 1/P);
end

% Author: Diellor Basha, 2026
