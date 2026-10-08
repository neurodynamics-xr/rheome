function P = tilepeaks(Y, G, opts)
% DETECT.TILEPEAKS  Peaks of a tile-level field: tiles whose value beats every adjacent tile.
%
%   P = rheome.detect.tilepeaks(Y, G)
%   P = rheome.detect.tilepeaks(Y, G, MinHeight=0.5, Relative=true)
%
% The tile twin of rheome.detect.peaks. Y is a field already pooled over tiles (rheome.geom.tilemean), so a peak is
% a whole tile -- many vertices at once -- and "neighbour" means a tile sharing a boundary (G.A > 0),
% not a mesh edge. Ties are broken by tile index, so a flat pair yields one peak.
%
% ⭐ THE PEAK IS A TILE, SO ITS POSITION IS QUANTISED TO THE TILING ON PURPOSE. The resolution is chosen
% by the tree depth, as the time resolution is chosen by the dyadic level, rather than imposed by the
% mesh; a trajectory is then a walk on the tile graph. The same rule as rheome.detect.peaks applies: a tile
% is a peak only if it beats its neighbours, never merely because it is some region's maximum.
%
% INPUTS
%   Y          [nTile x nT] tile values        G  rheome.geom.tiles struct (.A)
%   MinHeight  drop peaks below this (0)       Relative  true: fraction of each frame's maximum (false)
%
% OUTPUT (struct P)
%   .isPeak [nTile x nT] logical     .count [1 x nT]
%   .tiles  {1 x nT} peak tiles per frame, highest first     .values {1 x nT} their values
%
% See also: rheome.detect.tiletrack, rheome.geom.tilemean, rheome.geom.tiles, rheome.detect.peaks
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        G struct
        opts.MinHeight (1,1) double = 0
        opts.Relative (1,1) logical = false
    end
    nTile = size(Y, 1);  nT = size(Y, 2);
    [I, J] = find(G.A);                                   % both directions of each adjacency
    inc = sparse(1:numel(I), I, 1, numel(I), nTile);
    beaten = Y(J,:) > Y(I,:) | (Y(J,:) == Y(I,:) & J < I);
    isPk = (inc' * double(beaten)) == 0;
    thr = opts.MinHeight * ones(1, nT);
    if opts.Relative, thr = opts.MinHeight * max(Y, [], 1); end
    isPk = isPk & Y >= thr;
    P.isPeak = isPk;
    P.count = sum(isPk, 1);
    P.tiles = cell(1, nT);  P.values = cell(1, nT);
    for t = 1:nT
        k = find(isPk(:, t));  [vals, o] = sort(Y(k, t), 'descend');
        P.tiles{t} = k(o);  P.values{t} = vals;
    end
end

% Author: Diellor Basha, 2026
