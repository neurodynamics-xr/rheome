function c = tilecentres(P, tiles)
% CORES.TILECENTRES  The member vertex nearest each tile's centroid (nxr measure.vectors.tile_centres).
%
%   c = rheome.cores.tilecentres(P, tiles)
%
% P [N x 3] positions, tiles [N x 1] a code per vertex (0-based, negative off the tiles). c [max(tiles)+1 x 1]: the 1-based vertex of each tile
% that is nearest its members' centroid (ties: the first), 0 for a tile with no member.
%
% See also: rheome.cores.Pooling, rheome.cores.tilepool
%
% Author: Diellor Basha, 2026
    tiles = double(tiles(:));
    if isempty(tiles) || max(tiles) < 0, c = zeros(0, 1); return, end
    c = zeros(max(tiles) + 1, 1);
    for k = unique(tiles(tiles >= 0))'
        m = find(tiles == k);
        [~, j] = min(sum((P(m, :) - mean(P(m, :), 1)).^2, 2));
        c(k + 1) = m(j);
    end
end
