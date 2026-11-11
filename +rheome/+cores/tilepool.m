function pool = tilepool(kit, V, lab, support)
% CORES.TILEPOOL  The pooling of one ladder level's tiles: every support vertex's vectors are carried to its tile's centre (nxr Tiles).
%
%   pool = rheome.cores.tilepool(kit, V, lab, support)
%
% kit: rheome.geom.GeometryKit of the surface; V [N x 3] mm; lab [N x 1] tile codes of ONE level (a column of rheome.cores.ladder);
% support: 1-based vertices cores may snap to. A tile's centre is the support vertex nearest the centroid of its support vertices
% (rheome.cores.tilecentres); pool.Centre [ntile x 1] the centres (0: empty). Build it once a level and subject (one transport solve a tile).
%
% See also: rheome.cores.crossings, rheome.cores.Pooling
%
% Author: Diellor Basha, 2026
    t = -ones(size(V, 1), 1);  t(support) = double(lab(support));
    c = rheome.cores.tilecentres(V, t);
    a = zeros(size(V, 1), 1);  on = t >= 0;  a(on) = c(t(on) + 1);
    pool = rheome.cores.Pooling(kit, a);  pool.Centre = c;
end
