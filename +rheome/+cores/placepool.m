function pool = placepool(kit, V, support, places)
% CORES.PLACEPOOL  The pooling of the direction consistency: every support vertex is pooled at the nearest of the places (nxr place_pooling).
%
%   pool = rheome.cores.placepool(kit, V, support, places)
%
% kit: rheome.geom.GeometryKit; V [N x 3] mm; support: 1-based vertices; places: indices INTO the support (1-based). Each support vertex goes to the
% place nearest it (3-D); pool.PlaceOf [N x 1] is that place's position in `places` (0 off the support).
%
% See also: rheome.cores.placeconsistency, rheome.cores.Pooling
%
% Author: Diellor Basha, 2026
    Vs = V(support, :);
    of = knnsearch(Vs(places, :), Vs);
    a = zeros(size(V, 1), 1);  a(support) = support(places(of));
    pool = rheome.cores.Pooling(kit, a);
    pool.PlaceOf = zeros(size(V, 1), 1);  pool.PlaceOf(support) = of;
end
