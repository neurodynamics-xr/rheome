function Y = tilemean(G, X, a, stat)
% GEOM.TILEMEAN  Summarise a vertex map over each tile: the observation a tile makes of a field.
%
%   Y = rheome.geom.tilemean(G, X, a)              % area-weighted mean, [nTile x nT]
%   Y = rheome.geom.tilemean(G, X, a, 'median')    % per-tile median (unweighted)
%
% A tile is a set of vertices, so a tile observes a field by pooling it. The mean is one sparse
% product, Y = (P' * (a .* X)) ./ (P' * a), with P the membership from rheome.geom.tiles and a the vertex
% areas (row sums of the mass matrix) -- area, not vertex count, because the mesh is not uniform.
% Because it is linear it rolls up the tree exactly: a parent's mean is the area-weighted mean of its
% children's. The median does neither, and is offered for robustness against a few extreme vertices.
%
% ⚠⚠ A BAND-PASS MAP AVERAGES TOWARD ZERO OVER A TILE MUCH LARGER THAN ITS SCALE. A mexhat atom -- and
% the curl of a vortex -- has a positive core and a negative ring that integrate to nearly zero, so a
% tile that swallows both reports almost nothing, and one that straddles the ring can report the
% opposite sign. The tile has to be matched to the spatial band: that is the separable version of the
% joint tiling's diagonal. Planting measures where it breaks.
%
% INPUTS
%   G     rheome.geom.tiles struct (.P, .tileOf)      X  [nV x nT] vertex maps
%   a     [nV x 1] vertex areas, or [] for an unweighted mean
%   stat  'mean' (default) | 'median'
%
% OUTPUT
%   Y     [nTile x nT]
%
% See also: rheome.geom.tiles, rheome.detect.tilepeaks, rheome.ingest.groups
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(a), a = ones(size(X, 1), 1); end
    if nargin < 4, stat = 'mean'; end
    a = double(a(:));
    switch lower(stat)
        case 'mean'
            P = double(G.P);
            Y = (P' * (a .* X)) ./ (P' * a);
        case 'median'
            nTile = size(G.P, 2);  Y = zeros(nTile, size(X, 2));
            for k = 1:nTile, Y(k, :) = median(X(G.tileOf == k, :), 1); end
        otherwise
            error('geom:tilemean:stat', 'stat must be ''mean'' or ''median'', got ''%s''.', stat);
    end
end

% Author: Diellor Basha, 2026
