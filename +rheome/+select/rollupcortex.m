function R = rollupcortex(X, leafIds, K)
% SELECT.ROLLUPCORTEX  Sum finest-cell values up both trees: time by the store's levels, cortex by heap.
%
%   R = rheome.select.rollupcortex(X, leafIds, K)
%       X        [nLeaf x K(1) x p]  per-leaf values on the level-0 tiles (p stats or bands), SUMS
%       leafIds  [nLeaf x 1] whole-cortex heap ids of the leaves, all at one depth (rheome.geom.cortexnodes)
%       K        [1 x nLevel] tiles per time level, from db.grid.K (level 0 first)
%       R.ids    [nNode x 1] every node from the leaves up to the cortex root (1), ascending
%       R.lev    {1 x nLevel} [nNode x K(L+1) x p]
%
% The mergeable half of the store is SUMS, so a coarser cell is the sum of its finer ones on both axes:
%   time    level L tile k holds level-0 tiles k0 with floor((k0-1)/2^L)+1 = k -- the dyadic grid, with
%           the last tile of a level partial when its count is odd (75 -> 38 at level 6)
%   cortex  a node is the sum of its two children, parent = floor(id/2)
% Nothing is recomputed; every coarser cell is exact up to floating-point summation order.
%
% ⚠ ONLY SUMS MAY GO IN. An occupancy, a mean or a ratio does not roll up by addition -- store its
% numerator and denominator (on-area x time, area x time) and divide at read time.
%
% See also: rheome.flow.cortexfeatures, rheome.geom.cortexnodes, rheome.select.rows
%
% Author: Diellor Basha, 2026

    arguments
        X double
        leafIds (:,1) double {mustBeInteger, mustBePositive}
        K (1,:) double {mustBeInteger, mustBePositive}
    end
    [nLeaf, n0, p] = size(X);
    if nLeaf ~= numel(leafIds), error('select:rollupcortex:size', 'X has %d rows for %d leaf ids.', nLeaf, numel(leafIds)); end
    if n0 ~= K(1), error('select:rollupcortex:size', 'X has %d level-0 tiles; the grid has %d.', n0, K(1)); end
    d = unique(floor(log2(leafIds)));
    if numel(d) ~= 1, error('select:rollupcortex:depth', 'Leaves must share one depth.'); end
    % every ancestor, and a sparse leaf -> node membership (a node sums the leaves beneath it)
    ids = leafIds;  anc = leafIds;
    for s = 1:d
        anc = floor(anc / 2);  ids = [ids; anc]; %#ok<AGROW>
    end
    ids = unique(ids);
    rowsI = [];  colsI = [];  a = leafIds;
    for s = 0:d
        [~, loc] = ismember(a, ids);  rowsI = [rowsI; loc]; colsI = [colsI; (1:nLeaf)']; %#ok<AGROW>
        a = floor(a / 2);
    end
    Mn = sparse(rowsI, colsI, 1, numel(ids), nLeaf);
    R.ids = ids;  R.lev = cell(1, numel(K));
    for L = 0:numel(K) - 1
        kk = min(floor((0:n0-1)' / 2^L) + 1, K(L+1));
        St = sparse((1:n0)', kk, 1, n0, K(L+1));
        Y = zeros(numel(ids), K(L+1), p);
        for j = 1:p, Y(:,:,j) = Mn * X(:,:,j) * St; end
        R.lev{L+1} = Y;
    end
end

% Author: Diellor Basha, 2026
