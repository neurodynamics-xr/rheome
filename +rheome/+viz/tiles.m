function L = tiles(T, level, opts)
% VIZ.TILES  A tile level of a cortical pyramid as a viewer labelling: tiles by level on the cortex.
%
%   L = rheome.viz.tiles(T, level)                                  % T = rheome.geom.tree(S, ...)
%   L = rheome.viz.tiles({TL, TR}, 8, Vertices={iL, iR}, NV=nV)     % per-hemisphere trees (the atlas's ladder)
%   rheome.view(S, x, Labellings=L, Labelling=L.name)
%
% Each vertex gets the code of the tile that holds it at depth LEVEL (a leaf above that depth keeps its own tile), so
% the bundle carries the partition and the viewer paints it (setting 'labelling'). The tiles are numbered in table order,
% hemisphere after hemisphere; entries name each by its node and depth.
%
% INPUTS:
%   T      a rheome.geom.tree table, or a cell of them (one a hemisphere)
%   level  the depth shown (0: the root)
% OPTIONS:
%   Vertices  with a cell T: a cell of the vertex indices (1-based, into the whole surface) each tree's members index
%   NV        with a cell T: the whole surface's vertex count
%   Name      the labelling's name (default "tiles level <level>")
%
% OUTPUT:
%   L : struct .name .codes [nV x 1] int (-1: no tile) .entries (.code .name), for rheome.io.write.bundle 'Labellings'
%
% NOT YET (the viewer): tile OUTLINES with major/minor weights (cortex.marks 'tile outlines' and 'coarser-level
% boundaries') -- the viewer paints a labelling's fill only; filed with Product.
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', 0.07*V, 'Faces', F);
%   L = rheome.viz.tiles(rheome.geom.tree(S, MaxDepth=3), 2);
%   assert(numel(L.entries) == 4 && all(L.codes >= 0))
%
% See also: rheome.geom.tree, rheome.io.write.bundle, rheome.view
%
% Author: Diellor Basha, 2026
    arguments
        T
        level (1,1) double {mustBeInteger, mustBeNonnegative}
        opts.Vertices = {}
        opts.NV (1,1) double = 0
        opts.Name = ""
    end
    if istable(T), T = {T}; rows = {[]};
    else
        rows = opts.Vertices;
        if ~iscell(T) || ~iscell(rows) || numel(rows) ~= numel(T) || opts.NV < 1
            error('rheome:viz:tiles', 'A cell of trees needs Vertices= (one index vector each) and NV=.');
        end
    end
    nV = opts.NV;  if nV < 1, nV = numel(T{1}.members{1}); end
    codes = -ones(nV, 1, 'int32');  entries = struct('code', {}, 'name', {});
    for h = 1:numel(T)
        t = T{h};
        shown = t.depth == level | (t.is_leaf & t.depth < level);
        for k = find(shown)'
            m = t.members{k};
            if ~isempty(rows{h}), m = rows{h}(m); end
            c = numel(entries);
            codes(m) = c;
            entries(end+1) = struct('code', c, 'name', sprintf('node %d (depth %d)', t.node_id(k), t.depth(k))); %#ok<AGROW>
        end
    end
    name = string(opts.Name);  if name == "", name = "tiles level " + level; end
    L = struct('name', char(name), 'codes', codes, 'entries', entries);
end
