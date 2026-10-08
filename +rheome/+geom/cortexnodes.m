function C = cortexnodes(name, opts)
% GEOM.CORTEXNODES  One cortical position pyramid for both hemispheres, with whole-cortex node ids.
%
%   C = rheome.geom.cortexnodes('sub01')              % MaxDepth 7 per hemisphere
%   C = rheome.geom.cortexnodes(name, MaxDepth=6)
%   C.id(2, 5)                                         % whole-cortex id of left-tree heap node 5
%
% Every analysis builds its tiles per hemisphere (rheome.geom.tree on B.L.S / B.R.S), and a label written into
% the store needs ONE id per node across the whole cortex. The ids are the heap numbering of a whole-
% cortex tree: node 1 is the cortex, 2 the left hemisphere, 3 the right, and a node with heap id n at
% depth d of hemisphere h's own tree gets
%
%       id = n + (h - 1) * 2^d,     h = 2 (left), 3 (right)
%
% so parent = floor(id / 2) holds across the whole cortex and a node's hemisphere and depth are
% recoverable by arithmetic (depth = floor(log2(id)) - 1).
%
% ⚠ THE IDENTITY NEEDS A FULL TREE. Heap numbering holds only if every node splits down to MaxDepth
% (rheome.geom.tiles' .heap); this errors if either hemisphere's tree is not full, rather than issuing ids a
% later parent lookup would get wrong.
% ⚠ rheome.select.measure CANNOT CHECK cortical ids (the tree lives with the surface, not the recording), so
% this is the one function every writer and reader of cortical labels uses.
%
% OUTPUT (struct C)
%   .nodes   table, one row per node incl. the cortex root: node_id parent_id depth hemi n_vertices
%            is_leaf area diameter wavelength_mm centroid_x centroid_y centroid_z  (rheome.select.schema's
%            cortex_node columns, less recording_id)
%   .id      @(h, n) whole-cortex id of hemisphere h's heap node n (h = 2 left, 3 right)
%   .T       struct .L / .R: the per-hemisphere rheome.geom.tree tables
%   .S       struct .L / .R: the hemisphere surfaces
%
% See also: rheome.geom.tree, rheome.geom.tiles, rheome.select.measure, rheome.flow.labelalpha
%
% Author: Diellor Basha, 2026

    arguments
        name (1,1) string
        opts.MaxDepth (1,1) double {mustBeInteger, mustBePositive} = 7
    end
    B = rheome.load.bases(name);
    rows = {};  hs = ["L" "R"];
    C.T = struct();  C.S = struct();
    tot = 0;  cen = [0 0 0];  nvt = 0;
    for i = 1:2
        H = B.(char(hs(i)));  T = rheome.geom.tree(H.S, L=H.lbo.L, M=H.lbo.Mass, MaxDepth=opts.MaxDepth);
        nr = T.parent_id ~= 0;
        if ~all(T.parent_id(nr) == floor(T.node_id(nr) / 2)) || height(T) ~= 2^(opts.MaxDepth+1) - 1
            error('geom:cortexnodes:heap', ['The %s tree is not full to depth %d, so heap ids would ' ...
                'not identify parents. Lower MaxDepth.'], hs(i), opts.MaxDepth);
        end
        h = i + 1;  gid = T.node_id + (h - 1) .* 2.^T.depth;
        pid = floor(gid / 2);
        rows{end+1} = table(gid, pid, T.depth + 1, repmat(hs(i), height(T), 1), T.n_vertices, T.is_leaf, ...
            T.area, T.diameter, T.wavelength_mm, T.centroid(:,1), T.centroid(:,2), T.centroid(:,3), ...
            'VariableNames', {'node_id','parent_id','depth','hemi','n_vertices','is_leaf','area', ...
            'diameter','wavelength_mm','centroid_x','centroid_y','centroid_z'}); %#ok<AGROW>
        C.T.(char(hs(i))) = T;  C.S.(char(hs(i))) = H.S;
        tot = tot + T.area(1);  cen = cen + T.area(1) * T.centroid(1,:);  nvt = nvt + T.n_vertices(1);
    end
    root = table(1, 0, 0, "LR", nvt, false, tot, 2*sqrt(tot/pi), 1e3*2*sqrt(tot/pi), cen(1)/tot, ...
        cen(2)/tot, cen(3)/tot, 'VariableNames', rows{1}.Properties.VariableNames);
    C.nodes = sortrows(vertcat(root, rows{:}), 'node_id');
    C.id = @(h, n) n + (h - 1) .* 2.^floor(log2(n));
    if numel(unique(C.nodes.node_id)) ~= height(C.nodes)
        error('geom:cortexnodes:unique', 'Whole-cortex ids are not unique.');
    end
end

% Author: Diellor Basha, 2026
