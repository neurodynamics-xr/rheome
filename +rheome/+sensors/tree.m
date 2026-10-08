function T = tree(G, opts)
% SENSORS.TREE  Recursive spectral bisection of a sensor graph: the position pyramid.
%
%   T = rheome.sensors.tree(G)
%   T = rheome.sensors.tree(G, MinSize=1)
%
% Splits the array by the sign of the Fiedler vector of each part's own Laplacian, down to
% single sensors: a binary tree whose root is the whole array and whose leaves are the
% channels. On the CTF helmet the first cut follows left/right (|r| = 0.91 with the
% lateral coordinate), the second front/back, the third dorsal/lateral -- hemispheres,
% then quadrants, then lobe-sized patches of 18-49 sensors -- without being told any of
% it (docs/2026-09-22-ingest-notes.md).
%
% ⭐ A GROUP'S STATISTICS ARE SUMS AND MAXIMA OVER ITS SENSORS' ROWS, so the tree is a
% roll-up along the channel axis exactly as the dyadic grid is along time: a parent
% equals its children merged, and a query prunes down the tree from the root. No
% transform is involved; this axis is where "which part of the helmet" lives.
%
% ⚠ A disconnected part (Fiedler value ~ 0 with all members on one side) is split along
% its longest coordinate instead, so every internal node has two non-empty children.
%
% INPUTS:
%   G        rheome.sensors.graph struct (.L, .Vertices)
%   MinSize  stop splitting below this many sensors (1)
% OUTPUT (table T, one row per node, root first, children after parents, ids in that order):
%   node_id, parent_id (0 for the root), depth (root 0), n_sensors, is_leaf,
%   channel_id (the sensor for a leaf, 0 otherwise), members {row indices into the array},
%   centroid [x y z] m, diameter m, axis (the coordinate the cut followed, '' for leaves)
%
% See also: rheome.sensors.graph, rheome.ingest.groups, rheome.select.tree
%
% Author: Diellor Basha, 2026

    arguments
        G (1,1) struct
        opts.MinSize (1,1) double {mustBeInteger, mustBePositive} = 1
    end
    L = G.L;  P = G.Vertices;  n = size(P, 1);
    rows = {};
    queue = {struct('idx', (1:n)', 'parent', 0, 'depth', 0)};
    nextId = 1;
    while ~isempty(queue)
        nd = queue{1};  queue(1) = [];
        idx = nd.idx;  id = nextId;  nextId = nextId + 1;
        isLeaf = numel(idx) <= opts.MinSize;
        ax = '';
        if ~isLeaf
            [a, b, ax] = i_split(L, idx, P);
            queue{end+1} = struct('idx', a, 'parent', id, 'depth', nd.depth + 1); %#ok<AGROW>
            queue{end+1} = struct('idx', b, 'parent', id, 'depth', nd.depth + 1); %#ok<AGROW>
        end
        Q = P(idx, :);
        if numel(idx) > 1, [~, diam] = sen_metrics(Q); else, diam = 0; end      % no pdist: Statistics Toolbox
        ch = 0;  if numel(idx) == 1, ch = idx; end
        rows(end+1, :) = {id, nd.parent, nd.depth, numel(idx), isLeaf, ch, idx(:)', mean(Q, 1), diam, ax}; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'node_id','parent_id','depth','n_sensors','is_leaf','channel_id','members','centroid','diameter','axis'});
end

function [a, b, ax] = i_split(L, idx, P)
    Ls = full(L(idx, idx));  Ls = (Ls + Ls') / 2;
    [V, D] = eig(Ls);  [d, o] = sort(diag(D));  f = V(:, o(2));
    a = idx(f >= 0);  b = idx(f < 0);
    axn = {'x','y','z'};
    if d(2) > 1e-9 * max(d(end), 1) && ~isempty(a) && ~isempty(b)
        c = zeros(1, 3);  for k = 1:3, c(k) = abs(corr(f, P(idx, k))); end
        [~, k] = max(c);  ax = axn{k};
        return
    end
    % disconnected or degenerate: split along the longest coordinate
    ext = max(P(idx, :), [], 1) - min(P(idx, :), [], 1);  [~, k] = max(ext);  ax = axn{k};
    m = median(P(idx, k));
    a = idx(P(idx, k) >= m);  b = idx(P(idx, k) < m);
    if isempty(b), a = idx(1:floor(end/2));  b = idx(floor(end/2)+1:end); end
end
% Author: Diellor Basha, 2026
