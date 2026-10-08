function [R, cost] = tree(db, q)
% SELECT.TREE  Descend the sensor tree from the root, pruning with the bound, to the query's depth.
%
%   [R, cost] = rheome.select.tree(db, q)      q = rheome.select.query(..., Scope="group", Depth=d)
%
% The position-axis twin of rheome.select.frames: at the query's time level, the survivors at
% tree depth d are the nodes whose parent survived and whose bound statistic passes the
% parent test; at Depth the predicate itself is applied. Group rows are exact sums and
% maxima over member channels, so the bound is exact and tree(db, q) == scan(db, q) for
% the same query (tests/tSelectGroups.m). In SQL it is the same semi-join chain over
% sensor_group.parent_id instead of the tile's parent key.
%
% The two descents compose: frames on the root node finds WHEN, tree at those frames
% finds WHERE.
%
% Author: Diellor Basha, 2026

    t0 = tic;  b0 = db.cost('bytes');
    if ~strcmp(q.scope, 'group') || isempty(q.depth)
        error('select:tree:query', 'rheome.select.tree needs Scope="group" and a Depth.');
    end
    tr = db.tree;  nodes = db.groupNodes;
    depths = 0:q.depth;
    plan = table('Size', [0 4], 'VariableTypes', {'double','double','double','double'}, ...
                 'VariableNames', {'depth','candidates','examined','survivors'});
    % the scope for theta: every internal node at the target depth (as scan sees it)
    [Vq, sel, bandsAt] = sel_values(db, q.level, q.stat, [], q.band, 'group', [], q.depth);
    theta = sel_threshold(q.threshold, Vq);
    rule = q.bound.rule;
    if strcmp(rule, 'none')
        hits = sel_pred(Vq, q.op, theta);
        plan(1, :) = {q.depth, numel(Vq), numel(Vq), nnz(hits)};
        R = sel_result(db, q, q.level, hits, Vq, sel, bandsAt);
        cost = struct('executor', 'tree', 'theta', theta, 'tilesExamined', numel(Vq), 'tilesReturned', height(R), 'bytes', db.cost('bytes') - b0, 'seconds', toc(t0), 'plan', plan);
        return
    end
    K = db.grid.K(q.level+1);  nB = max(1, numel(bandsAt));
    survivors = containers.Map('KeyType', 'double', 'ValueType', 'any');       % node -> [K x nB] logical
    examined = 0;
    for d = depths
        atDepth = tr.node_id(tr.depth == d & ~tr.is_leaf)';
        atDepth = atDepth(ismember(atDepth, nodes));
        if d == 0, cand = atDepth;
        else, cand = atDepth(arrayfun(@(n) isKey(survivors, tr.parent_id(n)), atDepth)); end
        if isempty(cand), plan(end+1, :) = {d, 0, 0, 0}; break; end %#ok<AGROW>
        if d == q.depth
            [V, ~, ~] = sel_values(db, q.level, q.stat, [], q.band, 'group', cand, []);
            pass = sel_pred(V, q.op, theta);
        else
            B = sel_values(db, q.level, q.bound.stat, [], bandsAt, 'group', cand, []);
            switch rule
                case 'ge',      pass = B >= theta;
                case 'le',      pass = B <= theta;
                case 'ge_sqrt', pass = B.^2 >= theta;
            end
        end
        nSurv = 0;
        for i = 1:numel(cand)
            m = reshape(pass(:, i, :), K, nB);
            if d > 0, m = m & reshape(survivors(tr.parent_id(cand(i))), K, nB); end   % the parent's tiles
            if any(m, 'all'), survivors(cand(i)) = m;  nSurv = nSurv + 1; end
        end
        examined = examined + numel(cand) * K * nB;
        plan(end+1, :) = {d, numel(cand), numel(cand) * K * nB, nSurv}; %#ok<AGROW>
        if nSurv == 0, break; end
    end
    % assemble hits at the target depth in the scope's column order
    hits = false(K, numel(sel), nB);
    for i = 1:numel(sel)
        if isKey(survivors, sel(i)) && tr.depth(sel(i)) == q.depth
            hits(:, i, :) = reshape(survivors(sel(i)), K, 1, nB);
        end
    end
    R = sel_result(db, q, q.level, hits, Vq, sel, bandsAt);
    cost = struct('executor', 'tree', 'theta', theta, 'tilesExamined', examined, 'tilesReturned', height(R), ...
                  'bytes', db.cost('bytes') - b0, 'seconds', toc(t0), 'plan', plan);
end
% Author: Diellor Basha, 2026
