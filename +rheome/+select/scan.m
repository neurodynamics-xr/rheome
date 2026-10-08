function [R, cost] = scan(db, q)
% SELECT.SCAN  The reference executor: a full scan of the relation at the query's level.
%
%   [R, cost] = rheome.select.scan(db, q)
%
% Materialises the statistic for every tile in scope, resolves the threshold, applies the
% predicate and the duration rule. Slow by design; rheome.select.frames must equal it.
%
% Author: Diellor Basha, 2026

    t0 = tic;  before = db.cost('tiles');  b0 = db.cost('bytes');
    [V, chans, bandsAt] = sel_values(db, q.level, q.stat, q.channels, q.band, q.scope, q.groups, q.depth);
    theta = sel_threshold(q.threshold, V);
    hits = sel_pred(V, q.op, theta);
    R = sel_result(db, q, q.level, hits, V, chans, bandsAt);
    cost = struct('executor', 'scan', 'theta', theta, 'tilesExamined', db.cost('tiles') - before, ...
                  'tilesReturned', height(R), 'bytes', db.cost('bytes') - b0, 'seconds', toc(t0), 'plan', table());
end
% Author: Diellor Basha, 2026
