function [R, cost] = frames(db, q)
% SELECT.FRAMES  The pruned executor: descend the pyramid, discard subtrees the bound rules out.
%
%   [R, cost] = rheome.select.frames(db, q)
%
% From the top level down to q.level, the survivors at a level are the tiles whose parent
% survived and whose bound statistic passes the parent test (design §3.2); at q.level the
% predicate itself is applied to the survivors. Every stored statistic merges, so the
% bound is exact and no tile that satisfies the predicate is ever discarded:
% frames(db, q) == scan(db, q) (tests/tSelectFrames.m). Relative thresholds are resolved
% as in scan, over the scope at q.level, so both executors see one theta.
%
% cost.plan: one row per level visited -- candidates (children of survivors), examined,
% survivors. cost.tilesExamined counts the rows a database would touch through the
% semi-join chain; cost.bytes counts what the MATLAB store actually loaded.
%
% Author: Diellor Basha, 2026

    t0 = tic;  b0 = db.cost('bytes');
    g = db.grid;  Lq = q.level;  Ltop = g.Lmax;
    [Vq, chans, bandsAt] = sel_values(db, Lq, q.stat, q.channels, q.band, q.scope, q.groups, q.depth);   % the scope, and theta
    theta = sel_threshold(q.threshold, Vq);
    nC = numel(chans);  nB = max(1, numel(bandsAt));
    % bound values by level (band stats need the band's slot at each level)
    rule = q.bound.rule;
    plan = table('Size', [0 4], 'VariableTypes', {'double','double','double','double'}, ...
                 'VariableNames', {'level','candidates','examined','survivors'});
    if strcmp(rule, 'none')
        % nothing prunes: evaluate at the level directly (an index scan in the database)
        hits = sel_pred(Vq, q.op, theta);
        plan(1, :) = {Lq, numel(Vq), numel(Vq), nnz(hits)};
        R = sel_result(db, q, Lq, hits, Vq, chans, bandsAt);
        cost = i_cost(db, theta, plan, R, b0, t0, numel(Vq));
        return
    end
    examinedTotal = 0;
    surv = true(g.K(Ltop+1), nC, nB);                            % everything survives above the top
    for L = Ltop:-1:Lq
        K = g.K(L+1);
        if L == Ltop
            cand = true(K, nC, nB);
        else
            Kp = g.K(L+2);
            parentOf = min(floor(((1:K)' - 1) / 2) + 1, Kp);
            cand = surv(parentOf, :, :);
        end
        if L == Lq
            V = Vq;  pass = sel_pred(V, q.op, theta);
        else
            B = sel_values(db, L, q.bound.stat, chans, bandsAt, q.scope, chans, []);
            switch rule
                case 'ge',      pass = B >= theta;
                case 'le',      pass = B <= theta;
                case 'ge_sqrt', pass = B.^2 >= theta;             % envMax^2 bounds energy/n
            end
        end
        nCand = nnz(cand);
        surv = cand & pass;
        examinedTotal = examinedTotal + nCand;
        plan(end+1, :) = {L, nCand, nCand, nnz(surv)}; %#ok<AGROW>
        if ~any(surv, 'all')
            break
        end
    end
    if L > Lq, hits = false(g.K(Lq+1), nC, nB); else, hits = surv; end
    R = sel_result(db, q, Lq, hits, Vq, chans, bandsAt);
    cost = i_cost(db, theta, plan, R, b0, t0, examinedTotal);
end

function cost = i_cost(db, theta, plan, R, b0, t0, examined)
    cost = struct('executor', 'frames', 'theta', theta, 'tilesExamined', examined, ...
                  'tilesReturned', height(R), 'bytes', db.cost('bytes') - b0, 'seconds', toc(t0), 'plan', plan);
end
% Author: Diellor Basha, 2026
