function t = measure(db, rows, opts)
% SELECT.MEASURE  Write measurements onto pyramid nodes, in one transaction.
%
%   t = rheome.select.measure(db, rows, Kind="vortex_count")
%   t = rheome.select.measure(db, rows, Kind="rotation_sense", Scope="group", Source="flow")
%       rows  table with level, k and value, optionally unit_id (0 = none) and band_id (0 = none)
%       t     the txn row: txn_id, count, content_hash, created, and .existed
%
% ⭐ THE SECOND BOOKKEEPING MATRIX. The feature relations are dense and mergeable: one row
% per (unit, tile), a parent exactly the merge of its children, computed once and rolled up.
% This one is sparse and NOT mergeable: one row per (node, kind), holding whatever was
% measured where it was measured -- a vortex count, a rotation sense, a fitted slope, a
% track id. The pyramid is then just a graph to hang things on, which is all it needs to be
% for labelling, and the two matrices meet again in rheome.select.context.
%
% ⚠ A VALUE HERE SAYS NOTHING ABOUT ITS PARENT. Nothing rolls up, no bound holds, so a query
% on a measurement scans its level or uses an index; it cannot prune by descending the way
% rheome.select.frames does on the moments. Two levels carrying the same kind are two independent
% observations of overlapping spans -- never read one as a summary of the other.
%
% ⚠ ADDRESSED BY SCOPE AND UNIT, unlike `label`, which only knows channels. A measurement on
% a sensor-tree node (Scope="group", unit_id = node id) is the common case for flow, where
% the thing measured is a pattern over a patch, not a number on one sensor.
%
% ⭐ Scope="cortex" addresses a rheome.geom.tree node on the far side of the inverse (the cortex_node
% dimension), which is where a time-averaged magnitude, divergence or curl under a cortical tile
% belongs. Those are exactly the quantities that do NOT merge -- a curl on a patch is not the sum
% of the curls of its halves -- so this relation, not a feature relation, is their home.
% ⚠ A cortical measurement is only as good as the instrument under that node, measured by
% planting. Location recovers to 43-52 mm and SIZE DOES NOT RECOVER AT ALL below ~20 dB,
% so write a magnitude or a curl here and do not write a "blob size".
%
% ⚠ rheome.select.labelkinds is consulted on every write: a kind registered "unmeasurable" (e.g. a blob or
% envelope size) is REFUSED with its evidence in the message.
%
% Same transaction machinery as rheome.select.label: atomic write, idempotent by content hash,
% removed by rheome.select.rollback. Labels and measurements share the sidecar file and the
% transaction log, so one rollback cannot leave half a write behind.
%
% See also: rheome.select.measures, rheome.select.context, rheome.select.label, rheome.select.rollback
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        rows table
        opts.Kind   (1,1) string
        opts.Scope  (1,1) string {mustBeMember(opts.Scope, ["channel","group","cortex","none"])} = "none"
        opts.Source (1,1) string = "manual"
        opts.Author (1,1) string = ""
    end
    rid = string(db.recording_id);
    n = height(rows);
    need = {'level','k','value'};
    if ~all(ismember(need, rows.Properties.VariableNames))
        error('select:measure:rows', 'rows needs level, k and value.');
    end
    for c = {'unit_id','band_id'}
        if ~ismember(c{1}, rows.Properties.VariableNames), rows.(c{1}) = zeros(n, 1); end
    end
    g = db.grid;
    if any(rows.level < 0 | rows.level > g.Lmax) || any(rows.k < 1 | rows.k > g.K(rows.level + 1)')
        error('select:measure:tile', 'A row names a tile outside the grid.');
    end
    i_checkunits(db, opts.Scope, rows.unit_id);
    % ⚠ A KIND THE CALIBRATION REGISTRY HAS SHOWN TO BE MEANINGLESS IS REFUSED, not written. An
    % envelope or blob size reads the instrument's rendering (~150 mm) whatever the source's size, so a
    % row of it would be a confident wrong number in a column named for the brain. Unregistered kinds
    % are accepted as before; the registry is where a kind earns its calibration (rheome.select.labelkinds).
    reg = rheome.select.labelkinds(opts.Kind);
    if ~isempty(reg) && reg.calibration == "unmeasurable"
        error('select:measure:unmeasurable', ['Kind "%s" is registered as UNMEASURABLE through this ' ...
            'instrument (%s). It is not written.'], opts.Kind, reg.evidence);
    end

    rows = sortrows(rows(:, {'unit_id','level','k','band_id','value'}), {'unit_id','level','k','band_id'});
    content = jsonencode(struct('kind', opts.Kind, 'scope', opts.Scope, 'source', opts.Source, ...
                                'rows', table2struct(rows)));
    h = sel_hash(content);
    [L, T, M] = sel_labelio(db);
    hit = find(T.recording_id == rid & T.content_hash == string(h), 1);
    if ~isempty(hit)
        t = table2struct(T(hit, :));  t.existed = true;
        return
    end
    id = max([0; T.txn_id]) + 1;
    mid0 = max([0; M.measure_id]);
    new = table((mid0 + (1:n))', repmat(rid, n, 1), repmat(opts.Scope, n, 1), rows.unit_id, ...
                rows.level, rows.k, rows.band_id, repmat(opts.Kind, n, 1), rows.value, ...
                repmat(opts.Source, n, 1), repmat(id, n, 1), ...
                'VariableNames', M.Properties.VariableNames);
    tr = table(id, rid, opts.Kind, n, string(h), opts.Author, ...
               string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
               'VariableNames', T.Properties.VariableNames);
    sel_labelio(db, L, [T; tr], [M; new]);
    t = table2struct(tr);  t.existed = false;
end

% A unit id that names nothing is caught here rather than at the join, where it would look
% like a missing context row instead of a bad write.
function i_checkunits(db, scope, u)
    u = u(u ~= 0);
    if isempty(u), return; end
    switch scope
        case "channel"
            if any(u < 1 | u > db.meta.C)
                error('select:measure:unit', 'A channel id is outside 1..%d.', db.meta.C);
            end
        case "group"
            if isempty(db.groupNodes)
                error('select:measure:unit', 'This store has no sensor groups; Scope must be "channel" or "none".');
            end
            miss = setdiff(u, db.groupNodes);
            if ~isempty(miss)
                error('select:measure:unit', 'Node %d is not an internal node with rows.', miss(1));
            end
        case "cortex"
            % ⚠ MEMBERSHIP IS NOT CHECKABLE HERE, and saying so is better than pretending. The
            % cortex_node dimension is declared in rheome.select.schema and enforced by the Postgres FK in
            % rheome.select.ddl, but a MAT store holds no cortex_node array -- the tree lives with the
            % subject's surface (rheome.geom.tree over rheome.load.bases), not with the recording. So this checks
            % only that the ids are positive integers, and a wrong node id will be accepted.
            % ⭐ The guard that matters is therefore on the WRITER: build the ids from the same
            % rheome.geom.tree call that produced the values (rheome.flow.sweep does).
            if any(u ~= fix(u)) || any(u < 1)
                error('select:measure:unit', 'A cortex node id must be a positive integer.');
            end
        otherwise
            error('select:measure:unit', 'Scope "none" takes unit_id 0.');
    end
end
% Author: Diellor Basha, 2026
