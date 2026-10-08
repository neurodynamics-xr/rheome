function t = label(db, rows, opts)
% SELECT.LABEL  Append labels to tiles in one transaction: idempotent, atomic, rolled back by id.
%
%   t = rheome.select.label(db, rows, Kind="bad", Source="manual", Author="")
%       rows  table with level, k and optionally channel_id (0 = none), band_id (0 = none), value
%       t     the txn row (struct): txn_id, count, content_hash, created, and .existed
%
% One transaction: the txn row and the label rows are written together, atomically
% (select/private/sel_labelio). The content hash (MD5 of the sorted rows, kind and
% source) is unique per recording: submitting the same rows again returns the existing
% transaction without writing. In Postgres this is BEGIN; INSERT txn; INSERT label;
% COMMIT under the unique index txn_content (rheome.select.ddl).
%
% Labels are where bad segments, artefact marks and any later labelling go; features
% never move here (design §4).
%
% See also: rheome.select.labels, rheome.select.rollback
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        rows table
        opts.Kind   (1,1) string
        opts.Source (1,1) string = "manual"
        opts.Author (1,1) string = ""
    end
    rid = string(db.recording_id);
    n = height(rows);
    if ~all(ismember({'level','k'}, rows.Properties.VariableNames))
        error('select:label:rows', 'rows needs level and k.');
    end
    for c = {'channel_id','band_id','value'}
        if ~ismember(c{1}, rows.Properties.VariableNames)
            if strcmp(c{1}, 'value'), rows.(c{1}) = nan(n, 1); else, rows.(c{1}) = zeros(n, 1); end
        end
    end
    g = db.grid;
    if any(rows.level < 0 | rows.level > g.Lmax) || any(rows.k < 1 | rows.k > g.K(rows.level + 1)')
        error('select:label:tile', 'A row names a tile outside the grid.');
    end
    rows = sortrows(rows(:, {'channel_id','level','k','band_id','value'}), {'channel_id','level','k','band_id'});
    content = jsonencode(struct('kind', opts.Kind, 'source', opts.Source, 'rows', table2struct(rows)));
    h = sel_hash(content);
    [L, T] = sel_labelio(db);
    hit = find(T.recording_id == rid & T.content_hash == string(h), 1);
    if ~isempty(hit)
        t = table2struct(T(hit, :));  t.existed = true;
        return
    end
    id = max([0; T.txn_id]) + 1;
    lid0 = max([0; L.label_id]);
    new = table((lid0 + (1:n))', repmat(rid, n, 1), rows.channel_id, rows.level, rows.k, rows.band_id, ...
                repmat(opts.Kind, n, 1), rows.value, repmat(opts.Source, n, 1), repmat(id, n, 1), ...
                'VariableNames', L.Properties.VariableNames);
    tr = table(id, rid, opts.Kind, n, string(h), opts.Author, string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
               'VariableNames', T.Properties.VariableNames);
    [~, ~, M] = sel_labelio(db);
    sel_labelio(db, [L; new], [T; tr], M);           % measurements pass through untouched
    t = table2struct(tr);  t.existed = false;
end
% Author: Diellor Basha, 2026
