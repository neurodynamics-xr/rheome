function r = nxrsql(store, sql)
% NXRSQL  Rows of an nxr store's catalogue (dataset.sqlite), read-only, as a struct array ([] when none).
% In-process through MATLAB's Python (its standard sqlite3 module, a read-only URI): MATLAB's own sqlite needs the Database Toolbox.
% The catalogue must be schema 50 (PRAGMA user_version, nxr.subject@6.0), as nxr reads it: a 49 is refused with the command that migrates it,
% any other version with 're-convert it from its source'.
db = fullfile(store, "dataset.sqlite");
if ~isfile(db), error('io:read:nxr:notFound', 'No nxr store at %s (dataset.sqlite missing)', store); end
code = ["import sqlite3, json, urllib.parse"
        "cx = sqlite3.connect('file:' + urllib.parse.quote(db) + '?mode=ro', uri=True)"
        "v = cx.execute('PRAGMA user_version').fetchone()[0]"
        "cur = cx.execute(sql) if v == 50 else None"
        "cols = [d[0] for d in cur.description] if cur else []"
        "out = json.dumps({'version': v, 'rows': [dict(zip(cols, row)) for row in cur.fetchall()] if cur else []})"
        "cx.close()"];
try
    out = pyrun(code, "out", db=char(db), sql=char(sql));
catch e
    error('io:read:nxr:sqlite', 'reading %s through MATLAB''s Python (pyenv) failed: %s', db, e.message);
end
r = jsondecode(char(out));
if r.version ~= 50
    how = "re-convert it from its source";
    if r.version == 49, how = "run `nxr migrate --dataset <dataset folder>` (it keeps the old database and the arrays; add `--dry-run` to see what it will do first)"; end
    error('io:read:nxr:schema', '%s: schema version %d; Rheome reads version 50 only (as nxr does) - %s', db, r.version, how);
end
r = r.rows;
if isempty(r), r = []; end
end
