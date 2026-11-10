function r = nxrsql(store, sql)
% NXRSQL  Rows of an nxr store's catalogue (dataset.sqlite), read-only, as a struct array ([] when none).
% In-process through MATLAB's Python (its standard sqlite3 module, a read-only URI): MATLAB's own sqlite needs the Database Toolbox.
db = fullfile(store, "dataset.sqlite");
if ~isfile(db), error('io:read:nxr:notFound', 'No nxr store at %s (dataset.sqlite missing)', store); end
code = ["import sqlite3, json, urllib.parse"
        "cx = sqlite3.connect('file:' + urllib.parse.quote(db) + '?mode=ro', uri=True)"
        "cur = cx.execute(sql)"
        "cols = [d[0] for d in cur.description]"
        "out = json.dumps([dict(zip(cols, row)) for row in cur.fetchall()])"
        "cx.close()"];
try
    out = pyrun(code, "out", db=char(db), sql=char(sql));
catch e
    error('io:read:nxr:sqlite', 'reading %s through MATLAB''s Python (pyenv) failed: %s', db, e.message);
end
r = jsondecode(char(out));
if isempty(r), r = []; end
end
