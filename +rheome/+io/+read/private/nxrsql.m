function r = nxrsql(store, sql)
% NXRSQL  Rows of an nxr store's catalogue (dataset.sqlite), read-only, as a struct array ([] when none).
% Through the sqlite3 command-line tool: MATLAB's sqlite needs the Database Toolbox.
db = fullfile(store, "dataset.sqlite");
if ~isfile(db), error('io:read:nxr:notFound', 'No nxr store at %s (dataset.sqlite missing)', store); end
[st, out] = system(sprintf('sqlite3 -readonly -json "%s" "%s"', db, replace(sql, '"', '\"')));
if st ~= 0, error('io:read:nxr:sqlite', 'sqlite3 failed on %s: %s', db, strtrim(out)); end
if isempty(strtrim(out)), r = []; else, r = jsondecode(out); end
end
