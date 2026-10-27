function S = nxrsurface(store, subject, opts)
% IO.READ.NXRSURFACE  Read a subject's cortical surface and its medial wall from an nxr store (the editor's dataset folder).
%
%   S = rheome.io.read.nxrsurface(store, subject)
%   S = rheome.io.read.nxrsurface(store, subject, Surface="cortex_white")
%
% An nxr store is a folder holding dataset.sqlite (the catalogue) and one <subject>.nxr.zarr per subject. The surface is the
% subject's PRIMARY one (or the surface named by Surface=); its medial wall is the vertices its Desikan-Killiany atlas codes
% 'unassigned' (or 'unknown'), the wall rheome.geom.glue cuts. The lookups are nxr-convert's (eigen/hodge_dirac.py modes_input
% and medial_wall), so the same subject gives the same arrays.
%
% INPUTS:  store    the dataset folder (holds dataset.sqlite)     subject  the subject's name (its subject.name row)
% OPTIONS: Surface  a surface's name ("" = the primary)
%
% OUTPUT (struct S)
%   .Vertices [nV x 3] as stored (metres, float64)   .Faces [nF x 3] 1-based   .nV .nF
%   .Wall     [nV x 1] logical, the medial wall      .Comment the surface's name   .Subject   .SurfaceFile its positions' folder
%
% ⚠ The editor receives positions as float32: for a basis the editor will accept (rheome.geom.glue -> rheome.eigen.hodge_dirac,
%   as nxr-convert prepares it) round first, double(single(S.Vertices)).
%
% Reads the catalogue with the sqlite3 command-line tool (read-only; MATLAB's sqlite needs the Database Toolbox) and each
% array straight from its chunk file.
% shortcut: only uncompressed, single-chunk zarr v3 arrays (what nxr-convert writes for surfaces and atlases); any other
%   layout is refused, not misread. Upgrade to zarrread (R2024b+) when the toolbox moves past R2023b.
%
% See also: rheome.geom.glue, rheome.io.read.surface, rheome.eigen.hodge_dirac
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, opts.Surface (1,1) string = "", end
db = fullfile(store, "dataset.sqlite");
if ~isfile(db), error('io:read:nxrsurface:notFound', 'No nxr store at %s (dataset.sqlite missing)', store); end
q = @(s) "'" + replace(s, "'", "''") + "'";
pick = "m.is_primary = 1";
if opts.Surface ~= "", pick = "m.type = 'surface' AND m.name = " + q(opts.Surface); end
r = i_sql(db, "SELECT u.path AS store, m.name, g.positions_path AS positions, " + ...
    "(SELECT path FROM topology_cell WHERE topology_id = m.topology_id AND rank = 2) AS faces, " + ...
    "(SELECT s.array_path FROM selection s JOIN selection p ON p.id = s.member_of_id WHERE p.manifold_id = m.id AND p.name = m.name || '_atlases' AND s.name = 'Desikan-Killiany') AS atlas, " + ...
    "(SELECT e.code FROM selection s JOIN selection p ON p.id = s.member_of_id JOIN dictionary_entry e ON e.dictionary_id = s.dictionary_id " + ...
    " WHERE p.manifold_id = m.id AND p.name = m.name || '_atlases' AND s.name = 'Desikan-Killiany' AND e.name IN ('unassigned', 'unknown') ORDER BY e.name LIMIT 1) AS wall " + ...
    "FROM manifold m JOIN geometry g ON g.id = m.geometry_id JOIN subject u ON u.id = m.subject_id " + ...
    "WHERE u.name = " + q(subject) + " AND u.status = 'complete' AND " + pick);
if isempty(r)
    error('io:read:nxrsurface:noSurface', '%s: no %s in %s', subject, i_ifelse(opts.Surface == "", "primary surface", "surface '" + opts.Surface + "'"), store);
end
if isempty(r(1).atlas) || isempty(r(1).wall)
    error('io:read:nxrsurface:noWall', '%s: %s carries no Desikan-Killiany atlas with its medial wall - the hemispheres are glued by anatomy only', subject, r(1).name);
end
z = fullfile(store, r(1).store);
P = i_array(fullfile(z, r(1).positions));
F = i_array(fullfile(z, r(1).faces));
A = i_array(fullfile(z, r(1).atlas));
if size(P,2) ~= 3 || size(F,2) ~= 3 || numel(A) ~= size(P,1)
    error('io:read:nxrsurface:shape', '%s: positions %s, faces %s, atlas %d do not make one surface', subject, mat2str(size(P)), mat2str(size(F)), numel(A));
end
S = struct('Vertices', double(P), 'Faces', double(F) + 1, 'nV', size(P,1), 'nF', size(F,1), 'Wall', A(:) == r(1).wall, ...
    'Comment', r(1).name, 'Subject', char(subject), 'SurfaceFile', fullfile(z, r(1).positions));
end

function r = i_sql(db, sql)
[st, out] = system(sprintf('sqlite3 -readonly -json "%s" "%s"', db, replace(sql, '"', '\"')));
if st ~= 0, error('io:read:nxrsurface:sqlite', 'sqlite3 failed on %s: %s', db, strtrim(out)); end
if isempty(strtrim(out)), r = []; else, r = jsondecode(out); end
end

function x = i_array(folder)
% an uncompressed single-chunk zarr v3 array, C order, as an [n x m] (or [n x 1]) matrix of its own class
m = jsondecode(fileread(fullfile(folder, 'zarr.json')));
types = struct('float64', 'double', 'float32', 'single', 'int32', 'int32', 'uint32', 'uint32', 'int64', 'int64', 'uint8', 'uint8', 'int16', 'int16');
ok = isfield(m, 'zarr_format') && m.zarr_format == 3 && isfield(types, m.data_type) && numel(m.codecs) == 1 && ...
    strcmp(m.codecs(1).name, 'bytes') && isequal(m.chunk_grid.configuration.chunk_shape(:), m.shape(:)) && ...
    strcmp(m.chunk_key_encoding.name, 'default') && numel(m.shape) <= 2;
if ~ok, error('io:read:nxrsurface:layout', '%s: not an uncompressed single-chunk zarr v3 array (refused rather than misread)', folder); end
if isfield(m.codecs(1), 'configuration') && isfield(m.codecs(1).configuration, 'endian') && ~strcmp(m.codecs(1).configuration.endian, 'little')
    error('io:read:nxrsurface:layout', '%s: big-endian', folder);
end
key = strjoin(repmat("0", 1, numel(m.shape)), m.chunk_key_encoding.configuration.separator);
f = fullfile(folder, "c" + m.chunk_key_encoding.configuration.separator + key);
fid = fopen(f, 'r', 'ieee-le');
if fid < 0, error('io:read:nxrsurface:notFound', '%s: no chunk file (an all-fill array is not supported)', folder); end
c = onCleanup(@() fclose(fid));
t = types.(m.data_type);
x = fread(fid, inf, [t '=>' t]);
if numel(x) ~= prod(m.shape), error('io:read:nxrsurface:layout', '%s: %d values for shape %s', folder, numel(x), mat2str(m.shape)); end
if numel(m.shape) == 2, x = reshape(x, m.shape(2), m.shape(1))'; end
end

function v = i_ifelse(c, a, b)
if c, v = a; else, v = b; end
end
