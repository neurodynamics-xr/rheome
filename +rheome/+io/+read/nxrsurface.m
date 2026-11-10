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
%   .Sphere   [nV x 3] its registered sphere (the surface <name>_sphere, as stored), [] when none
%
% ⚠ The editor receives positions as float32: for a basis the editor will accept (rheome.geom.glue -> rheome.eigen.hodge_dirac,
%   as nxr-convert prepares it) round first, double(single(S.Vertices)).
%
% Reads the catalogue with the sqlite3 command-line tool (read-only; MATLAB's sqlite needs the Database Toolbox) and the
% arrays with rheome.io.read.zarr (uncompressed zarr v3; anything else is refused, not misread).
%
% See also: rheome.io.read.zarr, rheome.geom.glue, rheome.io.read.surface, rheome.eigen.hodge_dirac
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, opts.Surface (1,1) string = "", end
q = @(s) "'" + replace(s, "'", "''") + "'";
pick = "m.is_primary = 1";
if opts.Surface ~= "", pick = "m.type = 'surface' AND m.name = " + q(opts.Surface); end
r = nxrsql(store, "SELECT u.path AS store, m.name, g.positions_path AS positions, " + ...
    "(SELECT path FROM topology_cell WHERE topology_id = m.topology_id AND rank = 2) AS faces, " + ...
    "(SELECT s.array_path FROM selection s JOIN selection p ON p.id = s.member_of_id WHERE p.manifold_id = m.id AND p.name = m.name || '_atlases' AND s.name = 'Desikan-Killiany') AS atlas, " + ...
    "(SELECT e.code FROM selection s JOIN selection p ON p.id = s.member_of_id JOIN dictionary_entry e ON e.dictionary_id = s.dictionary_id " + ...
    " WHERE p.manifold_id = m.id AND p.name = m.name || '_atlases' AND s.name = 'Desikan-Killiany' AND e.name IN ('unassigned', 'unknown') ORDER BY e.name LIMIT 1) AS wall, " + ...
    "(SELECT gs.positions_path FROM manifold ms JOIN geometry gs ON gs.id = ms.geometry_id WHERE ms.subject_id = m.subject_id AND ms.name = m.name || '_sphere') AS sphere " + ...
    "FROM manifold m JOIN geometry g ON g.id = m.geometry_id JOIN subject u ON u.id = m.subject_id " + ...
    "WHERE u.name = " + q(subject) + " AND u.status = 'complete' AND " + pick);
if isempty(r)
    error('io:read:nxrsurface:noSurface', '%s: no %s in %s', subject, i_ifelse(opts.Surface == "", "primary surface", "surface '" + opts.Surface + "'"), store);
end
if isempty(r(1).atlas) || isempty(r(1).wall)
    error('io:read:nxrsurface:noWall', '%s: %s carries no Desikan-Killiany atlas with its medial wall - the hemispheres are glued by anatomy only', subject, r(1).name);
end
z = fullfile(store, r(1).store);
P = rheome.io.read.zarr(fullfile(z, r(1).positions));
F = rheome.io.read.zarr(fullfile(z, r(1).faces));
A = rheome.io.read.zarr(fullfile(z, r(1).atlas));
if size(P,2) ~= 3 || size(F,2) ~= 3 || numel(A) ~= size(P,1)
    error('io:read:nxrsurface:shape', '%s: positions %s, faces %s, atlas %d do not make one surface', subject, mat2str(size(P)), mat2str(size(F)), numel(A));
end
Sp = [];
if ~isempty(r(1).sphere), Sp = double(rheome.io.read.zarr(fullfile(z, r(1).sphere))); end
S = struct('Vertices', double(P), 'Faces', double(F) + 1, 'nV', size(P,1), 'nF', size(F,1), 'Wall', A(:) == r(1).wall, 'Sphere', Sp, ...
    'Comment', r(1).name, 'Subject', char(subject), 'SurfaceFile', fullfile(z, r(1).positions));
end

function v = i_ifelse(c, a, b)
if c, v = a; else, v = b; end
end
