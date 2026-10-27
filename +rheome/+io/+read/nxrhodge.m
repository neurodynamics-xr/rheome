function H = nxrhodge(store, subject, opts)
% IO.READ.NXRHODGE  Read a stored Hodge-Dirac basis (nxr-convert hodge) from an nxr store.
%
%   H = rheome.io.read.nxrhodge(store, subject)                   % the 'hodgeDirac' basis, every mode
%   H = rheome.io.read.nxrhodge(store, subject, Modes=200, Name="hodgeDirac")
%
% A basis is the points axis <name> (its positions the eigenvalues mu, ascending), the field <name>_kind (0 exact: the
% potential psi/sqrt(2); 2 coexact), the selection <name>_support (the primary-cortex vertices kept) and the operator
% <name>_eigenvectors [V x K]. When that operator's surface records a glued mesh (its producer_json 'mesh'), V is the glued
% surface's vertices (the support, then the bridge) and the surface is returned with it; otherwise V is the primary cortex,
% zero off the support. The lookups are the editor's (adapters/nxr hodgeStore.ts). The mass is not stored: rebuild it on the
% surface (rheome.operators.hodge_dirac).
%
% OPTIONS: Name the basis' name (default "hodgeDirac");  Modes the first k modes only (default: all)
%
% OUTPUT (struct H)
%   .mu [K x 1] double, .kind [K x 1], .P [V x K] single (one vertex potential a mode), .support [n x 1] 1-based
%   .glued  true when P is on a glued surface; then .Vertices [V x 3], .Faces [nF x 3] 1-based and .mesh (its key)
%   .Surface the name of the surface P is on;  .Name
%
% See also: rheome.io.read.nxrkernel, rheome.eigen.hodge_dirac, rheome.geom.glue
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, opts.Name (1,1) string = "hodgeDirac", opts.Modes double = [], end
q = @(s) "'" + replace(s, "'", "''") + "'";
n = opts.Name; who = "(SELECT id FROM subject WHERE name = " + q(subject) + " AND status = 'complete')";
r = nxrsql(store, "SELECT (SELECT path FROM subject WHERE id = " + who + ") AS store, " + ...
    "(SELECT g.positions_path FROM manifold m JOIN geometry g ON g.id = m.geometry_id WHERE m.subject_id = " + who + " AND m.name = " + q(n) + " AND m.type = 'points') AS mu, " + ...
    "(SELECT path FROM field WHERE subject_id = " + who + " AND name = " + q(n + "_kind") + ") AS kind, " + ...
    "(SELECT array_path FROM selection WHERE subject_id = " + who + " AND name = " + q(n + "_support") + ") AS support, " + ...
    "o.path AS P, o.n_rows, o.n_cols, t.name AS surface, json_extract(t.producer_json, '$.mesh') AS mesh, g.positions_path AS positions, " + ...
    "(SELECT path FROM topology_cell WHERE topology_id = t.topology_id AND rank = 2) AS faces " + ...
    "FROM operator o JOIN manifold t ON t.id = o.to_manifold_id LEFT JOIN geometry g ON g.id = t.geometry_id " + ...
    "WHERE o.subject_id = " + who + " AND o.name = " + q(n + "_eigenvectors") + " AND o.kind = 'eigenvectors' AND o.status = 'complete'");
if isempty(r) || isempty(r.mu) || isempty(r.kind) || isempty(r.support)
    error('io:read:nxr:notFound', '%s: no stored Hodge-Dirac basis ''%s''', subject, n);
end
z = fullfile(store, r.store);
K = r.n_cols; if ~isempty(opts.Modes), K = min(opts.Modes, K); end
H = struct('mu', double(rheome.io.read.zarr(fullfile(z, r.mu), {[1 K]})), 'kind', double(rheome.io.read.zarr(fullfile(z, r.kind), {[1 K]})), ...
    'P', rheome.io.read.zarr(fullfile(z, r.P), {[], [1 K]}), 'support', double(rheome.io.read.zarr(fullfile(z, r.support))) + 1, ...
    'glued', ~isempty(r.mesh), 'Vertices', [], 'Faces', [], 'mesh', '', 'Surface', r.surface, 'Name', char(n));
if H.glued
    H.Vertices = double(rheome.io.read.zarr(fullfile(z, r.positions)));
    H.Faces = double(rheome.io.read.zarr(fullfile(z, r.faces))) + 1;
    H.mesh = char(r.mesh);
    if ~strcmp(rheome.geom.meshkey(H.Faces), H.mesh), error('io:read:nxr:layout', '%s: the glued faces do not hash to the recorded key %s', n, H.mesh); end
end
end
