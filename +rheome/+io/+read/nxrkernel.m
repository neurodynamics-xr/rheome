function K = nxrkernel(store, subject, kernel)
% IO.READ.NXRKERNEL  Read an imaging kernel from an nxr store, with its channels, the recording it applies to and its stored fused weights.
%
%   K = rheome.io.read.nxrkernel(store, subject, kernel)
%
% A kernel is an operator of kind 'inverse kernel': [3V x M] with a free orientation (x, y, z of each primary-cortex vertex,
% interleaved) or [V x M] with a fixed one, its M columns the channels <kernel>_channels names (0-based members of the
% recording's channel set). It applies to the raw recording as stored: the editor applies no whitener or noise covariance of
% its own. Its recording is the editor's choice (sourceMappingSource): a recording on its channel set, the one of its session
% first. Its FUSED WEIGHTS B (the editor's hodgeSynthesisWeights, nxr-convert fused_weights) are the operator of kind
% 'hodge-dirac weights' whose params name this kernel: [2K x M], rows 1..K the normal-current weights a, K+1..2K the flow
% weights b, both in the mode order of the basis they were fused with (whose eigenvalues mu and glued mesh key they record).
%
% OUTPUT (struct K)
%   .W          [n_rows x M] single, the kernel;  .Free  true for 3 rows a vertex;  .Unit  'A·m' for a current estimate, else ''
%   .Channels   [M x 1] the recording rows (1-based) its columns read
%   .Recording  the name of the recording it applies to ('' when none is stored);  .Session, .Name, .Id, .Path
%   .Weights    struct array of its fused weights (empty when none): .B [2K x M] single, .K, .mu [K x 1], .mesh (the glued
%               mesh key, rheome.geom.meshkey), .Name, .Path — one per basis it was fused with (its own surface's, a template's, ...)
%
% See also: rheome.io.read.nxrrecording, rheome.io.read.nxrhodge, rheome.geom.meshkey
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, kernel (1,1) string, end
q = @(s) "'" + replace(s, "'", "''") + "'";
o = nxrsql(store, "SELECT o.id, o.name, o.path, o.session, o.n_rows, o.n_cols, o.to_value_type, o.params_json, o.from_manifold_id, " + ...
    "s.array_path AS channels, u.path AS store FROM operator o JOIN subject u ON u.id = o.subject_id LEFT JOIN selection s ON s.id = o.from_selection_id " + ...
    "WHERE u.name = " + q(subject) + " AND u.status = 'complete' AND o.status = 'complete' AND o.kind = 'inverse kernel' AND o.path IS NOT NULL " + ...
    "AND o.to_manifold_id IS NOT NULL AND o.name = " + q(kernel));
if isempty(o), error('io:read:nxr:notFound', '%s: no inverse kernel ''%s''', subject, kernel); end
if isempty(o.channels), error('io:read:nxr:layout', 'kernel ''%s'': no channel selection', kernel); end
z = fullfile(store, o.store);
W = rheome.io.read.zarr(fullfile(z, o.path));
ch = double(rheome.io.read.zarr(fullfile(z, o.channels))) + 1;
if size(W, 2) ~= numel(ch), error('io:read:nxr:layout', 'kernel ''%s'': %d columns for %d channels', kernel, size(W, 2), numel(ch)); end
est = '';
if ~isempty(o.params_json), p = jsondecode(o.params_json); if isfield(p, 'inverse') && isfield(p.inverse, 'estimate'), est = p.inverse.estimate; end, end
r = nxrsql(store, "SELECT f.name, f.session FROM field f JOIN geometry_dimension d ON d.ordinal = 0 JOIN manifold p ON p.id = f.manifold_id AND p.geometry_id = d.geometry_id " + ...
    "WHERE f.subject_id = (SELECT subject_id FROM operator WHERE id = " + q(o.id) + ") AND f.kind = 'recording' AND f.status = 'complete' " + ...
    "AND d.manifold_id = " + q(o.from_manifold_id) + " ORDER BY f.rowid");
rec = '';
if ~isempty(r)
    k = find(strcmp({r.session}, char(i_str(o.session))), 1); if isempty(k), k = 1; end
    rec = r(k).name;
end
b = nxrsql(store, "SELECT name, path, params_json FROM operator WHERE kind = 'hodge-dirac weights' AND status = 'complete' AND path IS NOT NULL " + ...
    "AND json_extract(params_json, '$.kernel') = " + q(o.id));
Wt = struct('B', {}, 'K', {}, 'mu', {}, 'mesh', {}, 'Name', {}, 'Path', {});
for i = 1:numel(b)
    p = jsondecode(b(i).params_json);
    B = rheome.io.read.zarr(fullfile(z, b(i).path));
    if size(B, 2) ~= numel(ch) || size(B, 1) ~= 2 * p.K || numel(p.mu) ~= p.K
        error('io:read:nxr:layout', '%s: [%d x %d] for K = %d and %d channels', b(i).name, size(B, 1), size(B, 2), p.K, numel(ch));
    end
    Wt(end+1) = struct('B', B, 'K', p.K, 'mu', p.mu(:), 'mesh', i_str(p.mesh), 'Name', b(i).name, 'Path', fullfile(z, b(i).path)); %#ok<AGROW>
end
K = struct('W', W, 'Free', strcmp(i_str(o.to_value_type), 'vector3'), 'Unit', i_ifelse(strcmp(est, 'current'), 'A·m', ''), ...
    'Channels', ch, 'Recording', rec, 'Session', i_str(o.session), 'Name', o.name, 'Id', o.id, 'Path', fullfile(z, o.path), 'Weights', Wt);
end

function v = i_str(x)
if isempty(x), v = ''; else, v = char(x); end
end

function v = i_ifelse(c, a, b)
if c, v = a; else, v = b; end
end
