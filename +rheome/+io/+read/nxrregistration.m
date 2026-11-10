function G = nxrregistration(store, subject, template)
% IO.READ.NXRREGISTRATION  Read a subject's stored registration to a template from an nxr dataset: T and the record it was built with.
%
%   G = rheome.io.read.nxrregistration(store, subject)               % the subject's only registration
%   G = rheome.io.read.nxrregistration(store, subject, template)     % the one to this template (its name, without the '@')
%
% The registration is the operator hodgeDirac_to_<template> (kind 'registration', stored on the subject by `nxr register build`):
% T [K_template x K_subject], the sphere registration written in the two Hodge-Dirac mode bases (rheome.register.mode_registration),
% from the subject's stored basis (rheome.io.read.nxrhodge) to the template's. Followed by its rows, as the Cortical Flow Editor's
% hodgeStore reads it (hodgeTransferStoreSource) and nxr's store_mode_registration writes it.
%
% OUTPUT (struct G)
%   .T         [K_template x K_subject] double (stored as float32)
%   .Template  the template's name;  .Name  the operator's name;  .Id, .Path
%   .Sampling  "barycentric", "nearest" or "shepard";  .Neighbours, .Power  (shepard only, else [])
%   .Flow      how a flow mode is carried ("fields": as flow fields of unit strength, rheome.inverse.hodge_register with mu)
%   .From, .To structs of the two packages the registration was built from: .Package (the sha256 of the modes), .K, .mu [K x 1],
%              .Mesh (the glued mesh key); .To also .Kind [K x 1] (0 exact, 2 coexact)
%
% See also: rheome.register.mode_registration, rheome.inverse.hodge_register, rheome.io.read.nxrhodge
%
% Author: Diellor Basha, 2026
arguments, store (1,1) string, subject (1,1) string, template (1,1) string = "", end
q = @(s) "'" + replace(s, "'", "''") + "'";
who = "(SELECT id FROM subject WHERE name = " + q(subject) + " AND status = 'complete')";
o = nxrsql(store, "SELECT o.id, o.name, o.path, o.n_rows, o.n_cols, o.params_json, u.path AS store FROM operator o JOIN subject u ON u.id = o.subject_id " + ...
    "JOIN manifold a ON a.id = o.from_manifold_id WHERE o.subject_id = " + who + " AND o.kind = 'registration' AND o.status = 'complete' AND o.path IS NOT NULL " + ...
    "AND substr(o.name, 1, 14) = 'hodgeDirac_to_' AND a.name = 'hodgeDirac' AND a.type = 'points' ORDER BY o.name");
if isempty(o), error('io:read:nxr:notFound', '%s: no stored registration (nxr register build)', subject); end
names = string({o.name});
if template ~= ""
    k = find(names == "hodgeDirac_to_" + erase(template, "@"), 1);
    if isempty(k), error('io:read:nxr:notFound', '%s: no registration to ''%s'' (stored: %s)', subject, template, strjoin(extractAfter(names, 14), ', ')); end
elseif numel(o) > 1
    error('io:read:nxr:ambiguous', '%s: registrations to %s; name one', subject, strjoin(extractAfter(names, 14), ', '));
else
    k = 1;
end
o = o(k); p = jsondecode(o.params_json);
z = fullfile(store, o.store);
T = double(rheome.io.read.zarr(fullfile(z, o.path)));
if ~isequal(size(T), [o.n_rows o.n_cols]) || size(T, 1) ~= p.to.K || size(T, 2) ~= p.from.K
    error('io:read:nxr:layout', '%s: T is %s, its record says %d template and %d subject modes', o.name, mat2str(size(T)), p.to.K, p.from.K);
end
G = struct('T', T, 'Template', char(extractAfter(string(o.name), 14)), 'Name', o.name, 'Id', o.id, 'Path', fullfile(z, o.path), ...
    'Sampling', string(p.sampling), 'Neighbours', i_opt(p, 'neighbours'), 'Power', i_opt(p, 'power'), 'Flow', string(i_opt(p, 'flow', "")), ...
    'From', i_package(p.from), 'To', i_package(p.to));
end

function v = i_opt(p, f, d)
if nargin < 3, d = []; end
if isfield(p, f) && ~isempty(p.(f)), v = p.(f); else, v = d; end
end

function s = i_package(x)
s = struct('Package', string(i_opt(x, 'package', "")), 'K', x.K, 'mu', x.mu(:), 'Mesh', string(i_opt(x, 'mesh', "")));
if isfield(x, 'kind'), s.Kind = x.kind(:); end
end
