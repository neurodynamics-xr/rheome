function d = of(x, opts)
% DOMAIN.OF  A domain descriptor for a mesh, a graph or a point set: what it is, how big, whose.
%
%   d = rheome.domain.of(S)                       % a surface struct (rheome.io.read.surface)
%   d = rheome.domain.of(S, Name="sub01", Kind="complex")
%   d = rheome.domain.of(G)                       % a sensor graph (rheome.sensors.graph)
%
% ⭐ WHY A DESCRIPTOR AND NOT THE MESH. Two fields are on the same domain when their
% domains have the same ID, and an operator is applicable when its declared input domain
% kind and element counts match. Neither question needs the coordinates, so the descriptor
% is small enough to travel with a field, to be compared, and to be written into a store.
%
% OUTPUT (struct d):
%   .kind      'complex' | 'graph' | 'points' (rheome.domain.kinds)
%   .nV .nE .nF   element counts; nE is derived from the faces when the mesh has no edges
%   .id        8 hex characters of MD5 over the TOPOLOGY (kind, counts, sorted edges) --
%              two hemispheres of the same cortex differ, a re-import of the same surface
%              does not
%   .name      a human handle ('' when the caller does not give one)
%   .source    the file the mesh came from, when the struct carries one (provenance)
%
% ⚠ THE ID IS OF THE TOPOLOGY, NOT THE GEOMETRY. An inflated and a white surface of the
% same subject share connectivity and therefore share an ID: a field defined on one IS
% defined on the other, which is exactly the property the flow code relies on when it
% visualises on an inflated mesh. Say it with .name when the geometry matters.
%
% See also: rheome.domain.kinds, rheome.domain.index, rheome.fieldtype.validate
%
% Author: Diellor Basha, 2026

    arguments
        x
        opts.Name (1,1) string = ""
        opts.Kind (1,1) string = ""
    end
    if isstruct(x) && isfield(x, 'kind') && isfield(x, 'id')
        d = x;  return                                        % already a descriptor
    end
    if ~isstruct(x)
        error('domain:of:input', 'Give a surface or graph struct; use rheome.domain.index for a line.');
    end

    V = i_get(x, 'Vertices', []);  F = i_get(x, 'Faces', []);
    nV = i_get(x, 'nV', size(V, 1));
    nF = i_get(x, 'nF', size(F, 1));
    if opts.Kind ~= ""
        kind = char(opts.Kind);
    elseif ~isempty(F)
        kind = 'complex';
    elseif ~isempty(i_get(x, 'W', [])) || ~isempty(i_get(x, 'Laplacian', []))
        kind = 'graph';
    else
        kind = 'points';
    end
    if ~ismember(string(kind), rheome.domain.kinds().kind)
        error('domain:of:kind', 'Unknown domain kind ''%s''.', kind);
    end

    E = i_edges(F);
    d = struct('kind', kind, 'nV', double(nV), 'nE', size(E, 1), 'nF', double(nF), ...
               'id', '', 'name', char(opts.Name), 'source', char(string(i_get(x, 'SurfaceFile', ''))));
    d.id = i_hash(sprintf('%s|%d|%d|%d|', kind, d.nV, d.nE, d.nF), E);
end

% The undirected edge set of a triangle list, sorted, unique: the topology an ID is of.
function E = i_edges(F)
    if isempty(F), E = zeros(0, 2); return; end
    E = [F(:, [1 2]); F(:, [2 3]); F(:, [3 1])];
    E = unique(sort(double(E), 2), 'rows');
end

function h = i_hash(head, E)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(head));
    if ~isempty(E), md.update(typecast(uint32(E(:)'), 'uint8')); end
    dg = typecast(md.digest(), 'uint8');
    h = lower(reshape(dec2hex(dg(1:4), 2)', 1, []));
end

function v = i_get(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
% Author: Diellor Basha, 2026
