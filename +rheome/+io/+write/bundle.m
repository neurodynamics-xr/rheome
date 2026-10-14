function manifest = bundle(folder, surf, result, varargin)
% IO.WRITE.BUNDLE  Write a result on a surface as a rheome.bundle@1 folder (what rheome.view hands the viewer).
%
%   manifest = rheome.io.write.bundle(folder, surf, result)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Fs', fs, 'Name', name, 'Settings', s)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Potentials', "J")   % + J's Hodge potentials phi, psi
%
% The folder holds raw little-endian arrays and one manifest.json, written with fwrite and nothing else:
%   positions.bin  float32 [3 x nV]  (the vertices, meters)      faces.bin  int32 [3 x nF], 0-based
%   field<k>.bin   float32, the field's [nV x T] ('real') or [3nV x T] ('vector', interleaved x1,y1,z1,...) as is
%   positions_<name>.bin  float32 [3 x nV], each further surface of the same mesh (inflated, sphere ...), optional,
%                         listed in the manifest's 'surfaces'
% The manifest may carry 'settings', the viewer's settings to open with (rheome.view writes them).
% Every file is .bin: MATLAB's uihtml server answers 404 for an extension it does not know (.f32, .i32; R2026a).
% The reader is the Cortical Flow viewer (nxr-cortical-flow editor/src/viewer/bundle.ts; format and protocol in its
% README.md); it finds the hemispheres itself, as the mesh's two largest connected components.
%
% INPUTS:
%   folder : the bundle folder (created; files in it are overwritten)
%   surf   : surface struct with .Vertices [nV x 3] (meters) and .Faces [nF x 3] (1-based), as rheome.io.read.surface,
%            OR a struct of such surfaces of one mesh, by name, e.g. struct('pial', S, 'inflated', Si); the first is
%            the default, the others are written as the manifest's 'surfaces' (same nV; the first's faces)
%   result : a [nV x T] or [3nV x T] array, OR a struct whose fields are such arrays (one bundle field each, by name),
%            e.g. struct('activity', I, 'velocity', v) with v as rheome.dynamics.opticalflow_scalar returns it.
%            The type follows the row count: nV rows 'real', 3nV rows 'vector'. Every field has the same T.
%
% OPTIONS (name/value):
%   'Fs'    sampling rate of the T columns, Hz (default 1)
%   'Name'  the bundle's name (default 'Rheome result')
%   'Settings'  struct of viewer settings, written as the manifest's 'settings' (default: none)
%   'Potentials'  the vector fields whose Hodge potentials are written too: names (string/cellstr), true (every
%               vector field) or false (default). Each adds two real fields, <v>_phi (scalar potential, sources and
%               sinks) and <v>_psi (stream function, vortices), from rheome.differential.helmholtz on the first
%               surface (its .Hemi, if any, sets the per-hemisphere solves), tagged "derivedFrom": "<v>" and
%               "form": "phi" | "psi" in the manifest (additive to rheome.bundle@1). The viewer cannot solve them.
%   'Paths'  PATHS ON THE SURFACE (a trajectory: a geodesic traced over the cortex), the manifest's 'paths' (viewer README
%            'paths[]'): a struct array with
%              name      the path's name (the viewer's 'paths' setting shows paths by name)                   REQUIRED
%              vertices  its points as vertices, [n x 1] 1-based                                      one of these two
%              faces, weights   or points of triangles: faces [n x 1] (1-based rows of the first surface's Faces, the
%                        corners as given) and barycentric weights [n x 3] (>= 0, summing to 1)
%              color     a figure-token colour's name ("atlas-event", "neutral-ink" ...), a neutral or "#rrggbb"
%                        (default: the viewer's, atlas-event)
%              width     pt (default: the viewer's, the tokens' line.event_pt)
%              arrow     true: an arrowhead at the end (default false)
%              geodesic  true: the vertices are waypoints the viewer joins by geodesics (default false)
%            Each is written as path<k>_v.bin (int32, 0-based) or path<k>_f.bin (int32, 0-based) + path<k>_w.bin
%            (float32 [3 x n]).
%
% OUTPUT:
%   manifest : full path of the manifest.json written.
%
% See also: rheome.view
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Fs', 1, @(x) isscalar(x) && x > 0);
    p.addParameter('Name', 'Rheome result', @(x) ischar(x) || isstring(x));
    p.addParameter('Settings', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Potentials', false);
    p.addParameter('Paths', [], @(x) isempty(x) || isstruct(x));
    p.parse(varargin{:});
    o = p.Results;

    if isfield(surf, 'Vertices'), sname = {'default'}; surfs = {surf};
    else, sname = fieldnames(surf)'; surfs = struct2cell(surf)';
    end
    if isempty(surfs) || ~all(cellfun(@(x) isstruct(x) && isfield(x, 'Vertices'), surfs))
        error('rheome:bundle:surface', 'A surface struct, or a struct of surface structs by name, is expected.');
    end
    surf = surfs{1};
    V = surf.Vertices; F = surf.Faces; nV = size(V, 1);
    if size(V, 2) ~= 3 || size(F, 2) ~= 3, error('rheome:bundle:surface', 'Vertices and Faces must be [n x 3].'); end
    if any(F(:) < 1 | F(:) > nV | F(:) ~= round(F(:))), error('rheome:bundle:faces', 'Faces must be 1-based indices into Vertices.'); end
    if isnumeric(result), result = struct('result', result); end
    names = fieldnames(result);
    if isempty(names), error('rheome:bundle:empty', 'The result has no field.'); end
    T = size(result.(names{1}), 2);
    pot = potentials(o.Potentials, names, result, nV);

    if ~isfolder(folder), mkdir(folder); end
    put(folder, 'positions.bin', V', 'single');
    put(folder, 'faces.bin', F' - 1, 'int32');
    sm = struct('positions', 'positions.bin', 'faces', 'faces.bin', 'unit', 'm', 'name', sname{1});
    if numel(surfs) > 1
        alt = struct('name', sname(2:end), 'positions', strcat('positions_', sname(2:end), '.bin'));
        for k = 2:numel(surfs)
            if ~isequal(size(surfs{k}.Vertices), [nV 3]), error('rheome:bundle:surface', '%s: %d vertices, want nV = %d.', sname{k}, size(surfs{k}.Vertices, 1), nV); end
            put(folder, alt(k-1).positions, surfs{k}.Vertices', 'single');
        end
    end
    fl = cell(1, numel(names));
    for k = 1:numel(names)
        x = result.(names{k});
        if ~isnumeric(x) || ~isreal(x) || ~ismatrix(x), error('rheome:bundle:field', '%s: a real numeric matrix is expected.', names{k}); end
        if size(x, 1) == nV, type = 'real';
        elseif size(x, 1) == 3*nV, type = 'vector';
        else, error('rheome:bundle:rows', '%s: %d rows, want nV = %d (real) or 3nV = %d (vector).', names{k}, size(x, 1), nV, 3*nV);
        end
        if size(x, 2) ~= T, error('rheome:bundle:time', '%s: %d columns, want T = %d (the first field''s).', names{k}, size(x, 2), T); end
        file = sprintf('field%d.bin', k);
        put(folder, file, x, 'single');
        fl{k} = struct('name', names{k}, 'type', type, 'values', file, 'time', true);
    end
    for k = 1:numel(pot)
        S = surf;
        if ~isfield(S, 'VertNormals') || isempty(S.VertNormals)
            fg = rheome.operators.face_gradient(V, F);  S.VertNormals = fg.W * fg.FaceNormal;
        end
        H = rheome.differential.helmholtz(double(result.(pot{k})), S);
        for f = {'phi', 'Phi'; 'psi', 'Psi'}'
            name = [pot{k} '_' f{1}];
            if any(strcmp(name, names)), error('rheome:bundle:potentials', '%s: the result already has a field of that name.', name); end
            n = numel(fl) + 1;  file = sprintf('field%d.bin', n);
            put(folder, file, H.(f{2}), 'single');
            fl{n} = struct('name', name, 'type', 'real', 'values', file, 'time', true, 'derivedFrom', pot{k}, 'form', f{1});
        end
    end
    m = struct('format', 'rheome.bundle@1', 'name', char(o.Name), ...
        'surface', sm, 'time', struct('samples', T, 'fs', o.Fs, 'origin', 0), 'fields', {fl});
    if numel(surfs) > 1, m.surfaces = num2cell(alt); end   % a JSON array even for one
    if ~isempty(o.Paths), m.paths = writePaths(folder, o.Paths, nV, size(F, 1)); end
    if ~isempty(o.Settings), m.settings = o.Settings; end
    manifest = fullfile(folder, 'manifest.json');
    fid = fopen(manifest, 'w');
    if fid < 0, error('rheome:bundle:write', 'Cannot write %s.', manifest); end
    fwrite(fid, jsonencode(m)); fclose(fid);
end

function pot = potentials(x, names, result, nV)
% the vector fields named by 'Potentials' (true: every one)
    vec = names(cellfun(@(f) size(result.(f), 1) == 3*nV, names));
    if (islogical(x) || isnumeric(x)) && isscalar(x)
        pot = {};  if x, pot = vec'; end
        return
    end
    if ~(ischar(x) || iscellstr(x) || isstring(x)), error('rheome:bundle:potentials', 'Potentials: names, true or false.'); end
    pot = cellstr(x);
    bad = setdiff(pot, vec);
    if ~isempty(bad), error('rheome:bundle:potentials', 'Potentials: %s is not a vector field of the result.', strjoin(bad, ', ')); end
end

function P = writePaths(folder, a, nV, nF)
% the paths, checked (as the viewer reads them: bundle.ts readPath) and written; a cell, a JSON array even for one
    id = 'rheome:bundle:paths';
    bad = setdiff(fieldnames(a), {'name', 'vertices', 'faces', 'weights', 'color', 'width', 'arrow', 'geodesic'});
    if ~isempty(bad), error(id, 'Paths: %s is not a path field.', strjoin(bad, ', ')); end
    P = cell(1, numel(a));
    for k = 1:numel(a)
        b = a(k);
        if ~isfield(b, 'name') || isempty(b.name), error(id, 'Paths(%d): a name.', k); end
        q = struct('name', char(b.name));
        if isfield(b, 'vertices') && ~isempty(b.vertices)
            v = double(b.vertices(:));
            if numel(v) < 2 || any(v < 1 | v > nV | v ~= round(v)), error(id, '%s: vertices, >= 2 of them in 1..%d.', q.name, nV); end
            q.vertices = sprintf('path%d_v.bin', k);  put(folder, q.vertices, v - 1, 'int32');
        elseif isfield(b, 'faces') && isfield(b, 'weights') && ~isempty(b.faces)
            f = double(b.faces(:));  w = double(b.weights);
            if numel(f) < 2 || any(f < 1 | f > nF | f ~= round(f)), error(id, '%s: faces, >= 2 of them in 1..%d.', q.name, nF); end
            if ~isequal(size(w), [numel(f) 3]) || any(w(:) < -1e-6) || any(abs(sum(w, 2) - 1) > 1e-6)
                error(id, '%s: weights [n x 3], barycentric (>= 0, each row summing to 1).', q.name);
            end
            q.faces = sprintf('path%d_f.bin', k);  put(folder, q.faces, f - 1, 'int32');
            q.weights = sprintf('path%d_w.bin', k);  put(folder, q.weights, w', 'single');
        else
            error(id, '%s: vertices, or faces and weights.', q.name);
        end
        if isfield(b, 'color') && ~isempty(b.color)
            c = char(b.color);
            if isempty(regexp(c, '^(#[0-9a-fA-F]{6}|[a-z][a-z0-9-]*)$', 'once')), error(id, '%s: color, a colour token''s name or #rrggbb.', q.name); end
            q.color = c;
        end
        if isfield(b, 'width') && ~isempty(b.width)
            if ~(isscalar(b.width) && b.width > 0 && b.width <= 20), error(id, '%s: width in pt, 0 < width <= 20.', q.name); end
            q.width = double(b.width);
        end
        for t = {'arrow', 'geodesic'}
            if isfield(b, t{1}) && ~isempty(b.(t{1})), q.(t{1}) = logical(b.(t{1})); end
        end
        P{k} = q;
    end
end

function put(folder, name, x, cls)
    f = fullfile(folder, name);
    fid = fopen(f, 'w', 'ieee-le');
    if fid < 0, error('rheome:bundle:write', 'Cannot write %s.', f); end
    fwrite(fid, x, cls); fclose(fid);
end
