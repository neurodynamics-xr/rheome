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
%   codes_<k>.bin  int32 [nV], labelling k's code a vertex (< 0: none)          optional, the manifest's 'labellings'
%   points_<k>_{vertices,frames,signs}.bin  int32 [n], point set k (0-based)      optional, the manifest's 'points'
%   path_<k>_vertices.bin  int32 [n] (0-based), OR path_<k>_faces.bin int32 [n] (0-based) and path_<k>_weights.bin
%                          float32 [3 x n] (a point's three barycentric weights after another)   optional, the manifest's 'paths'
%   selection_<k>.bin  int32 [n], vertex set k (0-based)                            optional, the manifest's 'selections'
% The manifest may carry 'settings', the viewer's settings to open with (rheome.view writes them).
% Every file is .bin: MATLAB's uihtml server answers 404 for an extension it does not know (.f32, .i32; R2026a).
% The reader is the Cortical Flow viewer (the Cortical Flow Editor's viewer/bundle.ts; format and protocol in its
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
%   'Labellings'  struct array (or cell of structs), one partition each: .name, .codes [nV x 1] integers (< 0: none),
%               optional .entries struct array (.code .name .color '#rrggbb'). rheome.viz.tiles makes one from a
%               rheome.geom.tree level. The viewer paints it in the colour layer's place (setting 'labelling').
%   'Points'    struct array (or cell of structs), one point set each: .name, .vertices [n x 1] (1-BASED), optional
%               .frames [n x 1] (1-based frame each point is shown at; 0: every frame; none: every frame), .signs
%               [n x 1] (+1 / -1, drawn as plus / minus; none: +1), .color '#rrggbb'. rheome.viz.points makes one
%               from rheome.detect.criticalPoints / rheome.flow.critical.
%   'Paths'     struct array (or cell of structs), one PATH ON THE SURFACE each (drawn as a depth-tested line; shown by the
%               viewer setting 'paths', default all): .name, and its points ON the mesh, either .vertices [n x 1] (1-BASED,
%               n >= 2) or .faces [n x 1] (1-BASED rows of the first surface's Faces) with .weights [n x 3] (barycentric:
%               >= 0, summing to 1; the corners in Faces' order), so the path follows every embedding. Optional:
%               .color (a figure-token colour's name, e.g. 'atlas-event' (the viewer's default), 'neutral-ink', a neutral
%               such as 'ink', or '#rrggbb'), .width (pt, 0 < width <= 20), .geodesic (true: the vertices are WAYPOINTS
%               the viewer joins by geodesics on the folded surface; vertices only), .arrow (true: an arrowhead at the
%               end), .ticks (a geodesic path as a time line: struct with .every_mm (> 0) or .count (>= 2), one of
%               them, and optionally .labels, texts). The viewer's README (the Cortical Flow Editor's viewer) says how
%               each is drawn.
%   'Selections' struct array (or cell of structs), one VERTEX SET each (a tile, a mask, a support; shown by the
%               viewer setting 'selections', default all): .name, .vertices [n x 1] (1-BASED) or a logical [nV] mask;
%               optional .style 'fill' (the viewer's default: a blended fill) | 'outline' (its edge in its colour) |
%               'extent' (a dashed ink line, never a fill), .color (as a path's), .width (pt, 0 < width <= 20).
%   Names are non-empty and unique within Paths and within Selections; keys not listed above are refused. Only what
%   is given is written.
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
    p.addParameter('Labellings', []);
    p.addParameter('Points', []);
    p.addParameter('Paths', []);
    p.addParameter('Selections', []);
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
    if ~isempty(o.Labellings), m.labellings = labellings(folder, o.Labellings, nV); end
    if ~isempty(o.Points), m.points = points(folder, o.Points, nV, T); end
    if ~isempty(o.Paths), m.paths = paths(folder, o.Paths, nV, size(F, 1)); end
    if ~isempty(o.Selections), m.selections = selections(folder, o.Selections, nV); end
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

function out = labellings(folder, L, nV)
% each labelling's codes written, its manifest entry returned (a cell: a JSON array even for one)
    L = asCell(L, 'Labellings');  out = cell(1, numel(L));
    for k = 1:numel(L)
        x = L{k};  need(x, {'name', 'codes'}, 'Labellings');
        c = double(x.codes(:));
        if numel(c) ~= nV || any(c ~= round(c)), error('rheome:bundle:labellings', '%s: codes must be %d integers (one a vertex).', char(x.name), nV); end
        file = sprintf('codes_%d.bin', k);  put(folder, file, c, 'int32');
        out{k} = struct('name', char(x.name), 'codes', file);
        if isfield(x, 'entries') && ~isempty(x.entries), out{k}.entries = num2cell(x.entries(:)'); end
    end
end

function out = points(folder, P, nV, T)
% each point set's vertices (0-based), frames (0-based, -1 every) and signs written, its manifest entry returned
    P = asCell(P, 'Points');  out = cell(1, numel(P));
    for k = 1:numel(P)
        x = P{k};  need(x, {'name', 'vertices'}, 'Points');
        v = double(x.vertices(:));  n = numel(v);
        if any(v < 1 | v > nV | v ~= round(v)), error('rheome:bundle:points', '%s: vertices must be 1-based indices into the %d vertices.', char(x.name), nV); end
        pre = sprintf('points_%d_', k);
        put(folder, [pre 'vertices.bin'], v - 1, 'int32');
        out{k} = struct('name', char(x.name), 'vertices', [pre 'vertices.bin']);
        if isfield(x, 'frames') && ~isempty(x.frames)
            f = double(x.frames(:));
            if isscalar(f), f = repmat(f, n, 1); end
            if numel(f) ~= n || any(f < 0 | f > T | f ~= round(f)), error('rheome:bundle:points', '%s: frames must be %d integers in 0 (every frame) .. %d.', char(x.name), n, T); end
            put(folder, [pre 'frames.bin'], f - 1, 'int32');   % 0 (every) -> -1
            out{k}.frames = [pre 'frames.bin'];
        end
        if isfield(x, 'signs') && ~isempty(x.signs)
            g = sign(double(x.signs(:)));
            if numel(g) ~= n || any(g == 0), error('rheome:bundle:points', '%s: signs must be %d values of +1 or -1.', char(x.name), n); end
            put(folder, [pre 'signs.bin'], g, 'int32');
            out{k}.signs = [pre 'signs.bin'];
        end
        if isfield(x, 'color') && ~isempty(x.color), out{k}.color = char(x.color); end
    end
end

function out = paths(folder, P, nV, nF)
% each path's points written (vertices, or faces and weights; 0-based), its manifest entry returned
    id = 'rheome:bundle:paths';
    P = asCell(P, 'Paths');  out = cell(1, numel(P));  seen = {};
    for k = 1:numel(P)
        x = P{k};  need(x, {'name'}, 'Paths');
        nm = itemName(x.name, 'Paths', seen);  seen{end+1} = nm; %#ok<AGROW>
        known(x, {'name', 'vertices', 'faces', 'weights', 'color', 'width', 'geodesic', 'arrow', 'ticks'}, nm, 'Paths');
        pre = sprintf('path_%d_', k);
        e = struct('name', nm);
        hasV = given(x, 'vertices');
        if hasV == (given(x, 'faces') || given(x, 'weights')) || (~hasV && ~(given(x, 'faces') && given(x, 'weights')))
            error(id, '%s: give vertices, or faces and weights (not both).', nm);
        end
        if hasV
            v = double(x.vertices(:));
            if numel(v) < 2 || any(v < 1 | v > nV | v ~= round(v)), error(id, '%s: vertices must be 2 or more 1-based indices into the %d vertices.', nm, nV); end
            put(folder, [pre 'vertices.bin'], v - 1, 'int32');
            e.vertices = [pre 'vertices.bin'];
        else
            f = double(x.faces(:));  n = numel(f);
            if n < 2 || any(f < 1 | f > nF | f ~= round(f)), error(id, '%s: faces must be 2 or more 1-based indices into the %d faces.', nm, nF); end
            w = double(x.weights);
            if ~isequal(size(w), [n 3]) || any(~(w(:) >= -1e-3)) || any(~(abs(sum(w, 2) - 1) < 1e-3))
                error(id, '%s: weights must be [n x 3], one row a face (n = %d), barycentric (>= 0, summing to 1).', nm, n);
            end
            put(folder, [pre 'faces.bin'], f - 1, 'int32');
            put(folder, [pre 'weights.bin'], w', 'single');
            e.faces = [pre 'faces.bin'];  e.weights = [pre 'weights.bin'];
        end
        if given(x, 'color'), e.color = colour(x.color, nm, id); end
        if given(x, 'width'), e.width = lineWidth(x.width, nm, id); end
        geo = false;
        if given(x, 'geodesic')
            geo = flag(x.geodesic, nm, 'geodesic', id);  e.geodesic = geo;
            if geo && ~hasV, error(id, '%s: a geodesic path is given by vertices (its waypoints).', nm); end
        end
        if given(x, 'arrow'), e.arrow = flag(x.arrow, nm, 'arrow', id); end
        if given(x, 'ticks')
            if ~geo, error(id, '%s: ticks are for a geodesic path (geodesic: true).', nm); end
            e.ticks = ticks(x.ticks, nm, id);
        end
        out{k} = e;
    end
end

function out = selections(folder, X, nV)
% each vertex set's vertices written (0-based), its manifest entry returned
    id = 'rheome:bundle:selections';
    X = asCell(X, 'Selections');  out = cell(1, numel(X));  seen = {};
    for k = 1:numel(X)
        x = X{k};  need(x, {'name', 'vertices'}, 'Selections');
        nm = itemName(x.name, 'Selections', seen);  seen{end+1} = nm; %#ok<AGROW>
        known(x, {'name', 'vertices', 'style', 'color', 'width'}, nm, 'Selections');
        v = x.vertices;
        if islogical(v) && numel(v) == nV, v = find(v(:));
        elseif islogical(v), v = -1;          % a mask is [nV]
        else, v = double(v(:));
        end
        if isempty(v) || any(v < 1 | v > nV | v ~= round(v)), error(id, '%s: vertices must be 1-based indices into the %d vertices, or a mask [nV].', nm, nV); end
        file = sprintf('selection_%d.bin', k);
        put(folder, file, v - 1, 'int32');
        e = struct('name', nm, 'vertices', file);
        if given(x, 'style')
            st = x.style;
            if ~((ischar(st) || (isstring(st) && isscalar(st))) && any(strcmp(char(st), {'fill', 'outline', 'extent'})))
                error(id, '%s: style must be fill, outline or extent.', nm);
            end
            e.style = char(st);
        end
        if given(x, 'color'), e.color = colour(x.color, nm, id); end
        if given(x, 'width'), e.width = lineWidth(x.width, nm, id); end
        out{k} = e;
    end
end

function t = given(x, f)
    t = isfield(x, f) && ~isempty(x.(f));
end

function nm = itemName(x, name, seen)
% a path's or vertex set's name: a non-empty text, once
    if ~((ischar(x) && (isrow(x) || isempty(x))) || (isstring(x) && isscalar(x))) || strlength(string(x)) == 0
        error(['rheome:bundle:' lower(name)], '%s: each name is a non-empty text.', name);
    end
    nm = char(x);
    if any(strcmp(nm, seen)), error(['rheome:bundle:' lower(name)], '%s: %s is named twice.', name, nm); end
end

function known(x, keys, nm, name)
    bad = setdiff(fieldnames(x), keys);
    if ~isempty(bad), error(['rheome:bundle:' lower(name)], '%s: %s: not a key of %s (%s).', nm, strjoin(bad, ', '), name, strjoin(keys, ', ')); end
end

function c = colour(x, nm, id)
% a figure-token colour's name (a neutral among them) or #rrggbb, as the viewer reads it
    if ~((ischar(x) && isrow(x)) || (isstring(x) && isscalar(x))) || isempty(regexp(char(x), '^(#[0-9a-fA-F]{6}|[a-z][a-z0-9-]*)$', 'once'))
        error(id, '%s: color must be a colour token''s name or #rrggbb.', nm);
    end
    c = char(x);
end

function w = lineWidth(x, nm, id)
    if ~(isnumeric(x) && isreal(x) && isscalar(x) && x > 0 && x <= 20), error(id, '%s: width must be a line width in pt, 0 < width <= 20.', nm); end
    w = double(x);
end

function b = flag(x, nm, what, id)
    if ~(isscalar(x) && (islogical(x) || (isnumeric(x) && (x == 0 || x == 1)))), error(id, '%s: %s must be true or false.', nm, what); end
    b = logical(x);
end

function t = ticks(x, nm, id)
% a geodesic's ticks: every_mm or count (one of them), optionally labels
    msg = '%s: ticks: every_mm (mm > 0) or count (2 or more), one of them, and optionally labels (texts).';
    if ~(isstruct(x) && isscalar(x)) || ~isempty(setdiff(fieldnames(x), {'every_mm', 'count', 'labels'})) || isfield(x, 'every_mm') == isfield(x, 'count')
        error(id, msg, nm);
    end
    if isfield(x, 'every_mm')
        e = x.every_mm;
        if ~(isnumeric(e) && isreal(e) && isscalar(e) && e > 0), error(id, msg, nm); end
        t = struct('every_mm', double(e));
    else
        c = x.count;
        if ~(isnumeric(c) && isreal(c) && isscalar(c) && c >= 2 && c == round(c)), error(id, msg, nm); end
        t = struct('count', double(c));
    end
    if isfield(x, 'labels')
        L = x.labels;
        if ischar(L) && (isrow(L) || isempty(L)), L = {L}; end
        if isstring(L), L = cellstr(L); end
        if ~iscellstr(L), error(id, msg, nm); end
        t.labels = reshape(L, 1, []);
    end
end

function c = asCell(x, name)
    if isstruct(x), c = num2cell(x(:)'); elseif iscell(x), c = x(:)'; else, c = {}; end
    if isempty(c) || ~all(cellfun(@(e) isstruct(e) && isscalar(e), c)), error(['rheome:bundle:' lower(name)], '%s: a struct array or a cell of structs is expected.', name); end
end

function need(x, keys, name)
    miss = keys(~isfield(x, keys));
    if ~isempty(miss), error(['rheome:bundle:' lower(name)], '%s: each needs %s.', name, strjoin(miss, ', ')); end
end

function put(folder, name, x, cls)
    f = fullfile(folder, name);
    fid = fopen(f, 'w', 'ieee-le');
    if fid < 0, error('rheome:bundle:write', 'Cannot write %s.', f); end
    fwrite(fid, x, cls); fclose(fid);
end
