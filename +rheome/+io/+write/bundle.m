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
%   'Labellings'  struct array (or cell of structs), one partition each: .name, .codes [nV x 1] integers (< 0: none),
%               optional .entries struct array (.code .name .color '#rrggbb'). rheome.viz.tiles makes one from a
%               rheome.geom.tree level. The viewer paints it in the colour layer's place (setting 'labelling').
%   'Points'    struct array (or cell of structs), one point set each: .name, .vertices [n x 1] (1-BASED), optional
%               .frames [n x 1] (1-based frame each point is shown at; 0: every frame; none: every frame), .signs
%               [n x 1] (+1 / -1, drawn as plus / minus; none: +1), .color '#rrggbb'. rheome.viz.points makes one
%               from rheome.detect.criticalPoints / rheome.flow.critical.
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
