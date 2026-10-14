function manifest = bundle(folder, surf, result, varargin)
% IO.WRITE.BUNDLE  Write a result on a surface as a rheome.bundle@1 folder (what rheome.view hands the viewer).
%
%   manifest = rheome.io.write.bundle(folder, surf, result)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Fs', fs, 'Name', name, 'Settings', s)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Potentials', "J")   % + J's Hodge potentials phi, psi
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Objects', "meg-helmet-306")  % + a scene object
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
%   'Objects'   SCENE OBJECTS for illustrations (a helmet, a scalp, electrodes, a scanner), written as the manifest's
%               'objects' (static meshes the viewer draws around the cortex; they carry no field). A struct array, or
%               a string / cell array of them, each one object:
%                 Name       its name (letters, digits, - and _), the viewer's 'objects' setting names it      REQUIRED
%                 Vertices   [n x 3] in Unit, and Faces [m x 3] 1-based (as a surface struct)          OR Folder
%                 Folder     an object folder (positions.bin, faces.bin, object.json), or the name of one in the
%                            media library (~/workspace/library/media/illustrations/objects/<name>); a string
%                            alone is a Folder. Its object.json gives Name, Unit, Transform, Role and Color, which
%                            the struct's fields override
%                 Unit       "m" | "cm" | "mm" (default "m"; the cortex is in m)
%                 Transform  4 x 4 affine object -> the surface's frame, in the surface's unit (default identity).
%                            A media object's default is its to_scs (fsaverage RAS -> SCS: a rotation only; an SCS
%                            cortex whose origin is not fsaverage's needs its own Transform)
%                 Role       "skin" | "anatomy" | "device" | "device2" | "line" (default the viewer's: "device")
%                 Color      "#rrggbb" or a neutral token name, in place of the role's colour
%                 Visible    shown when the viewer opens (default true)
%               Each is written as <Name>_p.bin (float32 [n x 3]) and <Name>_f.bin (int32 [m x 3], 0-based).
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
    p.addParameter('Objects', []);
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
    obj = objects(o.Objects);

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
    if ~isempty(obj), m.objects = writeObjects(folder, obj); end
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

function obj = objects(x)
% the scene objects, checked and completed (a media folder read), before anything is written
    if isempty(x), obj = {}; return; end
    if ischar(x) || isstring(x), x = cellstr(x); end
    if isstruct(x), x = num2cell(x); end
    if ~iscell(x), error('rheome:bundle:objects', 'Objects: a struct array, or a string / cell array of them.'); end
    obj = cell(1, numel(x));
    for k = 1:numel(x)
        a = x{k};
        if ischar(a) || isstring(a), a = struct('Folder', char(a)); end
        if ~isstruct(a) || ~isscalar(a), error('rheome:bundle:objects', 'Objects{%d}: a struct or a folder name.', k); end
        bad = setdiff(fieldnames(a), {'Name', 'Vertices', 'Faces', 'Folder', 'Unit', 'Transform', 'Role', 'Color', 'Visible'});
        if ~isempty(bad), error('rheome:bundle:objects', 'Objects{%d}: %s is not an object field.', k, strjoin(bad, ', ')); end
        e = struct('name', '', 'unit', 'm');
        if isfield(a, 'Folder')
            [a, e] = fromFolder(a, k);
        else
            if ~all(isfield(a, {'Vertices', 'Faces'})), error('rheome:bundle:objects', 'Objects{%d}: Vertices and Faces, or a Folder.', k); end
            V = a.Vertices;  F = a.Faces;
            if ~isnumeric(V) || size(V, 2) ~= 3 || isempty(V) || ~all(isfinite(V(:))), error('rheome:bundle:objects', 'Objects{%d}: Vertices [n x 3], finite.', k); end
            if ~isnumeric(F) || size(F, 2) ~= 3 || isempty(F) || any(F(:) < 1 | F(:) > size(V, 1) | F(:) ~= round(F(:)))
                error('rheome:bundle:objects', 'Objects{%d}: Faces [m x 3], 1-based indices into its Vertices.', k);
            end
            a.data = {V', 'single'; F' - 1, 'int32'};
        end
        if isfield(a, 'Name'), e.name = char(a.Name); end
        if isempty(regexp(e.name, '^[A-Za-z0-9_-]+$', 'once')), error('rheome:bundle:objects', 'Objects{%d}: Name "%s": letters, digits, - and _.', k, e.name); end
        if isfield(a, 'Unit'), e.unit = char(a.Unit); end
        if ~any(strcmp(e.unit, {'m', 'cm', 'mm'})), error('rheome:bundle:objects', '%s: Unit m, cm or mm.', e.name); end
        if isfield(a, 'Transform')
            M = double(a.Transform);
            if isequal(size(M), [4 4]), M = reshape(M', 1, []); end   % row-major, as the viewer reads it
            if numel(M) ~= 16 || ~all(isfinite(M(:))), error('rheome:bundle:objects', '%s: Transform a 4 x 4 affine.', e.name); end
            e.transform = M(:)';
        end
        if isfield(a, 'Role')
            e.role = char(a.Role);
            if ~any(strcmp(e.role, {'skin', 'anatomy', 'device', 'device2', 'line'})), error('rheome:bundle:objects', '%s: Role skin, anatomy, device, device2 or line.', e.name); end
        end
        if isfield(a, 'Color'), e.color = char(a.Color); end
        if isfield(a, 'Visible'), e.visible = logical(a.Visible); end
        obj{k} = struct('entry', e, 'data', {a.data});
    end
    n = cellfun(@(o) o.entry.name, obj, 'UniformOutput', false);
    if numel(unique(n)) < numel(n), error('rheome:bundle:objects', 'Objects: each Name once.'); end
end

function [a, e] = fromFolder(a, k)
% an object folder (the media library's, by name): its object.json's manifest entry and its two files, checked
    d = char(a.Folder);
    if ~isfolder(d), d = fullfile(char(java.lang.System.getProperty('user.home')), 'workspace', 'library', 'media', 'illustrations', 'objects', d); end
    if ~isfile(fullfile(d, 'object.json')), error('rheome:bundle:objects', 'Objects{%d}: no object folder %s.', k, char(a.Folder)); end
    j = jsondecode(fileread(fullfile(d, 'object.json')));
    me = struct();  if isfield(j, 'manifest_entry'), me = j.manifest_entry; end
    e = struct('name', j.name, 'unit', 'm');
    if isfield(j, 'unit'), e.unit = j.unit; end
    if isfield(me, 'transform') && ~isempty(me.transform), a = withDefault(a, 'Transform', me.transform); end
    for f = {'role', 'Role'; 'color', 'Color'}'
        if isfield(me, f{1}), a = withDefault(a, f{2}, me.(f{1})); elseif isfield(j, f{1}), a = withDefault(a, f{2}, j.(f{1})); end
    end
    P = rd(fullfile(d, 'positions.bin'), 'single');  Fc = rd(fullfile(d, 'faces.bin'), 'int32');
    if mod(numel(P), 3) || mod(numel(Fc), 3) || isempty(Fc) || any(Fc < 0 | Fc >= numel(P) / 3) || ~all(isfinite(P))
        error('rheome:bundle:objects', '%s: positions.bin / faces.bin are not float32 [n x 3] and 0-based int32 [m x 3] into them.', d);
    end
    a.data = {P, 'single'; Fc, 'int32'};
end

function a = withDefault(a, f, v)
    if ~isfield(a, f), a.(f) = v; end
end

function x = rd(f, cls)
    fid = fopen(f, 'r', 'ieee-le');
    if fid < 0, error('rheome:bundle:objects', 'Cannot read %s.', f); end
    x = fread(fid, inf, ['*' cls]);  fclose(fid);
end

function list = writeObjects(folder, obj)
% each object's two files, <name>_p.bin and <name>_f.bin, and its manifest entry (a JSON array even for one)
    list = cell(1, numel(obj));
    for k = 1:numel(obj)
        e = obj{k}.entry;
        e.positions = [e.name '_p.bin'];  e.faces = [e.name '_f.bin'];
        put(folder, e.positions, obj{k}.data{1, 1}, obj{k}.data{1, 2});
        put(folder, e.faces, obj{k}.data{2, 1}, obj{k}.data{2, 2});
        list{k} = e;
    end
end

function put(folder, name, x, cls)
    f = fullfile(folder, name);
    fid = fopen(f, 'w', 'ieee-le');
    if fid < 0, error('rheome:bundle:write', 'Cannot write %s.', f); end
    fwrite(fid, x, cls); fclose(fid);
end
