function manifest = bundle(folder, surf, result, varargin)
% IO.WRITE.BUNDLE  Write a result on a surface as a rheome.bundle@1 folder (what rheome.view hands the viewer).
%
%   manifest = rheome.io.write.bundle(folder, surf, result)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Fs', fs, 'Name', name, 'Settings', s)
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Potentials', "J")   % + J's Hodge potentials phi, psi
%   manifest = rheome.io.write.bundle(folder, surf, result, 'Tracks', tracks)     % + time-tracks under the stage
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
%   'Tracks'  TIME-TRACKS, the tracks rheome.show.timetracks takes (a struct array or cell of structs: Label, Unit
%             REQUIRED; Kind, Spec (the plant, rheome.show.track's), Values (detected; NaN none), Time (s, default the
%             bundle's time axis), Driver), written as the manifest's 'tracks' with the viewer's lowercase keys. The
%             viewer evaluates each plant over the bundle's time axis with the Editor's trackValue (default: none)
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
    p.addParameter('Tracks', []);
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
    tr = tracks(o.Tracks);

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
    if ~isempty(tr), m.tracks = tr; end
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

function tr = tracks(in)
% the time-tracks as the viewer reads them (tracks.ts): lowercase keys, arrays always arrays (a cell: jsonencode writes a
% 1-element vector or a 1x1 struct bare), NaN as null; each plant checked by rheome.show.track, as timetracks does
    tr = {};
    if isempty(in), return; end
    if isstruct(in), in = num2cell(in); end
    if ~iscell(in), error('rheome:bundle:tracks', 'Tracks: a struct array or a cell of structs (rheome.show.timetracks).'); end
    for i = 1:numel(in)
        a = in{i};  at = sprintf('Tracks(%d)', i);
        if ~isstruct(a) || ~isscalar(a), error('rheome:bundle:tracks', '%s: a struct.', at); end
        bad = setdiff(fieldnames(a), {'Label', 'Unit', 'Kind', 'Spec', 'Values', 'Time', 'Driver'});
        if ~isempty(bad), error('rheome:bundle:tracks', '%s: %s is not a track field.', at, strjoin(bad, ', ')); end
        if ~isfield(a, 'Label') || isempty(a.Label) || ~isfield(a, 'Unit') || isempty(a.Unit)
            error('rheome:bundle:tracks', '%s: Label and Unit are required (a lane says its quantity and unit).', at);
        end
        t = struct('label', char(a.Label), 'unit', char(a.Unit));
        if isfield(a, 'Kind') && ~isempty(a.Kind)
            t.kind = char(a.Kind);
            if ~any(strcmp(t.kind, {'value', 'angle', 'sign'})), error('rheome:bundle:tracks', '%s.Kind: value, angle or sign.', at); end
        end
        has = @(f) isfield(a, f) && ~isempty(a.(f));
        if has('Spec'), rheome.show.track(a.Spec, 0);  t.spec = spec(a.Spec); end   % refuses a bad plant here
        if has('Values')
            t.values = num2cell(double(a.Values(:)'));
            if has('Time')
                if numel(a.Time) ~= numel(a.Values), error('rheome:bundle:tracks', '%s: %d values and %d times.', at, numel(a.Values), numel(a.Time)); end
                t.time = num2cell(double(a.Time(:)'));
            end
        end
        if ~isfield(t, 'spec') && ~isfield(t, 'values'), error('rheome:bundle:tracks', '%s: a Spec (the plant) or Values (detected).', at); end
        if has('Driver'), t.driver = char(a.Driver); end
        tr{end+1} = t; %#ok<AGROW>
    end
end

function s = spec(a)
% rheome.show.track's spec with the viewer's keys: keys [{x, v, seg}], start, to, profile, osc {rate, depth, phase, shape}
    s = struct();
    if isfield(a, 'Keys') && ~isempty(a.Keys)
        K = a.Keys;
        if isnumeric(K), K = struct('x', num2cell(K(:, 1)), 'v', num2cell(K(:, 2))); end
        k = cell(1, numel(K));
        for j = 1:numel(K)
            k{j} = struct('x', K(j).x, 'v', K(j).v);
            if isfield(K, 'seg') && ~isempty(K(j).seg), k{j}.seg = char(K(j).seg); end
        end
        s.keys = k;
    end
    for f = {'Start', 'To', 'Profile'}
        if isfield(a, f{1}) && ~isempty(a.(f{1})), s.(lower(f{1})) = a.(f{1}); end
    end
    if isfield(s, 'profile'), s.profile = char(s.profile); end
    if isfield(a, 'Osc') && ~isempty(a.Osc)
        o = a.Osc;  s.osc = struct('rate', o.rate, 'depth', o.depth);
        if isfield(o, 'phase') && ~isempty(o.phase), s.osc.phase = o.phase; end
        if isfield(o, 'shape') && ~isempty(o.shape), s.osc.shape = char(o.shape); end
    end
end

function put(folder, name, x, cls)
    f = fullfile(folder, name);
    fid = fopen(f, 'w', 'ieee-le');
    if fid < 0, error('rheome:bundle:write', 'Cannot write %s.', f); end
    fwrite(fid, x, cls); fclose(fid);
end
