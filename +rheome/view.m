function out = view(varargin)
% VIEW  Show a result on the cortex in Rheome's 3D viewer, from one call; rheome.show when there is no viewer.
%
%   rheome.view(surf, result)
%   rheome.view(result)                         % result.surface carries the surface
%   rheome.view(surf, struct('activity', I, 'velocity', v), Fs=100, Layers={{"activity"}, {"velocity", "x", "arrows"}})
%   rheome.view(surf, struct('J', J), Potentials="J", Layers={{"J_psi"}, {"J", "x", "arrows"}})   % J's stream function
%   out = rheome.view(..., Name=Value)
%
% OPTIONAL: nothing else in Rheome needs the viewer, and Rheome works without it (R2023b and later).
% The result is written as a rheome.bundle@1 folder (rheome.io.write.bundle) and opened:
%   R2026a+       inside MATLAB, a uifigure + uihtml on the viewer's viewer.html (the Cortical Flow viewer,
%                 nxr-cortical-flow editor/dist-viewer; its README.md has the format and protocol). The viewer folder is RHEOME_VIEWER, default prefdir/rheome-viewer, installed
%                 there from the copy shipped in <rheome>/viewer when that exists; the bundle is written beneath
%                 it, because uihtml serves only that folder's tree. The figure is reused by the next call.
%   R2023b-R2025b the Cortical Flow desktop app, by its URL (corticalflow://view?bundle=<folder>).
% With no viewer, no WebGPU (the viewer's 'no-webgpu' error) or no app that handles corticalflow://, the chosen
% field is drawn at the chosen time with rheome.show.surface (real) or rheome.show.vectors (vector).
%
% INPUTS:
%   surf   : surface struct (.Vertices meters, .Faces 1-based), a Brainstorm surface .mat path, or a struct of
%            surfaces of one mesh by name (struct('pial', S, 'inflated', Si, 'sphere', Ss); the first is the default)
%   result : [nV x T] or [3nV x T] array, or a struct of such arrays (see rheome.io.write.bundle); in the
%            one-argument form, a struct with a 'surface' field as surf above and the arrays
%
% OPTIONS (name/value) -- the viewer's settings; one not given is left to the viewer's default:
%   'Layers'       what is drawn, in draw order: each layer is ONE field, ONE form of it and ONE representation
%                  (the viewer's 'layers', its README 'Layers'). A cell of {field, form, as} cells (form and as
%                  optional), or a struct array (or cell of structs) with the viewer's layer keys:
%                    field      a field's name (required): one of the result's, or a potential (see Potentials)
%                    form       real field: "f" (default); vector field: "abs" |v|, "nx" n.v, "div", "curl" (scalar
%                               forms), "x" v, "xt" the tangential part (vector forms; default "x", or "abs" when
%                               'as' is a scalar representation)
%                    as         scalar form: "colour" (default) | "isolines" | "halfMax"; vector form: "arrows"
%                    visible    true | false;  colormap  a colormap name;  clim  [lo hi] or "auto";
%                    threshold  [min max], NaN an open side;  opacity  0-1;  bands  isolines, 1-100;
%                    density    arrows, (0, 1];  scale  arrows, >= 0 (m; 0: 4 % of the surface);
%                    normalize  arrows, true | false;  color  arrows, a neutral token ("paper", "ink" ...) or #rrggbb
%                  A vector field is never painted itself: paint a scalar form of it, or draw it as arrows.
%                  Not given: the first real field in colour and the first vector field as arrows.
%   'Potentials'   vector fields whose Hodge potentials are written too, as real fields <v>_phi and <v>_psi
%                  (rheome.io.write.bundle 'Potentials'): names, true (every vector field) or false (default)
%   'Colormap'     a colormap name (e.g. "viridis", "coolwarm", "magma", "turbo"): the colour layer's
%   'CLim'         [lo hi] colour limits, in the field's units (default: 1st-99th percentile): the colour layer's
%   'Threshold'    [min max]: paint only where min <= x <= max; NaN for an open side (e.g. [0.5 NaN]): the colour layer's
%   'Field'        DEPRECATED (one release; use Layers): a real field painted, a vector field as arrows
%   'Arrows'       DEPRECATED (one release; use Layers): true | false, the arrows of the first vector field, or none
%   'Surface'      the surface shown, a name of surf's struct (e.g. "inflated")
%   'Hemisphere'   "both" | "left" | "right" (the mesh's two largest connected components)
%   'View'         camera: "left" | "right" | "top" | "bottom" | "front" | "back"
%   'ArrowDensity' fraction of vertices with an arrow, (0, 1]: every arrows layer's density
%   'ArrowScale'   arrow length at the top of its range, in the surface's units (m); 0: 4 % of the surface's
%                  size: every arrows layer's scale
%   'Time'         the time shown, s from the first sample (the nearest frame)
%   'Play'         true | false: play on open
%   'Speed'        playback speed, x real time (> 0)
% and how it is opened:
%   'Fs'      sampling rate of the T columns, Hz (default 1)
%   'Name'    the bundle's name (default 'Rheome result')
%   'Viewer'  'auto' (default) | 'editor' (uihtml) | 'app' | 'show' | 'none' (write the bundle only)
%   'Timeout' seconds to wait for the viewer's 'ready' (default 60)
%
% OUTPUT:
%   out : struct .viewer (the path taken: 'editor' | 'app' | 'show' | 'none'), .manifest (the bundle's
%         manifest.json, '' for 'show'), .settings (the viewer's settings sent, as in the manifest's 'settings';
%         the layer options are sent as 'layers', a cell of layer structs, set whole),
%         .figure (the uifigure or figure, [] for 'app' and 'none'), and on the 'editor' path:
%         .ready      what the viewer said when it opened the bundle (its fields, surfaces, every setting)
%         .set(Name=Value, ...)      changes the open view (the options above) without reopening it; a layer
%                                    option given acts on the layers last sent
%         .screenshot(file, scale)   writes the view as a PNG (scale, default 1: the window's pixels)
%
% The viewer gets the settings twice: as the manifest's 'settings', and with the bundle in uihtml's Data
% (struct('bundle', <manifest, relative>, 'settings', settings)); it answers with the event 'ready' or 'error'
% (data.code 'no-webgpu' when it cannot draw: rheome.show instead; 'settings' when it refused one: a warning,
% the rest applied; 'bundle' or 'boot': an error).
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', 0.07*V, 'Faces', F);
%   I = exp(-sum((V - [0 0 1]).^2, 2) / 0.1) .* cos(linspace(0, 2*pi, 20));
%   out = rheome.view(S, I, Fs=10, CLim=[-1 1], View="top", Viewer="show");
%   assert(strcmp(out.viewer, 'show'));  close(out.figure)
%
% See also: rheome.io.write.bundle, rheome.show.surface, rheome.show.vectors
%
% Author: Diellor Basha, 2026

    [surf, result, args] = splitInputs(varargin);
    o = options(args, true);
    if ischar(surf) || isstring(surf), surf = rheome.io.read.surface(surf); end
    if isnumeric(result), result = struct('result', result); end
    T = size(result.(firstField(result)), 2);
    s = settings(o, surf, result, T, []);
    out = struct('viewer', lower(char(o.Viewer)), 'manifest', '', 'settings', s, 'figure', []);

    vdir = viewerFolder();
    page = fullfile(vdir, 'viewer.html');
    if strcmp(out.viewer, 'auto')
        if ~isMATLABReleaseOlderThan('R2026a') && isfile(page), out.viewer = 'editor'; else, out.viewer = 'app'; end
    end

    if ~strcmp(out.viewer, 'show')
        % one bundle at a time: the last one is open and loaded, older ones are dropped (new name: no stale cache)
        root = fullfile(vdir, 'bundles');
        if isfolder(root), rmdir(root, 's'); end
        id = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'));
        out.manifest = rheome.io.write.bundle(fullfile(root, id), surf, result, 'Fs', o.Fs, 'Name', o.Name, 'Settings', s, ...
            'Potentials', o.Potentials);
    end

    switch out.viewer
        case 'editor'
            if ~isfile(page), error('rheome:view:noviewer', 'No viewer at %s (set RHEOME_VIEWER).', page); end
            [ok, out.figure, h, out.ready] = viaEditor(page, ['bundles/' id '/manifest.json'], s, o.Timeout);
            if ok && ~isempty(h)
                setappdata(h, 'rheomeLayers', []);
                if isfield(s, 'layers'), setappdata(h, 'rheomeLayers', s.layers); end
                out.set = @(varargin) setView(h, setfield(options(varargin, false), 'Fs', o.Fs), surf, result, T); %#ok<SFLD>
                out.screenshot = @(file, varargin) screenshot(h, file, o.Timeout, varargin{:});
            elseif ~ok
                out.viewer = 'show';
            end
        case 'app'
            if ~viaApp(fileparts(out.manifest)), out.viewer = 'show'; end
    end
    if strcmp(out.viewer, 'show'), out.figure = viaShow(surf, result, s); end
end

function [surf, result, args] = splitInputs(in)
% rheome.view(surf, result, ...) or rheome.view(result, ...) with result.surface
    if isempty(in), error('rheome:view:input', 'A result is expected.'); end
    a = in{1};
    if numel(in) >= 2 && (ischar(a) || isstring(a) || (isstruct(a) && ~isfield(a, 'surface')))
        surf = a; result = in{2}; args = in(3:end);
    elseif isstruct(a) && isfield(a, 'surface')
        surf = a.surface; result = rmfield(a, 'surface'); args = in(2:end);
    else
        result = {};
    end
    if ~(isnumeric(result) || isstruct(result))
        error('rheome:view:input', 'Call rheome.view(surf, result, ...) or rheome.view(result, ...) with result.surface.');
    end
end

function o = options(args, opening)
% the name/value options; opening: also how it is opened (a .set call takes the settings only)
    p = inputParser;
    p.FunctionName = 'rheome.view';
    for k = {'Layers', 'Field', 'Arrows', 'Colormap', 'CLim', 'Threshold', 'Surface', 'Hemisphere', 'View', 'ArrowDensity', ...
            'ArrowScale', 'Time', 'Play', 'Speed'}
        p.addParameter(k{1}, []);
    end
    if opening
        p.addParameter('Fs', 1);
        p.addParameter('Potentials', false);
        p.addParameter('Name', 'Rheome result');
        p.addParameter('Viewer', 'auto', @(x) any(strcmpi(x, {'auto', 'editor', 'app', 'show', 'none'})));
        p.addParameter('Timeout', 60);
    end
    p.parse(args{:});
    o = p.Results;
    if ~opening, o.Potentials = false; end
end

function f = firstField(result)
    f = fieldnames(result);
    if isempty(f), error('rheome:view:input', 'The result has no field.'); end
    f = f{1};
end

function s = settings(o, surf, result, T, base)
% the options given, checked and named as the viewer's settings (README.md 'The settings'); none given: its defaults.
% base: the layers last sent (a .set call), {} or [] on opening
    s = struct();
    S1 = surf;
    if ~isfield(S1, 'Vertices'), S1 = surf.(firstField(surf)); end
    [names, types] = fieldTypes(result, size(S1.Vertices, 1), o.Potentials);
    L = layers(o, names, types, base);
    if ~isempty(L), s.layers = L; end
    if ~isempty(o.Surface)
        if isfield(surf, 'Vertices'), names = {'default'}; else, names = fieldnames(surf); end
        s.surface = pick(o.Surface, names, 'Surface');
    end
    if ~isempty(o.Hemisphere), s.hemisphere = pick(o.Hemisphere, {'both', 'left', 'right'}, 'Hemisphere'); end
    if ~isempty(o.View), s.view = pick(o.View, {'left', 'right', 'top', 'bottom', 'front', 'back'}, 'View'); end
    if ~isempty(o.Time)
        t = o.Time;
        if ~isscalar(t) || ~isfinite(t) || t < 0 || t > (T - 1) / o.Fs, bad('Time', sprintf('in [0, %g] s', (T - 1) / o.Fs)); end
        s.frame = min(T - 1, round(t * o.Fs));
    end
    if ~isempty(o.Play), s.playing = flag(o.Play, 'Play'); end
    if ~isempty(o.Speed), s.fps = positive(o.Speed, 'Speed', false) * o.Fs; end
end

function [names, types] = fieldTypes(result, nV, pot)
% the fields a layer may name and their types: the result's, then the potentials the bundle adds (<v>_phi, <v>_psi)
    names = fieldnames(result)';
    types = repmat({'real'}, size(names));
    types(cellfun(@(f) size(result.(f), 1) == 3*nV, names)) = {'vector'};
    vec = names(strcmp(types, 'vector'));
    if (islogical(pot) || isnumeric(pot)) && isscalar(pot)
        if ~pot, vec = {}; end
    else
        vec = cellstr(pot);   % checked by rheome.io.write.bundle
    end
    for v = vec(:)'
        names = [names {[v{1} '_phi'], [v{1} '_psi']}]; %#ok<AGROW>
        types = [types {'real', 'real'}]; %#ok<AGROW>
    end
end

function L = layers(o, names, types, base)
% the layers to send (a cell of structs, set whole), or [] when no layer option is given
    L = [];
    if all(cellfun(@isempty, {o.Layers, o.Field, o.Arrows, o.Colormap, o.CLim, o.Threshold, o.ArrowDensity, o.ArrowScale}))
        return
    end
    if ~isempty(o.Layers)
        L = o.Layers;
        if isstruct(L), L = num2cell(L); end
        if ~iscell(L), bad('Layers', 'a cell of {field, form, as} cells or a struct array of layers'); end
        if ~isempty(L) && ~iscell(L{1}) && ~isstruct(L{1}), L = {L}; end   % one layer {field, form, as}
        L = cellfun(@(x) layer(x, names, types), L(:)', 'UniformOutput', false);
    elseif ~isempty(base)
        L = base;
    else   % the viewer's defaults: the first real field in colour, the first vector field as arrows
        L = {};
        r = find(strcmp(types, 'real'), 1);  v = find(strcmp(types, 'vector'), 1);
        if ~isempty(r), L{end+1} = struct('field', names{r}, 'form', 'f', 'as', 'colour'); end
        if ~isempty(v), L{end+1} = struct('field', names{v}, 'form', 'x', 'as', 'arrows'); end
    end
    scal = @(L) cellfun(@(x) ~strcmp(x.as, 'arrows'), L);
    if ~isempty(o.Field)   % the old 'scalar' / 'vector': a real field in the scalar layers, a vector one as the arrows
        warning('rheome:view:deprecated', 'Field is deprecated (removed in the next release); use Layers, e.g. Layers={{"%s"}}.', char(o.Field));
        f = pick(o.Field, names, 'Field');
        if strcmp(types{strcmp(names, f)}, 'vector')
            k = find(~scal(L), 1);
            if isempty(k), L{end+1} = struct('field', f, 'form', 'x', 'as', 'arrows'); else, L{k}.field = f; end
        elseif ~any(scal(L))
            L = [{struct('field', f, 'form', 'f', 'as', 'colour')} L];
        else
            for k = find(scal(L)), L{k}.field = f; L{k}.form = 'f'; end
        end
    end
    if ~isempty(o.Arrows)
        warning('rheome:view:deprecated', 'Arrows is deprecated (removed in the next release); use Layers, a layer {field, "x", "arrows"}.');
        if ~flag(o.Arrows, 'Arrows'), L = L(scal(L));
        elseif all(scal(L))
            v = find(strcmp(types, 'vector'), 1);
            if isempty(v), bad('Arrows', 'a vector field ([3nV x T]) in the result'); end
            L{end+1} = struct('field', names{v}, 'form', 'x', 'as', 'arrows');
        end
    end
    as = @(r) find(cellfun(@(x) any(strcmp(x.as, r)), L));
    L = look(L, as({'colour'}), 'colour', {'colormap', o.Colormap, 'Colormap'; 'threshold', o.Threshold, 'Threshold'}, names, types);
    L = look(L, as({'colour', 'isolines'}), 'colour or isolines', {'clim', o.CLim, 'CLim'}, names, types);
    L = look(L, find(~scal(L)), 'arrows', {'density', o.ArrowDensity, 'ArrowDensity'; 'scale', o.ArrowScale, 'ArrowScale'}, names, types);
end

function L = look(L, ks, what, kv, names, types)
% the look options given (key, value, option name), set on the layers ks
    for i = 1:size(kv, 1)
        if isempty(kv{i, 2}), continue; end
        if isempty(ks), bad(kv{i, 3}, sprintf('a %s layer to act on', what)); end
        for k = ks
            x = L{k};  x.(kv{i, 1}) = kv{i, 2};
            L{k} = layer(x, names, types);
        end
    end
end

function y = layer(x, names, types)
% one layer, checked: {field, form, as} or a struct with the viewer's layer keys (an empty one: not given)
    if iscell(x)
        k = {'field', 'form', 'as'};
        if isempty(x) || numel(x) > 3, bad('Layers', 'a layer {field, form, as}'); end
        x = cell2struct(x(:), k(1:numel(x)), 1);
    end
    if ~isstruct(x) || ~isscalar(x) || ~isfield(x, 'field'), bad('Layers', 'a layer with a field'); end
    keys = fieldnames(x)';
    other = setdiff(keys, {'field', 'form', 'as', 'visible', 'colormap', 'clim', 'threshold', 'opacity', 'bands', ...
        'density', 'scale', 'normalize', 'color'}, 'stable');
    if ~isempty(other), bad('Layers', ['the viewer''s layer keys, not ' strjoin(other, ', ')]); end
    f = pick(x.field, names, 'Layers field');
    scalarAs = {'colour', 'isolines', 'halfMax'};
    if strcmp(types{strcmp(names, f)}, 'vector'), forms = {'abs', 'nx', 'div', 'curl', 'x', 'xt'}; else, forms = {'f'}; end
    as = '';  if isfield(x, 'as') && ~isempty(x.as), as = pick(x.as, [scalarAs {'arrows'}], 'Layers as'); end
    if isfield(x, 'form') && ~isempty(x.form), form = pick(x.form, forms, ['Layers form of ' f]);
    elseif isscalar(forms), form = 'f';
    elseif isempty(as) || strcmp(as, 'arrows'), form = 'x';
    else, form = 'abs';
    end
    vecForm = any(strcmp(form, {'x', 'xt'}));
    if isempty(as), as = scalarAs{1}; if vecForm, as = 'arrows'; end, end
    if vecForm ~= strcmp(as, 'arrows')
        bad('Layers', sprintf('"%s" of %s drawn as %s: a vector form as arrows, a scalar form as colour, isolines or halfMax', form, f, as));
    end
    y = struct('field', f, 'form', form, 'as', as);
    for k = setdiff(keys, {'field', 'form', 'as'}, 'stable')
        v = x.(k{1});
        if isempty(v), continue; end   % a struct array's unset key: the viewer's default
        switch k{1}
            case {'visible', 'normalize'}, v = flag(v, k{1});
            case {'colormap', 'color'}
                if ~(ischar(v) || (isstring(v) && isscalar(v))), bad(k{1}, 'a name'); end
                v = char(v);
            case 'clim'
                if (ischar(v) || isstring(v)) && strcmpi(v, 'auto'), v = 'auto';
                else
                    v = double(v);
                    if numel(v) ~= 2 || ~all(isfinite(v)) || v(1) >= v(2), bad('clim', '[lo hi], lo < hi, or "auto"'); end
                    v = v(:)';
                end
            case 'threshold'
                if isstruct(v) && all(isfield(v, {'min', 'max'})), t = [v.min v.max]; else, t = double(v); end   % checked already
                if numel(t) ~= 2 || any(isinf(t)) || t(1) > t(2), bad('threshold', '[min max], min <= max, NaN for an open side'); end
                v = struct('min', t(1), 'max', t(2));        % NaN is JSON null: that side open
            case 'opacity', v = positive(v, 'opacity', true);  if v > 1, bad('opacity', 'in [0, 1]'); end
            case 'bands', v = positive(v, 'bands', false);  if v > 100 || v ~= round(v), bad('bands', 'an integer in 1-100'); end
            case 'density', v = positive(v, 'density', false);  if v > 1, bad('density', 'in (0, 1]'); end
            case 'scale', v = positive(v, 'scale', true);
        end
        y.(k{1}) = v;
    end
end

function setView(h, o, surf, result, T)
% .set: the settings given; a layer option acts on the layers last sent
    s = settings(o, surf, result, T, getappdata(h, 'rheomeLayers'));
    if isfield(s, 'layers'), setappdata(h, 'rheomeLayers', s.layers); end
    sendEventToHTMLSource(h, 'set', s);
end

function x = pick(x, allowed, name)
    if ~(ischar(x) || (isstring(x) && isscalar(x))), bad(name, ['one of: ' strjoin(allowed, ', ')]); end
    k = find(strcmpi(char(x), allowed), 1);
    if isempty(k), bad(name, ['one of: ' strjoin(allowed, ', ')]); end
    x = allowed{k};
end

function x = positive(x, name, zeroOk)
    if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x < 0 || (~zeroOk && x == 0), bad(name, 'a positive number'); end
    x = double(x);
end

function x = flag(x, name)
    if ~isscalar(x) || ~(islogical(x) || isnumeric(x)), bad(name, 'true or false'); end
    x = logical(x);
end

function bad(name, what)
    error('rheome:view:option', '%s: %s is expected.', name, what);
end

function [ok, fig, h, ready] = viaEditor(page, manifest, s, timeout)
% the viewer inside MATLAB; ok = false when it has no WebGPU; h the uihtml ([] when it did not answer)
    persistent hh
    if isempty(hh) || ~isvalid(hh) || ~strcmp(hh.HTMLSource, page)
        fig = uifigure('Name', 'Rheome - viewer', 'Position', [100 100 1400 900]);
        g = uigridlayout(fig, [1 1], 'Padding', 0);
        hh = uihtml(g, 'HTMLSource', page, 'HTMLEventReceivedFcn', @onEvent);
    end
    h = hh;  fig = ancestor(h, 'figure');  ready = [];
    figure(fig);
    h.UserData = struct();
    h.Data = struct('bundle', manifest, 'settings', s);   % a reused figure: DataChanged opens it in place of the last
    [name, d] = waitFor(h, {'ready', 'error'}, timeout);
    ok = true;
    switch name
        case 'ready', ready = d;
        case 'error'
            if isstruct(d) && isfield(d, 'code') && strcmp(d.code, 'no-webgpu')
                delete(fig); fig = []; h = []; ok = false;
                return
            end
            error('rheome:view:editor', 'The viewer could not open the bundle (%s): %s', field(d, 'code'), field(d, 'message'));
        otherwise
            if isvalid(h), warning('rheome:view:timeout', 'The viewer did not say ready within %g s; the figure is left open.', timeout); end
            h = [];
    end
end

function [name, d] = waitFor(h, names, timeout)
% the first of the named events the viewer sends (onEvent keeps each in UserData), '' at the timeout
    name = '';  d = [];  t = tic;
    while isvalid(h) && toc(t) < timeout
        k = find(isfield(h.UserData, names), 1);
        if ~isempty(k), name = names{k}; d = h.UserData.(name); return; end
        pause(0.1);
    end
end

function onEvent(src, ev)
% each event kept by name; a refused setting is a warning (the viewer applied the rest), never the open's answer
    d = ev.HTMLEventData;
    if strcmp(ev.HTMLEventName, 'error') && isstruct(d) && isfield(d, 'code') && any(strcmp(d.code, {'settings', 'screenshot'}))
        if strcmp(d.code, 'settings'), warning('rheome:view:settings', 'The viewer refused a setting: %s', field(d, 'message')); return; end
        src.UserData.screenshot = d;  return      % a failed screenshot answers the screenshot call
    end
    src.UserData.(ev.HTMLEventName) = d;
end

function x = field(d, name)
    x = '';
    if isstruct(d) && isfield(d, name), x = char(string(d.(name))); end
end

function screenshot(h, file, timeout, scale)
% the view as a PNG: the viewer answers 'screenshot' with a data URL
    if nargin < 4, scale = 1; end
    if ~isvalid(h), error('rheome:view:closed', 'The viewer window is closed.'); end
    if isfield(h.UserData, 'screenshot'), h.UserData = rmfield(h.UserData, 'screenshot'); end
    sendEventToHTMLSource(h, 'screenshot', struct('scale', scale));
    [name, d] = waitFor(h, {'screenshot'}, timeout);
    png = field(d, 'png');
    if isempty(name) || ~startsWith(png, 'data:image/png;base64,')
        error('rheome:view:screenshot', 'No screenshot from the viewer: %s', field(d, 'message'));
    end
    fid = fopen(file, 'w');
    if fid < 0, error('rheome:view:screenshot', 'Cannot write %s.', file); end
    fwrite(fid, matlab.net.base64decode(extractAfter(png, ',')));  fclose(fid);
end

function ok = viaApp(folder)
% the desktop app by its URL; ok = false when no app handles corticalflow://
    url = ['corticalflow://view?bundle=' urlencode(folder)];
    if ismac
        cmd = sprintf('open "%s"', url);
    elseif ispc
        [s, ~] = system('reg query HKCR\corticalflow');
        if s ~= 0, ok = false; return; end
        cmd = sprintf('start "" "%s"', url);
    else
        cmd = sprintf('xdg-open "%s"', url);
    end
    [s, ~] = system(cmd);
    ok = s == 0;
end

function fig = viaShow(surf, result, s)
% no viewer: the chosen field at the chosen time in a MATLAB figure, with the settings rheome.show has
    if ~isfield(surf, 'Vertices')
        names = fieldnames(surf); k = 1;
        if isfield(s, 'surface'), k = find(strcmp(names, s.surface)); end
        S1 = surf.(names{1}); surf = surf.(names{k}); surf.Faces = S1.Faces;
    end
    % the first visible colour layer, else the first visible arrows layer.
    % shortcut: a vector field's scalar form is drawn as its vectors, a potential as its vector field (show has no forms)
    f = firstField(result);  c = struct();
    if isfield(s, 'layers')
        vis = s.layers(cellfun(@(x) ~isfield(x, 'visible') || x.visible, s.layers));
        k = find(cellfun(@(x) strcmp(x.as, 'colour'), vis), 1);
        if isempty(k), k = find(cellfun(@(x) strcmp(x.as, 'arrows'), vis), 1); end
        if ~isempty(k), c = vis{k}; f = c.field; end
    end
    if ~isfield(result, f), f = extractBefore(f, max(strfind(f, '_'))); end   % <v>_phi, <v>_psi: v
    x = result.(f);
    t = 1;
    if isfield(s, 'frame'), t = s.frame + 1; end
    opt = {'Title', sprintf('%s, frame %d of %d (rheome.show)', f, t, size(x, 2))};
    if isfield(c, 'clim') && isnumeric(c.clim), opt = [opt {'CLim', c.clim}]; end
    if isfield(s, 'view')   % SCS: +x front, +y left, +z top; MATLAB's az 0 looks from -y
        opt = [opt {'View', struct('left', [180 0], 'right', [0 0], 'top', [0 90], 'bottom', [0 -90], ...
            'front', [90 0], 'back', [-90 0]).(s.view)}];
    end
    if isfield(c, 'colormap')
        if exist(c.colormap, 'file') == 2, opt = [opt {'Colormap', feval(c.colormap, 256)}]; end   % a MATLAB colormap name; else the default
    end
    if size(x, 1) == size(surf.Vertices, 1)
        ax = rheome.show.surface(surf, x(:, t), opt{:});
    else
        ax = rheome.show.vectors(surf, x(:, t), opt{:});
    end
    fig = ancestor(ax, 'figure');
end
