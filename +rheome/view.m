function out = view(varargin)
% VIEW  Show a result on the cortex in Rheome's 3D viewer, from one call; rheome.show when there is no viewer.
%
%   rheome.view(surf, result)
%   rheome.view(result)                         % result.surface carries the surface
%   rheome.view(surf, struct('activity', I, 'velocity', v), Fs=100, Field="activity", CLim=[0 1], Play=true)
%   out = rheome.view(..., Name=Value)
%
% OPTIONAL: nothing else in Rheome needs the viewer, and Rheome works without it (R2023b and later).
% The result is written as a rheome.bundle@1 folder (rheome.io.write.bundle) and opened:
%   R2026a+       inside MATLAB, a uifigure + uihtml on the viewer's matlab.html (the Cortical Flow viewer's
%                 MATLAB build). The viewer folder is RHEOME_VIEWER, default prefdir/rheome-viewer, installed
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
%   'Field'        the field shown (default: the first)
%   'Colormap'     a colormap name (e.g. "viridis", "coolwarm") or an [n x 3] RGB matrix in [0, 1]
%   'CLim'         [lo hi] colour limits, in the field's units
%   'Threshold'    values with |x| below it are not coloured (>= 0)
%   'Surface'      the surface shown, a name of surf's struct (e.g. "inflated")
%   'Hemisphere'   "both" | "left" | "right" (needs the surface's Structures split, .Hemi)
%   'View'         camera: "lateral" | "medial" | "dorsal" | "ventral" | "anterior" | "posterior" | [az el] degrees
%   'Arrows'       true | false: the vector field as arrows
%   'ArrowDensity' fraction of vertices with an arrow, (0, 1]
%   'ArrowScale'   arrow length factor (> 0)
%   'Overlays'     names of further fields drawn over the shown one (e.g. "velocity", "critical")
%   'Time'         the time shown, s from the first sample
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
%         manifest.json, '' for 'show'), .settings (the settings sent, as in the manifest's 'view'), .figure (the
%         uifigure or figure, [] for 'app' and 'none').
%
% The viewer gets the settings twice: as the manifest's 'view', and with the bundle in uihtml's Data
% (struct('bundle', <manifest, relative>, 'settings', settings)); it answers with the event 'ready' or 'error'
% (data.code 'no-webgpu' when it cannot draw).
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', 0.07*V, 'Faces', F);
%   I = exp(-sum((V - [0 0 1]).^2, 2) / 0.1) .* cos(linspace(0, 2*pi, 20));
%   out = rheome.view(S, I, Fs=10, CLim=[-1 1], View="dorsal", Viewer="show");
%   assert(strcmp(out.viewer, 'show'));  close(out.figure)
%
% See also: rheome.io.write.bundle, rheome.show.surface, rheome.show.vectors
%
% Author: Diellor Basha, 2026

    [surf, result, args] = splitInputs(varargin);
    p = inputParser;
    p.FunctionName = 'rheome.view';
    for k = settingNames(), p.addParameter(k{1}, []); end
    p.addParameter('Fs', 1);
    p.addParameter('Name', 'Rheome result');
    p.addParameter('Viewer', 'auto', @(x) any(strcmpi(x, {'auto', 'editor', 'app', 'show', 'none'})));
    p.addParameter('Timeout', 60);
    p.parse(args{:});
    o = p.Results;
    if ischar(surf) || isstring(surf), surf = rheome.io.read.surface(surf); end
    if isnumeric(result), result = struct('result', result); end
    T = size(result.(firstField(result)), 2);
    s = settings(o, surf, result, T);
    out = struct('viewer', lower(char(o.Viewer)), 'manifest', '', 'settings', s, 'figure', []);

    vdir = viewerFolder();
    page = fullfile(vdir, 'matlab.html');
    if strcmp(out.viewer, 'auto')
        if ~isMATLABReleaseOlderThan('R2026a') && isfile(page), out.viewer = 'editor'; else, out.viewer = 'app'; end
    end

    if ~strcmp(out.viewer, 'show')
        % one bundle at a time: the last one is open and loaded, older ones are dropped (new name: no stale cache)
        root = fullfile(vdir, 'bundles');
        if isfolder(root), rmdir(root, 's'); end
        id = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'));
        out.manifest = rheome.io.write.bundle(fullfile(root, id), surf, result, 'Fs', o.Fs, 'Name', o.Name, 'View', s);
    end

    switch out.viewer
        case 'editor'
            if ~isfile(page), error('rheome:view:noviewer', 'No viewer at %s (set RHEOME_VIEWER).', page); end
            [ok, out.figure] = viaEditor(page, ['bundles/' id '/manifest.json'], s, o.Timeout);
            if ~ok, out.viewer = 'show'; end
        case 'app'
            if ~viaApp(fileparts(out.manifest)), out.viewer = 'show'; end
    end
    if strcmp(out.viewer, 'show'), out.figure = viaShow(surf, result, s, o.Fs); end
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

function n = settingNames()
    n = {'Field', 'Colormap', 'CLim', 'Threshold', 'Surface', 'Hemisphere', 'View', 'Arrows', 'ArrowDensity', ...
        'ArrowScale', 'Overlays', 'Time', 'Play', 'Speed'};
end

function f = firstField(result)
    f = fieldnames(result);
    if isempty(f), error('rheome:view:input', 'The result has no field.'); end
    f = f{1};
end

function s = settings(o, surf, result, T)
% the options given, checked and named as the viewer's settings (rheome.bundle@1 'view')
    s = struct();
    fields = fieldnames(result);
    if ~isempty(o.Field), s.field = pick(o.Field, fields, 'Field'); end
    if ~isempty(o.Colormap)
        c = o.Colormap;
        if isnumeric(c)
            if size(c, 2) ~= 3 || any(c(:) < 0 | c(:) > 1), bad('Colormap', 'an [n x 3] matrix in [0, 1]'); end
            s.colormap = c;
        elseif ischar(c) || (isstring(c) && isscalar(c)), s.colormap = char(c);
        else, bad('Colormap', 'a name or an [n x 3] matrix');
        end
    end
    if ~isempty(o.CLim)
        c = double(o.CLim);
        if numel(c) ~= 2 || ~all(isfinite(c)) || c(1) >= c(2), bad('CLim', '[lo hi], lo < hi'); end
        s.clim = c(:)';
    end
    if ~isempty(o.Threshold), s.threshold = positive(o.Threshold, 'Threshold', true); end
    if ~isempty(o.Surface)
        if isfield(surf, 'Vertices'), names = {'default'}; else, names = fieldnames(surf); end
        s.surface = pick(o.Surface, names, 'Surface');
    end
    if ~isempty(o.Hemisphere)
        s.hemisphere = pick(o.Hemisphere, {'both', 'left', 'right'}, 'Hemisphere');
        S1 = surf;
        if ~isfield(S1, 'Vertices'), S1 = surf.(firstField(surf)); end
        if ~strcmp(s.hemisphere, 'both') && ~(isfield(S1, 'Hemi') && numel(S1.Hemi) == 2)
            error('rheome:view:option', 'Hemisphere "%s" needs the surface''s Structures split (.Hemi).', s.hemisphere);
        end
    end
    if ~isempty(o.View)
        v = o.View;
        if isnumeric(v)
            if numel(v) ~= 2 || ~all(isfinite(v)), bad('View', '[az el] degrees'); end
            s.camera = double(v(:)');
        else
            s.camera = pick(v, {'lateral', 'medial', 'dorsal', 'ventral', 'anterior', 'posterior'}, 'View');
        end
    end
    if ~isempty(o.Arrows) || ~isempty(o.ArrowDensity) || ~isempty(o.ArrowScale)
        s.arrows = struct();
        if ~isempty(o.Arrows), s.arrows.visible = flag(o.Arrows, 'Arrows'); end
        if ~isempty(o.ArrowDensity)
            d = positive(o.ArrowDensity, 'ArrowDensity', false);
            if d > 1, bad('ArrowDensity', 'in (0, 1]'); end
            s.arrows.density = d;
        end
        if ~isempty(o.ArrowScale), s.arrows.scale = positive(o.ArrowScale, 'ArrowScale', false); end
    end
    if ~isempty(o.Overlays)
        ov = cellstr(o.Overlays);
        for k = 1:numel(ov), ov{k} = pick(ov{k}, fields, 'Overlays'); end
        s.overlays = ov(:)';
    end
    if ~isempty(o.Time)
        t = o.Time;
        if ~isscalar(t) || ~isfinite(t) || t < 0 || t > (T - 1) / o.Fs, bad('Time', sprintf('in [0, %g] s', (T - 1) / o.Fs)); end
        s.time = double(t);
    end
    if ~isempty(o.Play), s.play = flag(o.Play, 'Play'); end
    if ~isempty(o.Speed), s.speed = positive(o.Speed, 'Speed', false); end
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

function vdir = viewerFolder()
% RHEOME_VIEWER, else prefdir/rheome-viewer, installed from <rheome>/viewer when that ships a newer copy
    vdir = getenv('RHEOME_VIEWER');
    if ~isempty(vdir), return; end
    vdir = fullfile(prefdir, 'rheome-viewer');
    shipped = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'viewer');
    if isfile(fullfile(shipped, 'matlab.html')) && ~strcmp(stamp(shipped), stamp(vdir))
        if isfolder(vdir), rmdir(vdir, 's'); end
        copyfile(shipped, vdir);
    end
end

function v = stamp(d)
    v = '';
    f = fullfile(d, 'VERSION');
    if isfile(f), v = strtrim(fileread(f)); end
end

function [ok, fig] = viaEditor(page, manifest, s, timeout)
% the viewer inside MATLAB; ok = false when it has no WebGPU
    persistent h
    if isempty(h) || ~isvalid(h) || ~strcmp(h.HTMLSource, page)
        fig = uifigure('Name', 'Rheome - viewer', 'Position', [100 100 1400 900]);
        g = uigridlayout(fig, [1 1], 'Padding', 0);
        h = uihtml(g, 'HTMLSource', page, 'HTMLEventReceivedFcn', @onEvent);
    end
    fig = ancestor(h, 'figure');
    figure(fig);
    h.UserData = [];
    h.Data = struct('bundle', manifest, 'settings', s);   % a reused figure: DataChanged opens it in place of the last
    t = tic;
    while isempty(h.UserData) && toc(t) < timeout && isvalid(h), pause(0.1); end
    ok = true;
    if ~isvalid(h), return; end
    ev = h.UserData;
    if isempty(ev)
        warning('rheome:view:timeout', 'The viewer did not say ready within %g s; the figure is left open.', timeout);
    elseif strcmp(ev.name, 'error')
        d = ev.data;
        if isstruct(d) && isfield(d, 'code') && strcmp(d.code, 'no-webgpu')
            delete(fig); fig = []; ok = false;
            return
        end
        msg = '';
        if isstruct(d) && isfield(d, 'message'), msg = d.message; end
        error('rheome:view:editor', 'The viewer could not open the bundle: %s', msg);
    end
end

function onEvent(src, ev)
    if any(strcmp(ev.HTMLEventName, {'ready', 'error'}))
        src.UserData = struct('name', ev.HTMLEventName, 'data', ev.HTMLEventData);
    end
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

function fig = viaShow(surf, result, s, fs)
% no viewer: the chosen field at the chosen time in a MATLAB figure, with the settings rheome.show has
    if ~isfield(surf, 'Vertices')
        names = fieldnames(surf); k = 1;
        if isfield(s, 'surface'), k = find(strcmp(names, s.surface)); end
        S1 = surf.(names{1}); surf = surf.(names{k}); surf.Faces = S1.Faces;
    end
    f = firstField(result);
    if isfield(s, 'field'), f = s.field; end
    x = result.(f);
    t = 1;
    if isfield(s, 'time'), t = round(s.time * fs) + 1; end
    opt = {'Title', sprintf('%s, frame %d of %d (rheome.show)', f, t, size(x, 2))};
    if isfield(s, 'clim'), opt = [opt {'CLim', s.clim}]; end
    if isfield(s, 'camera')
        c = s.camera;
        if ischar(c)
            c = struct('lateral', [-90 0], 'medial', [90 0], 'dorsal', [0 90], 'ventral', [0 -90], ...
                'anterior', [180 0], 'posterior', [0 0]).(c);
        end
        opt = [opt {'View', c}];
    end
    if isfield(s, 'colormap')
        c = s.colormap;
        if ischar(c) && exist(c, 'file') == 2, c = feval(c, 256); end   % a MATLAB colormap name; else the default
        if isnumeric(c), opt = [opt {'Colormap', c}]; end
    end
    if size(x, 1) == size(surf.Vertices, 1)
        ax = rheome.show.surface(surf, x(:, t), opt{:});
    else
        ax = rheome.show.vectors(surf, x(:, t), opt{:});
    end
    fig = ancestor(ax, 'figure');
end
