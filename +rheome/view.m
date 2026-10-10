function out = view(surf, result, varargin)
% VIEW  Show a result on a surface in the Cortical Flow Editor, from one call; rheome.show when there is no Editor.
%
%   rheome.view(surf, result)
%   rheome.view(surf, struct('activity', I, 'velocity', v), 'Fs', 100)
%   out = rheome.view(..., name, value)
%
% The result is written as a rheome.bundle@1 folder (rheome.io.write.bundle) and opened in the Editor:
%   R2026a+       inside MATLAB: a uifigure + uihtml on the viewer's matlab.html (the Editor's MATLAB build). The
%                 viewer is a folder you install (RHEOME_VIEWER, default prefdir/rheome-viewer) and the bundle is
%                 written beneath it, because uihtml serves only that folder's tree. The figure is reused.
%   R2023b-R2025b the Cortical Flow desktop app, by its URL (corticalflow://view?bundle=<folder>); one window, reused.
% With no viewer, no WebGPU (the Editor's 'no-webgpu' error) or no app, the first field's first frame is drawn with
% rheome.show.surface (real) or rheome.show.vectors (vector). Nothing else in Rheome needs the viewer.
%
% INPUTS:
%   surf   : surface struct (.Vertices meters, .Faces 1-based) or a Brainstorm surface .mat path
%   result : [nV x T] or [3nV x T] array, or a struct of such arrays (see rheome.io.write.bundle)
%
% OPTIONS (name/value):
%   'Fs'      sampling rate of the T columns, Hz (default 1)
%   'Name'    the bundle's name (default 'Rheome result')
%   'Viewer'  'auto' (default) | 'editor' (uihtml) | 'app' | 'show' | 'none' (write the bundle only)
%   'Timeout' seconds to wait for the Editor's 'ready' (default 60)
%
% OUTPUT:
%   out : struct .viewer (the path taken: 'editor' | 'app' | 'show' | 'none'), .manifest (the bundle's
%         manifest.json, '' for 'show'), .figure (the uifigure or figure, [] for 'app' and 'none').
%
% See also: rheome.io.write.bundle, rheome.show.surface, rheome.show.vectors
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Fs', 1);
    p.addParameter('Name', 'Rheome result');
    p.addParameter('Viewer', 'auto', @(x) any(strcmpi(x, {'auto', 'editor', 'app', 'show', 'none'})));
    p.addParameter('Timeout', 60);
    p.parse(varargin{:});
    o = p.Results;
    if ischar(surf) || isstring(surf), surf = rheome.io.read.surface(surf); end
    if isnumeric(result), result = struct('result', result); end
    out = struct('viewer', lower(char(o.Viewer)), 'manifest', '', 'figure', []);

    vdir = getenv('RHEOME_VIEWER');
    if isempty(vdir), vdir = fullfile(prefdir, 'rheome-viewer'); end
    page = fullfile(vdir, 'matlab.html');
    if strcmp(out.viewer, 'auto')
        if ~isMATLABReleaseOlderThan('R2026a') && isfile(page), out.viewer = 'editor'; else, out.viewer = 'app'; end
    end

    if ~strcmp(out.viewer, 'show')
        % one bundle at a time: the last one is open and loaded, older ones are dropped (new name: no stale cache)
        root = fullfile(vdir, 'bundles');
        if isfolder(root), rmdir(root, 's'); end
        id = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'));
        out.manifest = rheome.io.write.bundle(fullfile(root, id), surf, result, 'Fs', o.Fs, 'Name', o.Name);
    end

    switch out.viewer
        case 'editor'
            if ~isfile(page), error('rheome:view:noviewer', 'No viewer at %s (set RHEOME_VIEWER).', page); end
            [ok, out.figure] = viaEditor(page, ['bundles/' id '/manifest.json'], o.Timeout);
            if ~ok, out.viewer = 'show'; end
        case 'app'
            if ~viaApp(fileparts(out.manifest)), out.viewer = 'show'; end
    end
    if strcmp(out.viewer, 'show'), out.figure = viaShow(surf, result); end
end

function [ok, fig] = viaEditor(page, manifest, timeout)
% the Editor inside MATLAB; ok = false when it has no WebGPU
    persistent h
    if isempty(h) || ~isvalid(h)
        fig = uifigure('Name', 'Rheome - Cortical Flow Editor', 'Position', [100 100 1400 900]);
        g = uigridlayout(fig, [1 1], 'Padding', 0);
        h = uihtml(g, 'HTMLSource', page, 'HTMLEventReceivedFcn', @onEvent);
    end
    fig = ancestor(h, 'figure');
    figure(fig);
    h.UserData = [];
    h.Data = struct('bundle', manifest);   % a reused figure: DataChanged opens it in place of the last
    t = tic;
    while isempty(h.UserData) && toc(t) < timeout && isvalid(h), pause(0.1); end
    ok = true;
    if ~isvalid(h), return; end
    ev = h.UserData;
    if isempty(ev)
        warning('rheome:view:timeout', 'The Editor did not say ready within %g s; the figure is left open.', timeout);
    elseif strcmp(ev.name, 'error')
        d = ev.data;
        if isstruct(d) && isfield(d, 'code') && strcmp(d.code, 'no-webgpu')
            delete(fig); fig = []; ok = false;
            return
        end
        msg = '';
        if isstruct(d) && isfield(d, 'message'), msg = d.message; end
        error('rheome:view:editor', 'The Editor could not open the bundle: %s', msg);
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

function fig = viaShow(surf, result)
% no Editor: the first field's first frame in a MATLAB figure
    names = fieldnames(result);
    x = result.(names{1});
    ttl = sprintf('%s, frame 1 of %d (rheome.show)', names{1}, size(x, 2));
    if size(x, 1) == size(surf.Vertices, 1)
        ax = rheome.show.surface(surf, x(:, 1), 'Title', ttl);
    else
        ax = rheome.show.vectors(surf, x(:, 1), 'Title', ttl);
    end
    fig = ancestor(ax, 'figure');
end
