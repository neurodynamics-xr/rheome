function out = snapshot(surf, result, varargin)
% SNAPSHOT  Print-ready cortical panels from the Cortical Flow viewer, headless, by the lab's figure tokens.
%
%   out = rheome.snapshot(surf, result, Shots=shots, Out=folder)
%   out = rheome.snapshot(surf, result, Shots=struct('Field', "div", 'Unit', "1/s", 'Label', "divergence"), Out="figs")
%
% Every cortical panel of a figure the same way: the viewer (nxr-cortical-flow editor/dist-viewer, the renderer of
% rheome.view) draws each view in headless Chrome at print resolution, and MATLAB lays the views out with the
% figure tokens (library/styles/tokens/figure-tokens.json): the colormap and limits by quantity, the colourbar with
% quantity, unit and reference, view labels, a geodesic scale bar on the folded surface, and for a time series one
% strip with the interval stated once. Each shot is written as <Out>/<Name>.svg (text, bars and colourbar editable;
% the views embedded at their rendered resolution), printed by Chrome to <Out>/<Name>.pdf (fonts embedded, the views
% at full resolution: MATLAB's own PDF export resamples images to screen pixels) and .png (at the journal's
% combination dpi), and the viewer's raw renders are kept as <Out>/<Name>/*.png (for boards composed elsewhere), each with
% its marks sidecar <render>.marks.json (every overlay drawn and its key: figure_check.py E1, rheome.viz.marks).
% THE PRIMARY API for printed cortex, with rheome.view (interactive); see rheome.viz.
%
% THE SHOT LIST -- Shots is a struct array or a cell array of structs, one panel each, with the fields:
%   Field      the real field painted (a field of result)                                  (default: the first real)
%   Unit       its unit, as printed on the colourbar ("a.u." only when it has none)        REQUIRED
%   Label      the quantity's name on the colourbar                                         (default: Field)
%   Reference  what a level is relative to ("4-min median"): "quantity (unit re reference)" (default: none)
%   Quantity   "signed" | "magnitude" | "power" | "phase" | "coupling" | "code": the token's colormap and limits
%              (default: "signed" when the field has both signs, else "magnitude")
%   CLim       [lo hi] in the field's units (default by quantity: signed/coupling +-p95 of |x|, magnitude [0 p99],
%              power [p1 p99], phase [-pi pi], code [min max], over the frames and hemispheres shown)
%   Colormap   a viewer colormap name, in place of the quantity's token
%   Time       s from the first sample; a vector of equally spaced times is a strip (default 0)
%   Frames     1-based frames, in place of Time
%   View       "standard" (L lateral, L medial, R medial, R lateral: the figure system's four), one of those four,
%              or a viewer camera "left" | "right" | "top" | "bottom" | "front" | "back" with Hemisphere
%              (default "standard"; a strip takes one view)
%   Hemisphere "both" | "left" | "right", with a viewer camera (default "both")
%   Surface    the surface shown, a name of surf's struct (default: the first)
%   Vector     a vector field drawn as arrows (default none); arrows show the same vector the colour shows
%   Settings   struct of further viewer settings, applied last (README of the viewer, 'The settings')
%   Name       the files' name (default "shot<k>")
%
% OPTIONS (name/value):
%   'Out'       folder written (created)                                                     REQUIRED
%   'Fs'        sampling rate of result's columns, Hz (default 1)
%   'Width'     the panel's width: "single" | "onehalf" | "double" (the token journal's widths) or mm
%               (default: "double" for four views or a strip, "single" for one view)
%   'ScaleBar'  the scale bar's length, mm along the folded surface; 0: none (default 20)
%   'Folded'    the name of the folded surface in surf, on which the bar is measured (default "pial" when surf has
%               it, else the first); on another surface shown, the legend says the distance is on the folded one
%   'Format'    printed from the SVG (always written): "pdf", "png", both (default) or [] (the SVG only)
%   'Tokens'    the figure tokens' JSON (default $NXR_FIGURE_TOKENS, else ~/workspace/library/styles/tokens/
%               figure-tokens.json)
%   'Viewer'    the viewer folder (default as rheome.view: $RHEOME_VIEWER, else prefdir/rheome-viewer)
%   'Chrome'    Chrome or Chromium (default $RHEOME_CHROME, else the usual install)
%   'Node'      Node.js 22+ (default $RHEOME_NODE, else the usual install)
%   'Labellings' partitions for the viewer (rheome.viz.tiles), shown by a shot's Settings.labelling
%   'Points'    point sets for the viewer (rheome.viz.points: critical points), shown by Settings.points (default all)
%   'Scale'     the viewer's least pixel ratio (default: the tokens' cortex.viewer.pixelRatio); raised per shot so the
%               views print at the journal's combination dpi
%
% OUTPUT:
%   out : struct array, one per shot: .name .svg .pdf .png .renders (the viewer's PNGs) .clim .colormap .quantity
%         .colorbar (its label) .views .times (s) .dpi (resolution of the views as printed) .scalebar (.mm, .view,
%         .length_mm measured along the drawn curve, .vertices) .legend (the sentence the figure legend needs:
%         views, surface, the bar's distance, the interval) .settings (each render's viewer settings)
%
% THE SCALE BAR is a curve on the surface, not a line on the image: the shortest path along the mesh edges of the
% folded surface (rheome.geom.edgegraph, no diagonals, so the curve lies on the mesh), cut at exactly 'ScaleBar' mm
% of that path's length, and drawn where it lies on the view, a visible stretch near the view's lower edge. Its
% length is the drawn curve's length on the folded surface; a mesh-edge path is 3-13 % longer than the geodesic
% between its ends (rheome.geom.edgegraph), so its ends are a little closer than the label says, never farther.
% The camera is read from each render (the viewer's stage is orthographic): the silhouette of the shown vertices is
% fitted to the PNG's alpha and checked (rheome:snapshot:camera when it does not fit).
%
% Needs the viewer (WebGPU), Chrome and Node; nothing else in Rheome does.
%
% Example:
%   S = rheome.load.surface('omega_sub0002');  x = S.Vertices(:, 1) * 1e3;
%   out = rheome.snapshot(S, struct('x', x), Out=tempname, ...
%       Shots=struct('Field', "x", 'Label', "anterior position", 'Unit', "mm"));
%   open(out.pdf)
%
% See also: rheome.view, rheome.io.write.bundle, rheome.geom.edgegraph
%
% Author: Diellor Basha, 2026

    p = inputParser;  p.FunctionName = 'rheome.snapshot';
    p.addParameter('Shots', {});
    p.addParameter('Out', '');
    p.addParameter('Fs', 1, @(x) isscalar(x) && x > 0);
    p.addParameter('Width', []);
    p.addParameter('ScaleBar', 20, @(x) isscalar(x) && x >= 0);
    p.addParameter('Folded', '');
    p.addParameter('Format', ["pdf" "png"]);
    p.addParameter('Tokens', '');
    p.addParameter('Viewer', '');
    p.addParameter('Chrome', '');
    p.addParameter('Node', '');
    p.addParameter('Scale', []);
    p.addParameter('Wait', 600);
    p.addParameter('Labellings', []);
    p.addParameter('Points', []);
    p.addParameter('Timeout', 60);
    p.parse(varargin{:});
    o = p.Results;
    if isempty(o.Out), error('rheome:snapshot:input', 'Out: the folder written is expected.'); end
    o.Out = char(o.Out);  o.Format = lower(string(o.Format));

    % the inputs, as rheome.view takes them
    if ischar(surf) || isstring(surf), surf = rheome.io.read.surface(surf); end
    if isnumeric(result), result = struct('result', result); end
    if isfield(surf, 'Vertices'), surfs = struct('default', surf); else, surfs = surf; end
    snames = fieldnames(surfs);
    F = surfs.(snames{1}).Faces;
    T = size(result.(firstOf(fieldnames(result))), 2);

    tok = tokens(o.Tokens);
    if isempty(o.Scale), o.Scale = tok.cortex.viewer.pixelRatio; end
    vdir = char(o.Viewer);
    if isempty(vdir), vdir = viewerFolder(); end
    if ~isfile(fullfile(vdir, 'viewer.html')), error('rheome:snapshot:noviewer', 'No viewer at %s (Viewer= or RHEOME_VIEWER).', vdir); end
    cmaps = viewerColormaps(vdir);

    shots = normalise(o.Shots, result, surfs, T, o.Fs, tok, cmaps);
    hemi = hemispheres(surfs.(snames{1}));
    folded = o.Folded;
    if isempty(folded), folded = snames{1}; if any(strcmp(snames, 'pial')), folded = 'pial'; end, end
    folded = pickName(folded, snames, 'Folded');

    % one bundle, one viewer run for every render
    work = tempname;  mkdir(work);  cleanWork = onCleanup(@() rmdir(work, 's'));
    if isfield(surf, 'Vertices'), bs = surf; else, bs = surfs; end
    rheome.io.write.bundle(fullfile(work, 'bundle'), bs, result, 'Fs', o.Fs, 'Name', 'rheome.snapshot', ...
        'Labellings', o.Labellings, 'Points', o.Points);
    if ~isfolder(o.Out), mkdir(o.Out); end
    renders = {};
    J = tok.journal.(tok.journal.default);
    for k = 1:numel(shots)
        shots(k).dir = fullfile(o.Out, shots(k).Name);
        if ~isfolder(shots(k).dir), mkdir(shots(k).dir); end
        % the pixel ratio that prints the views at the journal's dpi: a cortex fills about half the render's width
        n = numel(shots(k).renders);  W = panelWidth(o.Width, n, J);
        need = ((W - (n - 1) * tok.spacing.viewGap_mm) / n) / 25.4 * J.dpi.combination / (0.5 * tok.cortex.viewer.size_px(1));
        scale = max(o.Scale, ceil(10 * need) / 10);
        for r = 1:numel(shots(k).renders)
            renders{end+1} = struct('name', shots(k).renders(r).file, 'settings', shots(k).renders(r).settings, 'scale', scale); %#ok<AGROW>
        end
    end
    job = struct('viewer', vdir, 'bundle', fullfile(work, 'bundle'), 'out', work, 'chrome', chromePath(o.Chrome), ...
        'size', tok.cortex.viewer.size_px(:)', 'scale', o.Scale, 'wait', o.Wait, 'timeout', o.Timeout, 'shots', {renders});
    runDriver(job, nodePath(o.Node), work);

    graphCache = [];
    out = struct([]);  prints = {};
    for k = 1:numel(shots)
        s = shots(k);
        for r = 1:numel(s.renders)
            movefile(fullfile(work, [s.renders(r).file '.png']), fullfile(s.dir, [s.renders(r).label '.png']), 'f');
            mk = fullfile(work, [s.renders(r).file '.marks.json']);   % the render's overlays and keys (rheome.viz.marks)
            if isfile(mk), movefile(mk, fullfile(s.dir, [s.renders(r).label '.marks.json']), 'f'); end
        end
        [o1, graphCache, pj] = compose(s, surfs, F, hemi, folded, graphCache, tok, cmaps, o);
        out = [out o1]; %#ok<AGROW>
        if ~isempty(pj.pdf) || ~isempty(pj.png), prints{end+1} = pj; end %#ok<AGROW>
    end
    if ~isempty(prints)   % the SVGs printed by Chrome: MATLAB's PDF export resamples images to screen pixels
        runDriver(struct('chrome', job.chrome, 'print', {prints}), nodePath(o.Node), work);
    end
end

%% ---------------------------------------------------------------- the shots

function shots = normalise(in, result, surfs, T, fs, tok, cmaps)
% each shot checked, its defaults from the tokens, and its renders (one per view or frame) with full viewer settings
    if isempty(in), in = {struct()}; end
    if isstruct(in), in = num2cell(in); end
    if ~iscell(in), error('rheome:snapshot:input', 'Shots: a struct array or a cell array of structs is expected.'); end
    known = {'Field', 'Unit', 'Label', 'Reference', 'Quantity', 'CLim', 'Colormap', 'Time', 'Frames', 'View', ...
        'Hemisphere', 'Surface', 'Vector', 'Settings', 'Name'};
    fields = fieldnames(result);
    nV = size(surfs.(firstOf(fieldnames(surfs))).Vertices, 1);
    real = fields(cellfun(@(f) size(result.(f), 1) == nV, fields));
    vec = fields(cellfun(@(f) size(result.(f), 1) == 3*nV, fields));
    hemi = [];
    for k = numel(in):-1:1
        a = in{k};
        if ~isstruct(a) || ~isscalar(a), error('rheome:snapshot:input', 'Shots{%d}: a struct is expected.', k); end
        bad = setdiff(fieldnames(a), known);
        if ~isempty(bad), error('rheome:snapshot:shot', 'Shots{%d}: %s is not a shot field (%s).', k, strjoin(bad, ', '), strjoin(known, ', ')); end
        g = @(n, d) get(a, n, d);
        s.Name = char(g('Name', sprintf('shot%d', k)));
        if isempty(regexp(s.Name, '^[\w\-\.]+$', 'once')), error('rheome:snapshot:shot', 'Name: letters, digits, - _ . only (%s).', s.Name); end
        if isempty(real), error('rheome:snapshot:shot', 'The result has no real ([nV x T]) field to paint.'); end
        s.Field = pickName(g('Field', real{1}), real, 'Field');
        s.Unit = char(g('Unit', ''));
        if isempty(s.Unit), error('rheome:snapshot:unit', '%s: a colourbar needs a unit (Unit="a.u." only when there is none; say why in the legend).', s.Name); end
        s.Label = char(g('Label', s.Field));
        s.Reference = char(g('Reference', ''));
        x = result.(s.Field);

        % frames: 1-based columns; a strip is equally spaced
        if isfield(a, 'Frames') && isfield(a, 'Time'), error('rheome:snapshot:shot', '%s: Time or Frames, not both.', s.Name); end
        if isfield(a, 'Frames'), fr = double(a.Frames(:)');
        else, fr = round(double(g('Time', 0)) * fs + 1); fr = fr(:)';
        end
        if isempty(fr) || any(fr < 1 | fr > T | fr ~= round(fr)), error('rheome:snapshot:shot', '%s: frames in 1..%d (times in [0, %g] s).', s.Name, T, (T-1)/fs); end
        if numel(fr) > 2 && any(diff(fr, 2)), error('rheome:snapshot:shot', '%s: a strip needs equally spaced times (its interval is stated once).', s.Name); end
        if numel(fr) > 1 && any(diff(fr) <= 0), error('rheome:snapshot:shot', '%s: times must increase.', s.Name); end
        if numel(fr) > tok.dynamics.strip.frames_range(2)
            warning('rheome:snapshot:frames', '%s: %d frames; the tokens ask for %d-%d in a strip, and a video beyond.', s.Name, numel(fr), tok.dynamics.strip.frames_range);
        end
        s.frames = fr;  s.times = (fr - 1) / fs;

        % views: the figure system's four, one of them, or a viewer camera with a hemisphere
        v = char(g('View', 'standard'));
        views4 = {'L lateral', 'left', 'left'; 'L medial', 'left', 'right'; 'R medial', 'right', 'left'; 'R lateral', 'right', 'right'};
        if strcmpi(v, 'standard'), s.views = views4;
        elseif any(strcmpi(v, views4(:, 1))), s.views = views4(strcmpi(v, views4(:, 1)), :);
        else
            cam = pickName(v, {'left', 'right', 'top', 'bottom', 'front', 'back'}, 'View');
            h = pickName(g('Hemisphere', 'both'), {'both', 'left', 'right'}, 'Hemisphere');
            s.views = {sprintf('%s, %s', h, cam), h, cam};
        end
        if numel(fr) > 1 && size(s.views, 1) > 1, error('rheome:snapshot:shot', '%s: a strip takes one view (one view and one scale for the whole strip).', s.Name); end
        s.Surface = pickName(g('Surface', firstOf(fieldnames(surfs))), fieldnames(surfs), 'Surface');

        % quantity -> colormap and limits, by the tokens
        if isempty(hemi), hemi = hemispheres(surfs.(firstOf(fieldnames(surfs)))); end
        shown = ismember(hemi, find(ismember({'left', 'right'}, s.views(:, 2)))) | any(strcmp(s.views(:, 2), 'both'));
        xs = x(shown, fr);  xs = xs(isfinite(xs));
        q = char(g('Quantity', ''));
        if isempty(q), if any(xs < 0) && any(xs > 0), q = 'signed'; else, q = 'magnitude'; end, end
        s.Quantity = pickName(q, {'signed', 'magnitude', 'power', 'phase', 'coupling', 'code'}, 'Quantity');
        bq = tok.colour.byQuantity.(s.Quantity);
        tokenName = regexp(bq.colormap, '^[\w\-]+', 'match', 'once');
        s.Colormap = char(g('Colormap', tok.colour.colormap.(strrep(tokenName, '-', '_')).x_value));
        if ~isKey(cmaps, s.Colormap)
            error('rheome:snapshot:colormap', '%s: the viewer has no colormap ''%s'' (%s quantity); give Colormap= or add it to the viewer.', s.Name, s.Colormap, s.Quantity);
        end
        diverging = any(strcmp(s.Quantity, {'signed', 'coupling'}));
        c = g('CLim', []);
        if isempty(c)
            switch s.Quantity
                case {'signed', 'coupling'}, L = prctile(abs(xs), 95); c = [-L L];
                case 'magnitude', c = [0 prctile(xs, 99)];
                case 'power', c = prctile(xs, [1 99]);
                case 'phase', c = [-pi pi];
                case 'code', c = [min(xs) max(xs)];
            end
            if ~(c(2) > c(1)), c = c(1) + [-1 1] * max(1, abs(c(1))) * eps(1e6); end   % a constant field: any non-empty range
        end
        c = double(c(:)');
        if numel(c) ~= 2 || ~all(isfinite(c)) || c(1) >= c(2), error('rheome:snapshot:shot', '%s: CLim [lo hi], lo < hi.', s.Name); end
        if diverging && abs(c(1) + c(2)) > 1e-9 * max(abs(c)), error('rheome:snapshot:clim', '%s: a diverging scale is symmetric around 0 (got [%g %g]).', s.Name, c(1), c(2)); end
        s.CLim = c;  s.diverging = diverging;
        s.Vector = char(g('Vector', ''));
        if ~isempty(s.Vector), s.Vector = pickName(s.Vector, vec, 'Vector'); end
        extra = g('Settings', struct());
        if ~isstruct(extra), error('rheome:snapshot:shot', '%s: Settings is a struct of viewer settings.', s.Name); end

        % the viewer's settings, every one that a render depends on (each render is set whole)
        base = struct('controls', false, 'background', 'transparent', 'scaleBar', false, 'opacity', 1, 'threshold', NaN, ...   % NaN: JSON null, everywhere
            'labelling', '', 'points', {{}}, 'scalar', s.Field, 'vector', s.Vector, 'colormap', s.Colormap, 'clim', c, ...
            'surface', s.Surface, 'isolines', 0, 'halfMax', false, 'playing', false);
        if strcmp(s.Quantity, 'signed'), base.isolines = tok.cortex.viewer.isolines; end
        if ~isempty(s.Vector)
            base.arrowDensity = tok.cortex.viewer.arrowDensity;
            if diverging, base.arrowDensity = tok.cortex.viewer.arrowDensitySigned; end
            base.arrowColor = 'ink';
            if any(strcmp(s.Quantity, {'magnitude', 'power'})), base.arrowColor = 'paper'; end   % paper on a dark map
        end
        for f = fieldnames(extra)', base.(f{1}) = extra.(f{1}); end
        rs = struct('file', {}, 'label', {}, 'view', {}, 'hemisphere', {}, 'camera', {}, 'frame', {}, 'settings', {});
        for vi = 1:size(s.views, 1)
            for fi = 1:numel(fr)
                st = base;  st.hemisphere = s.views{vi, 2};  st.view = s.views{vi, 3};  st.frame = fr(fi) - 1;
                lab = regexprep(lower(s.views{vi, 1}), '[^a-z0-9]+', '-');
                if numel(fr) > 1, lab = sprintf('%s_f%d', lab, fr(fi)); end
                rs(end+1) = struct('file', sprintf('s%d_%d_%d', k, vi, fi), 'label', lab, 'view', s.views{vi, 1}, ...
                    'hemisphere', s.views{vi, 2}, 'camera', s.views{vi, 3}, 'frame', fr(fi), 'settings', st); %#ok<AGROW>
            end
        end
        s.renders = rs;  s.dir = '';
        shots(k) = s;
    end
end

function v = get(a, n, d)
    if isfield(a, n) && ~isempty(a.(n)), v = a.(n); else, v = d; end
end

function x = pickName(x, allowed, name)
    if ~(ischar(x) || (isstring(x) && isscalar(x))), error('rheome:snapshot:option', '%s: one of %s is expected.', name, strjoin(allowed, ', ')); end
    k = find(strcmpi(char(x), allowed), 1);
    if isempty(k), error('rheome:snapshot:option', '%s: ''%s'' is not one of %s.', name, char(x), strjoin(allowed, ', ')); end
    x = allowed{k};
end

function f = firstOf(c)
    f = c{1};
end

%% ---------------------------------------------------------------- tokens, colormaps, tools

function T = tokens(file)
% the lab's figure tokens (Graphics; library/styles/tokens/figure-tokens.json)
    if isstruct(file), T = file; return; end
    file = char(file);
    if isempty(file), file = getenv('NXR_FIGURE_TOKENS'); end
    if isempty(file), file = fullfile(char(java.lang.System.getProperty('user.home')), 'workspace', 'library', 'styles', 'tokens', 'figure-tokens.json'); end
    if ~isfile(file), error('rheome:snapshot:tokens', 'No figure tokens at %s (Tokens= or NXR_FIGURE_TOKENS).', file); end
    T = jsondecode(fileread(file));
end

function m = viewerColormaps(vdir)
% the viewer's own colour tables (its default library), so the colourbar is the render's colours
    d = fullfile(vdir, 'colormaps', 'default');
    j = jsondecode(fileread(fullfile(d, 'colormaps.json')));
    [img, ~, alpha] = imread(fullfile(d, 'colormaps.png')); %#ok<ASGLU>
    m = containers.Map();
    for k = 1:numel(j.colormaps)
        c = j.colormaps(k);  if iscell(c), c = c{1}; end
        m(c.name) = double(squeeze(img(c.index + 1, :, 1:3))) / 255;
    end
end

function c = chromePath(c)
    c = char(c);
    if isempty(c), c = getenv('RHEOME_CHROME'); end
    if isempty(c)
        for p = {'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '/Applications/Chromium.app/Contents/MacOS/Chromium', ...
                '/usr/bin/google-chrome', '/usr/bin/chromium', '/usr/bin/chromium-browser'}
            if isfile(p{1}), c = p{1}; break; end
        end
    end
    if isempty(c) || ~isfile(c), error('rheome:snapshot:chrome', 'No Chrome found (Chrome= or RHEOME_CHROME).'); end
end

function n = nodePath(n)
    n = char(n);
    if isempty(n), n = getenv('RHEOME_NODE'); end
    if isempty(n)
        for p = {'/opt/homebrew/bin/node', '/usr/local/bin/node', '/usr/bin/node'}
            if isfile(p{1}), n = p{1}; break; end
        end
    end
    if isempty(n)
        [st, w] = system('command -v node');
        if st == 0, n = strtrim(w); end
    end
    if isempty(n) || ~isfile(n), error('rheome:snapshot:node', 'No Node.js found (Node= or RHEOME_NODE).'); end
end

function runDriver(job, node, work)
    jf = fullfile(work, 'job.json');
    fid = fopen(jf, 'w');  fwrite(fid, jsonencode(job));  fclose(fid);
    drv = fullfile(fileparts(mfilename('fullpath')), 'private', 'snapshot_driver.mjs');
    [st, msg] = system(sprintf('"%s" "%s" "%s" 2>&1', node, drv, jf));
    if st ~= 0, error('rheome:snapshot:viewer', 'The viewer run failed: %s', strtrim(msg)); end
end

%% ---------------------------------------------------------------- the geometry: hemispheres, camera, visibility, bar

function h = hemispheres(S)
% 1 left, 2 right, 0 neither: the mesh's two largest connected components, the left the one with larger mean y (as
% the viewer finds them)
    F = double(S.Faces);  nV = size(S.Vertices, 1);
    g = graph([F(:,1); F(:,2); F(:,3)], [F(:,2); F(:,3); F(:,1)], [], nV);
    c = conncomp(g)';
    n = accumarray(c, 1);
    [~, o] = sort(n, 'descend');
    h = zeros(nV, 1);
    if numel(o) < 2, h(:) = 1; return; end
    a = c == o(1);  b = c == o(2);
    if mean(S.Vertices(a, 2)) >= mean(S.Vertices(b, 2)), h(a) = 1; h(b) = 2; else, h(a) = 2; h(b) = 1; end
end

function cam = fitCamera(alpha, X, camera)
% the orthographic camera of a render, read from it: the shown vertices' projection fitted to the PNG's alpha.
% Candidates: the four proper rotations about the view direction (right x up = toward the camera: no mirror).
    d = struct('left', [0 1 0], 'right', [0 -1 0], 'top', [0 0 1], 'bottom', [0 0 -1], 'front', [1 0 0], 'back', [-1 0 0]).(camera);
    ax = eye(3);  ax = ax(abs(ax * d') < 0.5, :);
    m = alpha > 0;
    [r, c] = find(m);
    if isempty(r), error('rheome:snapshot:camera', 'The render is empty.'); end
    box = [min(c) max(c) min(r) max(r)];
    best = [];  scores = [];
    for right = {ax(1,:), -ax(1,:), ax(2,:), -ax(2,:)}
        R = [right{1}; -cross(d, right{1})];          % rows: screen x (right), screen y (down) = -up, up = d x right
        uv = X * R';
        ext = max(uv) - min(uv);
        sx = (box(2) - box(1) + 1) / ext(1);  sy = (box(4) - box(3) + 1) / ext(2);
        s = (sx + sy) / 2;
        off = [box(1) - 0.5 box(3) - 0.5] - s * min(uv);
        px = s * uv + off;
        in = m(sub2ind(size(m), min(max(round(px(:, 2)), 1), size(m, 1)), min(max(round(px(:, 1)), 1), size(m, 2))));
        score = mean(in) - abs(sx / sy - 1);
        scores(end+1) = score; %#ok<AGROW>
        if isempty(best) || score > best.score, best = struct('R', R, 's', s, 'off', off, 'score', score, 'inside', mean(in), 'aspect', sx / sy, 'd', d); end
    end
    srt = sort(scores, 'descend');
    if best.inside < 0.97 || abs(best.aspect - 1) > 0.03 || srt(1) - srt(2) < 0.01
        error('rheome:snapshot:camera', 'The render''s camera could not be read (%.1f %% of vertices inside, aspect %.3f, margin %.3f).', ...
            100 * best.inside, best.aspect, srt(1) - srt(2));
    end
    cam = best;
end

function px = project(cam, X)
    px = cam.s * (X * cam.R') + cam.off;
end

function vis = visible(Q, d, V, F)
% exact: a point is visible when the ray from it toward the camera (+d) hits no triangle (Moller-Trumbore)
    A = V(F(:,1), :);  e1 = V(F(:,2), :) - A;  e2 = V(F(:,3), :) - A;
    pv = cross(repmat(d, size(F, 1), 1), e2, 2);  det = sum(e1 .* pv, 2);
    vis = true(size(Q, 1), 1);
    for k = 1:size(Q, 1)
        t0 = Q(k, :) - A;
        u = sum(t0 .* pv, 2) ./ det;
        qv = cross(t0, e1, 2);
        v = (qv * d') ./ det;
        t = sum(e2 .* qv, 2) ./ det;
        hit = abs(det) > 1e-18 & u >= 0 & v >= 0 & u + v <= 1 & t > 1e-4;   % beyond 0.1 mm: not its own triangles
        vis(k) = ~any(hit);
    end
end

function [bar, G] = scaleBar(Lmm, Sf, Vd, F, shown, cam, alpha, G)
% a curve on the surface: the mesh-edge shortest path on the folded surface, cut at Lmm of its length, at a visible
% stretch near the lower edge of the view; drawn at the shown surface's positions (Vd)
    bar = [];
    if isempty(G), G = rheome.geom.edgegraph(Sf, Diagonals=false); end
    L = Lmm / 1e3;
    px = project(cam, Vd);
    idx = find(shown);
    % rough front-most test for the candidates (exact test on the chosen curve below): the nearest to the camera in its cell
    depth = Vd * cam.d';
    elen = median(G.Edges.Weight) * cam.s;
    csz = max(4, 2 * elen);
    key = floor(px(idx, :) / csz);  [~, ~, ci] = unique(key, 'rows');
    front = accumarray(ci, depth(idx), [], @max);
    cand = idx(depth(idx) >= front(ci) - 1.5e-3);
    [r, c] = find(alpha > 0);
    target = [min(c) + 0.7 * (max(c) - min(c)), min(r) + 0.88 * (max(r) - min(r))];
    [~, o] = sort(vecnorm(px(cand, :) - target, 2, 2));
    for s0 = cand(o(1:min(20, end)))'
        D = distances(G, s0);
        ends = cand(D(cand) >= L & D(cand) <= 1.2 * L);
        if isempty(ends), continue; end
        v = px(ends, :) - px(s0, :);
        flat = vecnorm(v, 2, 2) / (L * cam.s);                 % chord on the image / its length: 1 = flat to the view
        ends = ends(flat >= 0.75);  v = v(flat >= 0.75, :);     % not foreshortened into a sulcus
        if isempty(ends), continue; end
        slope = abs(v(:, 2)) ./ max(vecnorm(v, 2, 2), eps) + 0.5 * (v(:, 1) > 0);   % horizontal, leftward first
        [~, eo] = sort(slope);
        for e = ends(eo(1:min(5, end)))'
            path = shortestpath(G, s0, e);
            w = vecnorm(diff(Sf.Vertices(path, :)), 2, 2);
            cw = [0; cumsum(w)];
            j = find(cw >= L, 1);
            f = (L - cw(j-1)) / w(j-1);
            P = [Vd(path(1:j-1), :); Vd(path(j-1), :) + f * (Vd(path(j), :) - Vd(path(j-1), :))];
            if all(visible(P, cam.d, Vd, F(all(shown(F), 2), :)))
                bar = struct('mm', Lmm, 'length_mm', 1e3 * (cw(j-1) + f * w(j-1)), 'vertices', path(1:j), 'px', project(cam, P));
                return
            end
        end
    end
end

%% ---------------------------------------------------------------- the panel

function [r, G, job] = compose(s, surfs, F, hemi, folded, G, tok, cmaps, o)
% the panel as an SVG in mm at print size (views embedded as PNG at their rendered resolution; text, bars and colourbar
% vector, editable), and the print job that turns it into PDF and PNG
    J = tok.journal.(tok.journal.default);
    ink = tok.colour.neutral.neutral_ink.x_value;  slate = tok.colour.neutral.neutral_slate.x_value;
    paper = tok.colour.neutral.neutral_paper.x_value;
    fam = strjoin(strcat('''', tok.type.family.x_value, ''''), ', ');
    nview = numel(s.renders);  strip = numel(s.frames) > 1;
    W = panelWidth(o.Width, nview, J);

    % the renders, cropped to what they show; one scale (mm per pixel) for the row
    img = cell(1, nview);  al = img;  crop = zeros(nview, 4);
    for i = 1:nview
        [img{i}, ~, al{i}] = imread(fullfile(s.dir, [s.renders(i).label '.png']));
        if isempty(al{i}), error('rheome:snapshot:alpha', 'The render has no alpha (background must be transparent).'); end
        [rr, cc] = find(al{i} > 0);
        pad = round(0.02 * size(al{i}, 2));
        crop(i, :) = [max(1, min(cc) - pad), min(size(al{i}, 2), max(cc) + pad), max(1, min(rr) - pad), min(size(al{i}, 1), max(rr) + pad)];
    end
    if strip, crop = repmat([min(crop(:, 1)) max(crop(:, 2)) min(crop(:, 3)) max(crop(:, 4))], nview, 1); end   % one box for every frame
    gap = tok.spacing.viewGap_mm;
    cw = (W - (nview - 1) * gap) / nview;
    k = min(cw ./ (crop(:, 2) - crop(:, 1) + 1));                 % mm per pixel
    rowH = k * max(crop(:, 4) - crop(:, 3) + 1);
    dpi = 25.4 / k;
    if dpi < J.dpi.combination, warning('rheome:snapshot:dpi', '%s: the views print at %.0f dpi (< %d); raise Scale or narrow Width.', s.Name, dpi, J.dpi.combination); end

    % the vertical layout, mm from the top
    pt = 25.4 / 72;
    fa = tok.type.annotation.size_pt;  fv = tok.type.viewLabel.size_pt;
    y0 = 0;  if strip, y0 = 1.8 * fa * pt; end
    y1 = y0 + rowH;
    below = 1.6 * fv * pt;  if strip, below = 2.4 * fa * pt + 1; end
    yBar = y1 + below + tok.spacing.colorbarGap_mm;  th = tok.spacing.colorbarThickness_mm;
    H = yBar + th + 0.6 + 1.25 * (tok.type.colorbarTick.size_pt + tok.type.colorbarLabel.size_pt) * pt + 0.4;

    svg = {sprintf(['<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%gmm" height="%gmm" ' ...
        'viewBox="0 0 %g %g" font-family="%s">'], W, H, W, H, fam), ...
        '<style>@page{size:auto;margin:0}</style>', sprintf('<rect width="%g" height="%g" fill="%s"/>', W, H, paper)};
    T = @(x, y, str, sz, col, anchor, extra) sprintf('<text x="%.3f" y="%.3f" font-size="%.4f" fill="%s" text-anchor="%s"%s>%s</text>', ...
        x, y, sz * pt, col, anchor, extra, esc(str));

    r = struct('name', s.Name, 'svg', '', 'pdf', '', 'png', '', 'renders', {fullfile(s.dir, strcat({s.renders.label}, '.png'))}, ...
        'clim', s.CLim, 'colormap', s.Colormap, 'quantity', s.Quantity, 'colorbar', '', 'views', {unique({s.renders.view}, 'stable')}, ...
        'times', s.times, 'dpi', dpi, 'scalebar', [], 'legend', '', 'settings', {{s.renders.settings}});
    Vd = surfs.(s.Surface).Vertices;
    for i = 1:nview
        c = crop(i, :);  wpx = c(2) - c(1) + 1;  hpx = c(4) - c(3) + 1;
        x = (i - 1) * (cw + gap) + (cw - k * wpx) / 2;  y = y0 + (rowH - k * hpx) / 2;
        svg{end+1} = sprintf('<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" preserveAspectRatio="none" xlink:href="data:image/png;base64,%s"/>', ...
            x, y, k * wpx, k * hpx, png64(img{i}(c(3):c(4), c(1):c(2), :), al{i}(c(3):c(4), c(1):c(2)))); %#ok<AGROW>
        cx = x + k * wpx / 2;
        if strip
            svg{end+1} = T(cx, y0 - 0.25 * fa * pt, sprintf('t = %s', timeText(s.times(i))), fa, ink, 'middle', ''); %#ok<AGROW>
        else
            svg{end+1} = T(cx, y1 + 1.1 * fv * pt, s.renders(i).view, fv, slate, 'middle', ''); %#ok<AGROW>
        end
        % one geodesic scale bar per row: on the first view where a visible stretch is found
        if o.ScaleBar > 0 && isempty(r.scalebar)
            h = s.renders(i).hemisphere;  shown = hemi > 0;
            if ~strcmp(h, 'both'), shown = hemi == find(strcmp({'left', 'right'}, h)); end
            cam = fitCamera(al{i}, Vd(shown, :), s.renders(i).camera);
            [bar, G] = scaleBar(o.ScaleBar, surfs.(folded), Vd, F, shown, cam, al{i}, G);
            if ~isempty(bar)
                q = [x y] + k * (bar.px - [c(1) c(3)] + 0.5);          % pixel centres -> mm on the page
                pts = strjoin(ptsText(q), ' ');
                wl = tok.scalebar.space.width_pt * pt;
                svg{end+1} = sprintf(['<polyline points="%s" fill="none" stroke="%s" stroke-width="%.4f" stroke-linejoin="round" stroke-linecap="round"/>' ...
                    '<polyline points="%s" fill="none" stroke="%s" stroke-width="%.4f" stroke-linejoin="round" stroke-linecap="round"/>'], ...
                    pts, paper, wl + 2 * 0.8 * pt, pts, ink, wl); %#ok<AGROW>   % ink over a paper halo
                lab = sprintf('%g mm', o.ScaleBar);  lxy = [mean(q([1 end], 1)), min(q(:, 2)) - 0.7];
                for ang = (0:7) * pi / 4   % the halo: eight paper copies 0.8 pt out (a stroked text prints as a Type 3 font)
                    svg{end+1} = T(lxy(1) + 0.8 * pt * cos(ang), lxy(2) + 0.8 * pt * sin(ang), lab, fa, paper, 'middle', ''); %#ok<AGROW>
                end
                svg{end+1} = T(lxy(1), lxy(2), lab, fa, ink, 'middle', ''); %#ok<AGROW>
                bar.view = s.renders(i).view;
                r.scalebar = bar;
            end
        end
    end
    if o.ScaleBar > 0 && isempty(r.scalebar), warning('rheome:snapshot:scalebar', '%s: no visible stretch of %g mm for the scale bar; none drawn.', s.Name, o.ScaleBar); end
    if strip
        ya = y1 + 1.0;  hl = 1.2;  x0 = 0.05 * W;  xe = 0.95 * W;
        svg{end+1} = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="%s" stroke-width="%.4f"/><polygon points="%.3f,%.3f %.3f,%.3f %.3f,%.3f" fill="%s"/>', ...
            x0, ya, xe - hl, ya, ink, tok.line.axis_pt * pt, xe, ya, xe - hl, ya - 0.45 * hl, xe - hl, ya + 0.45 * hl, ink);
        svg{end+1} = T(W / 2, ya + 0.6 + fa * pt, sprintf('%st = %s', char(916), timeText(s.times(2) - s.times(1))), fa, ink, 'middle', '');
    end

    % the colourbar: the viewer's own colours, quantity (unit re reference), both ends ticked, 0 on a diverging bar
    cmap = cmaps(s.Colormap);
    len = tok.spacing.colorbarLength_frac * W;  xb = (W - len) / 2;
    svg{end+1} = sprintf('<image x="%.3f" y="%.3f" width="%.3f" height="%.3f" preserveAspectRatio="none" xlink:href="data:image/png;base64,%s"/>', ...
        xb, yBar, len, th, png64(uint8(round(255 * reshape(cmap, 1, [], 3))), []));
    svg{end+1} = sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="none" stroke="%s" stroke-width="%.4f"/>', xb, yBar, len, th, slate, tok.line.outline_pt * pt);
    switch s.Quantity
        case {'signed', 'coupling'}, ticks = [s.CLim(1) 0 s.CLim(2)];  labels = arrayfun(@num, ticks, 'UniformOutput', false);
        case 'phase', ticks = [-pi 0 pi];  labels = {[char(8722) char(960)], '0', char(960)};
        otherwise, ticks = [s.CLim(1) mean(s.CLim) s.CLim(2)];  labels = arrayfun(@num, ticks, 'UniformOutput', false);
    end
    tl = tok.spacing.tickLength_pt * pt;  ft = tok.type.colorbarTick.size_pt;
    for j = 1:numel(ticks)
        tx = xb + len * (ticks(j) - s.CLim(1)) / diff(s.CLim);
        svg{end+1} = sprintf('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="%s" stroke-width="%.4f"/>', tx, yBar + th, tx, yBar + th + tl, ink, tok.line.tick_pt * pt); %#ok<AGROW>
        svg{end+1} = T(tx, yBar + th + tl + 0.3 + 0.8 * ft * pt, labels{j}, ft, ink, 'middle', ''); %#ok<AGROW>
    end
    if isempty(s.Reference), r.colorbar = sprintf('%s (%s)', s.Label, s.Unit); else, r.colorbar = sprintf('%s (%s re %s)', s.Label, s.Unit, s.Reference); end
    svg{end+1} = T(W / 2, H - 0.4 - 0.2 * tok.type.colorbarLabel.size_pt * pt, r.colorbar, tok.type.colorbarLabel.size_pt, ink, 'middle', '');
    svg{end+1} = '</svg>';

    % what the legend must say
    sname = s.Surface;
    if strcmp(s.Surface, folded)
        if strcmp(sname, 'default'), sname = 'folded cortical'; else, sname = sprintf('folded (%s)', sname); end
    end
    lg = sprintf('%s view%s of the %s surface', strjoin(r.views, ', '), plural(numel(r.views)), sname);
    if ~isempty(r.scalebar)
        if strcmp(s.Surface, folded), lg = sprintf('%s; scale bar %g mm along the folded surface', lg, o.ScaleBar);
        else, lg = sprintf('%s; scale bar %g mm, distances on the folded surface', lg, o.ScaleBar);
        end
    end
    if strip, lg = sprintf('%s; %d frames, %st = %s', lg, numel(s.frames), char(916), timeText(s.times(2) - s.times(1))); end
    r.legend = [lg '.'];

    base = fullfile(o.Out, s.Name);
    r.svg = [base '.svg'];
    fid = fopen(r.svg, 'w');
    if fid < 0, error('rheome:snapshot:write', 'Cannot write %s.', r.svg); end
    fwrite(fid, unicode2native(strjoin(svg, newline), 'UTF-8'));  fclose(fid);
    if any(o.Format == "pdf"), r.pdf = [base '.pdf']; end
    if any(o.Format == "png"), r.png = [base '.png']; end
    job = struct('svg', r.svg, 'pdf', r.pdf, 'png', r.png, 'width_mm', W, 'height_mm', H, 'dpi', J.dpi.combination);
end

function W = panelWidth(W, nview, J)
% mm: given, or the journal's double column for a row of views, single for one
    if isempty(W), if nview > 1, W = 'double'; else, W = 'single'; end, end
    if ~isnumeric(W), W = J.width_mm.(pickName(W, fieldnames(J.width_mm), 'Width')); end
end

function c = ptsText(q)
    c = arrayfun(@(i) sprintf('%.3f,%.3f', q(i, 1), q(i, 2)), 1:size(q, 1), 'UniformOutput', false);
end

function b = png64(rgb, alpha)
% an image as base64 PNG (alpha kept when given)
    f = [tempname '.png'];  c = onCleanup(@() delete(f));
    if isempty(alpha), imwrite(rgb, f); else, imwrite(rgb, f, 'Alpha', alpha); end
    fid = fopen(f, 'r');  b = matlab.net.base64encode(fread(fid, inf, '*uint8')');  fclose(fid);
end

function t = esc(t)
    t = strrep(strrep(strrep(t, '&', '&amp;'), '<', '&lt;'), '>', '&gt;');
end

function t = timeText(s)
    if abs(s) < 1, t = sprintf('%g ms', round(s * 1e3, 3, 'significant')); else, t = sprintf('%g s', round(s, 4, 'significant')); end
end

function t = num(x)
    t = strrep(sprintf('%.3g', x), '-', char(8722));
end

function s = plural(n)
    s = '';  if n > 1, s = 's'; end
end
