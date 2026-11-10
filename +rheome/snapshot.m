function out = snapshot(surf, result, varargin)
% SNAPSHOT  Print-ready cortical panels from the Cortical Flow viewer, headless, by the lab's figure tokens.
%
%   out = rheome.snapshot(surf, result, Shots=shots, Out=folder)
%   out = rheome.snapshot(surf, result, Shots=struct('Field', "div", 'Unit', "1/s", 'Label', "divergence"), Out="figs")
%
% Every cortical panel of a figure the same way: the Cortical Flow viewer (nxr-cortical-flow editor/dist-viewer, the
% renderer of rheome.view) draws AND COMPOSES each panel in headless Chrome, by the figure tokens
% (library/styles/tokens/figure-tokens.json): the colormap and limits by quantity, the colourbar with quantity, unit
% and reference, view labels, a geodesic scale bar on the folded surface, selections, paths, an inset, and for a time
% series one strip with the interval stated once. The panel is the viewer's print figure (decision D4 of nxr's
% docs/specs/2026-11-09-nxr-viewer-twin-design.md): THE KEY BAND ON TOP (the colourbar, then a key of every mark: the
% path, its onset, its ticks, its hidden stretches, the selections, the bar), then the views. rheome.snapshot only
% checks the shots, writes the bundle and the job, runs the viewer's own driver (<viewer>/driver/snapshot-driver.mjs,
% job format nxr-snapshot-job@1) and reads back its report; nxr.snapshot (Python) sends the same job byte for byte.
% Each shot is written as <Out>/<Name>.svg (text, bars and colourbar editable; the views embedded at their rendered
% resolution), printed by Chrome to <Out>/<Name>.pdf (fonts embedded, the views at full resolution) and .png (at the
% journal's combination dpi), with <Out>/<Name>.marks.json (every overlay drawn and its key: figure_check.py E1,
% rheome.viz.marks) and <Out>/<Name>.report.json (the viewer's report, 0-based); the viewer's raw renders are kept as
% <Out>/<Name>/<view>.png (a strip's frames <view>_f<frame>.png, the frame 0-based; an inset's render
% <view>_inset.png), for boards composed elsewhere. <Out>/job.json is the job sent to the viewer.
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
%   Colormap   a viewer colormap name (or a figure-token colormap), in place of the quantity's token
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
%   Inset      a MULTISCALE INSET (default none): a region of one view re-rendered by the viewer with its camera on
%              the region (and, if given, with a finer field), cut to the region's shape on the cortex and laid out
%              in its own column right of the views, the region outlined on the view and around the inset and
%              joined by a leader, "xm" under it. Not a stretched crop: it prints at >= the journal's dpi.
%              A struct with
%                Vertices  the region: indices or a logical [nV x 1] (any shape: a tile, a parcel, a patch)
%                          (default: the shot's first Selection)
%                Label     what the inset shows, under it ("tile 12, level 5")                        REQUIRED
%                Field     the real field painted in the inset (default the shot's Field), e.g. the field at a
%                          deeper level; it shares the panel's colourbar (same quantity, unit and limits)
%                Zoom      the inset's pixel ratio over the views' (default 4, >= 2: rule I1 of the tokens' inset;
%                          the viewer refuses an inset under 2x the view's pixels per mm)
%                View      the view the region is taken from (default the first of the shot's views)
%                Shape     "tile" (default): cut by the region's own mask and outlined by its contour (a dyadic tile,
%                          a parcel); "disk": the same, for a geodesic disk; "square": a square around the region
%                          (tokens annotation.inset). With no Vertices, the shot's first Selection is the region
%                          (and a disk selection makes a disk inset)
%              A strip takes no inset.
%   Selection  SELECTED regions (default none), each a state of the surface (tokens cortex.selection): a blended
%              fill in its axis colour at the tokens' alpha, outlined as a tile; or, Dashed, an EXTENT (a support
%              or mask: a dashed ink line, no fill). A struct (array) with
%                Vertices  the region: indices or a logical [nV x 1] (a tile, a parcel)          one of these two
%                Centre, RadiusMM   a GEODESIC DISK: every vertex within RadiusMM of Centre along the folded surface
%                Label     its name for the key ("tile 13 (level 4), selected")                       REQUIRED
%                Axis      "space" (default) | "time" | "joint": the colour (atlas-space, atlas-time, atlas-joint-cell)
%                Dashed    true for an extent (default false)
%              A Centre (also allowed with Vertices) is drawn as the tokens' centre: one size, filled with the colour.
%   Paths      GEODESIC paths on the folded surface (default none; tokens cortex.paths): ink over a paper halo, an
%              arrowhead at the end, its stretches hidden by the surface thin and dotted. A struct (array) with
%                From, To  vertices: the viewer's exact geodesic between them (geometry-central) on the folded surface
%                Label     what it is ("path of the event")                                            REQUIRED
%              (Points with How, a path the caller traced, is refused for now: the viewer draws only the geodesics
%              between From and To.)
%
% OPTIONS (name/value):
%   'Out'       folder written (created)                                                     REQUIRED
%   'Fs'        sampling rate of result's columns, Hz (default 1)
%   'Width'     the panel's width: "single" | "onehalf" | "double" (the token theme's widths) or mm, from the theme
%               page's least width to its double column (default: "double" for several views, a strip or an inset,
%               "single" for one view)
%   'ScaleBar'  the scale bar's length, mm along the folded surface; 0: none (default 20)
%   'Folded'    the name of the folded surface in surf, on which distances are measured (default "pial" when surf has
%               it, else the first); on another surface shown, the legend says the distance is on the folded one
%   'Format'    printed from the SVG (always written): "pdf", "png", both (default) or [] (the SVG only)
%   'Tokens'    the figure tokens' JSON (default $NXR_FIGURE_TOKENS, else $HOME/workspace/library/styles/tokens/
%               figure-tokens.json if it exists, else the viewer's vendored viewer/figure/figure-tokens.json). The viewer
%               is served the file used, and its sha256 and path go into every <Name>.marks.json (tokens_sha256, tokens_path).
%   'Viewer'    the viewer folder (default as rheome.view: $RHEOME_VIEWER, else prefdir/rheome-viewer); it must ship
%               driver/snapshot-driver.mjs (a viewer build of nxr-cortical-flow ad9d824, 2026-11-09, or later)
%   'Chrome'    Chrome or Chromium (default $RHEOME_CHROME, else the usual install)
%   'Node'      Node.js 22+ (default $RHEOME_NODE, else the usual install)
%   'Labellings' partitions for the viewer (rheome.viz.tiles), shown by a shot's Settings.labelling
%   'Points'    point sets for the viewer (rheome.viz.points: critical points), shown by Settings.points (default all)
%   'Wait'      ms after the viewer opens the bundle (default 600); 'Timeout' s per shot (default 60)
%   'Scale'     ignored (kept so older calls run): the viewer raises its pixel ratio itself so the views print at the
%               journal's combination dpi
%
% OUTPUT:
%   out : struct array, one per shot: .name .svg .pdf .png .renders (the viewer's raw PNGs of the views, without the
%         inset's) and the viewer's report, its indices made 1-based: .clim .colormap .quantity .colorbar (its label)
%         .views .times (s) .dpi (resolution of the views as printed) .scalebar (.mm, .length_mm measured along the
%         drawn curve, .vertices, .px, .view) or [] .legend (the sentence the figure legend needs: views, surface, the
%         bar's distance, the interval, the inset, the marks) .inset (.label .view .magnification (x, the inset's
%         length per the view's) .dpi .vertices .shape) or [] .selections (.label .kind ("tile" | "disk") .axis .colour
%         .dashed .vertices) .paths (.label .how .points [n x 3] m .length_mm along the folded surface, .geodesic_mm
%         the distance between its ends, .px and .visible: one cell per render) .build (the viewer's build)
%   PIXELS (.scalebar.px, .paths.px) are in the raw render's pixels (<Out>/<Name>/<view>.png): x right, y down, the
%   origin at the render's top-left CORNER (pixel (r, c) of imread covers [c-1, c] x [r-1, r]).
%
% THE SCALE BAR is a curve on the surface, not a line on the image: the viewer's exact geodesic on the folded surface,
% cut at exactly 'ScaleBar' mm of its length, placed by Rheome's rule (an anchor among the front-most vertices near
% the lower right of the silhouette, ends 1-1.2 L away, flat to the view, horizontal and leftward first) where every
% point of it is seen; one bar a row. Limits, visibility (a ray cast toward the camera) and the bar are the viewer's
% (nxr-cortical-flow editor/src/viewer/print/limits.ts, visible.ts, scalebar.ts), by the rules Rheome used.
%
% Needs the viewer (WebGPU), Chrome and Node; nothing else in Rheome does.
%
% Example:
%   S = rheome.load.surface('omega_sub0002');  x = S.Vertices(:, 1) * 1e3;
%   out = rheome.snapshot(S, struct('x', x), Out=tempname, ...
%       Shots=struct('Field', "x", 'Label', "anterior position", 'Unit', "mm"));
%   open(out.pdf)
%
% See also: rheome.view, rheome.io.write.bundle, rheome.viz.marks
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
    T = size(result.(firstOf(fieldnames(result))), 2);

    work = tempname;  mkdir(work);  cleanWork = onCleanup(@() rmdir(work, 's'));
    [tok, tokFile, tokSha, tokShown] = tokens(o.Tokens, work);
    vdir = char(o.Viewer);
    if isempty(vdir), vdir = viewerFolder(); end
    if ~isfile(fullfile(vdir, 'viewer.html')), error('rheome:snapshot:noviewer', 'No viewer at %s (Viewer= or RHEOME_VIEWER).', vdir); end
    drv = fullfile(vdir, 'driver', 'snapshot-driver.mjs');
    if ~isfile(drv)
        error('rheome:snapshot:driver', ['The viewer at %s has no snapshot driver (driver/snapshot-driver.mjs): rheome.snapshot needs ' ...
            'a Cortical Flow viewer build of nxr-cortical-flow ad9d824 (2026-11-09) or later.'], vdir);
    end
    cmaps = viewerColormaps(vdir);

    folded = o.Folded;
    if isempty(folded), folded = snames{1}; if any(strcmp(snames, 'pial')), folded = 'pial'; end, end
    folded = pickName(folded, snames, 'Folded');
    J = tok.theme.(tok.theme.default).page;
    shots = normalise(o.Shots, result, surfs, T, o.Fs, tok, cmaps);

    % one bundle, one job, one viewer run for every shot
    if isfield(surf, 'Vertices'), bs = surf; else, bs = surfs; end
    rheome.io.write.bundle(fullfile(work, 'bundle'), bs, result, 'Fs', o.Fs, 'Name', 'rheome.snapshot', ...
        'Labellings', o.Labellings, 'Points', o.Points);
    list = cell(1, numel(shots));
    for k = 1:numel(shots), list{k} = struct('name', shots(k).Name, 'shot', viewerShot(shots(k), folded, o, J)); end
    job = struct();
    job.format = 'nxr-snapshot-job@1';  job.viewer = vdir;  job.bundle = 'bundle';  job.out = 'out';
    job.chrome = chromePath(o.Chrome);  job.tokens = tokFile;  job.tokens_sha256 = tokSha;  job.tokens_path = tokShown;
    job.size = double(tok.cortex.viewer.size_px(:)');  job.wait = o.Wait;  job.timeout = o.Timeout;  job.shots = list;
    node = nodePath(o.Node);
    jf = runDriver(job, 'job.json', node, drv, work);
    if ~isfolder(o.Out), mkdir(o.Out); end
    copyfile(jf, fullfile(o.Out, 'job.json'));

    % the files into Out, and out from each shot's report
    done = jsondecode(fileread(fullfile(work, 'out', 'driver.json')));
    ds = done.shots;  if isstruct(ds), ds = num2cell(ds); end
    out = struct([]);  prints = {};
    for k = 1:numel(shots)
        s = shots(k);  d = ds{k};  from = fullfile(work, 'out');
        labels = cellstr(d.renders);  labels = labels(:)';
        if ~isempty(d.warnings)
            for w = cellstr(d.warnings)', warnOf(s.Name, w{1}); end
        end
        rdir = fullfile(o.Out, s.Name);
        if ~isfolder(rdir), mkdir(rdir); end
        for j = 1:numel(labels), movefile(fullfile(from, s.Name, [labels{j} '.png']), fullfile(rdir, [labels{j} '.png']), 'f'); end
        base = fullfile(o.Out, s.Name);
        for ext = {'.svg', '.marks.json', '.report.json'}, movefile(fullfile(from, [s.Name ext{1}]), [base ext{1}], 'f'); end
        r = fromReport(jsondecode(fileread([base '.report.json'])), s, base, fullfile(rdir, strcat(labels(1:s.ntiles), '.png')), o.Format);
        out = [out r]; %#ok<AGROW>
        if ~isempty(r.pdf) || ~isempty(r.png)
            [wmm, hmm] = svgSize(r.svg);
            prints{end+1} = struct('svg', r.svg, 'pdf', r.pdf, 'png', r.png, 'width_mm', wmm, 'height_mm', hmm, 'dpi', J.dpi.combination); %#ok<AGROW>
        end
    end
    if ~isempty(prints)   % the SVGs printed by Chrome: MATLAB's PDF export resamples images to screen pixels
        runDriver(struct('chrome', job.chrome, 'print', {prints}), 'print.json', node, drv, work);
    end
end

%% ---------------------------------------------------------------- the shots

function shots = normalise(in, result, surfs, T, fs, tok, cmaps)
% each shot checked, as the viewer will draw it (it is refused here, by name, rather than by the viewer)
    if isempty(in), in = {struct()}; end
    if isstruct(in), in = num2cell(in); end
    if ~iscell(in), error('rheome:snapshot:input', 'Shots: a struct array or a cell array of structs is expected.'); end
    known = {'Field', 'Unit', 'Label', 'Reference', 'Quantity', 'CLim', 'Colormap', 'Time', 'Frames', 'View', ...
        'Hemisphere', 'Surface', 'Vector', 'Settings', 'Name', 'Inset', 'Selection', 'Paths'};
    fields = fieldnames(result);
    S1 = surfs.(firstOf(fieldnames(surfs)));
    nV = size(S1.Vertices, 1);
    real = fields(cellfun(@(f) size(result.(f), 1) == nV, fields));
    vec = fields(cellfun(@(f) size(result.(f), 1) == 3*nV, fields));
    hemi = [];
    for k = numel(in):-1:1
        a = in{k};
        if ~isstruct(a) || ~isscalar(a), error('rheome:snapshot:input', 'Shots{%d}: a struct is expected.', k); end
        bad = setdiff(fieldnames(a), known);
        if ~isempty(bad), error('rheome:snapshot:shot', 'Shots{%d}: %s is not a shot field (%s).', k, strjoin(bad, ', '), strjoin(known, ', ')); end
        g = @(n, d) get(a, n, d);
        s = struct();
        s.Name = char(g('Name', sprintf('shot%d', k)));
        if isempty(regexp(s.Name, '^[\w\-\.]+$', 'once')), error('rheome:snapshot:shot', 'Name: letters, digits, - _ . only (%s).', s.Name); end
        if isempty(real), error('rheome:snapshot:shot', 'The result has no real ([nV x T]) field to paint.'); end
        s.Field = pickName(g('Field', real{1}), real, 'Field');
        s.Unit = char(g('Unit', ''));
        if isempty(s.Unit), error('rheome:snapshot:unit', '%s: a colourbar needs a unit (Unit="a.u." only when there is none; say why in the legend).', s.Name); end
        s.Label = char(g('Label', s.Field));
        s.Reference = char(g('Reference', ''));

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
        s.frames = fr;

        % views: the figure system's four, one of them, or a viewer camera with a hemisphere
        v = char(g('View', 'standard'));
        views4 = {'L lateral', 'left', 'left'; 'L medial', 'left', 'right'; 'R medial', 'right', 'left'; 'R lateral', 'right', 'right'};
        s.Hemisphere = '';
        if strcmpi(v, 'standard'), s.views = views4;  s.view = 'standard';
        elseif any(strcmpi(v, views4(:, 1))), s.views = views4(strcmpi(v, views4(:, 1)), :);  s.view = s.views{1, 1};
        else
            cam = pickName(v, {'left', 'right', 'top', 'bottom', 'front', 'back'}, 'View');
            h = pickName(g('Hemisphere', 'both'), {'both', 'left', 'right'}, 'Hemisphere');
            s.views = {sprintf('%s, %s', h, cam), h, cam};  s.view = cam;  s.Hemisphere = h;
        end
        if numel(fr) > 1 && size(s.views, 1) > 1, error('rheome:snapshot:shot', '%s: a strip takes one view (one view and one scale for the whole strip).', s.Name); end
        s.Surface = pickName(g('Surface', firstOf(fieldnames(surfs))), fieldnames(surfs), 'Surface');

        % the colour scale checked: the quantity's colormap is the viewer's, a diverging CLim symmetric. The class when
        % none is given is the viewer's rule (both signs over the frames and hemispheres shown: signed), for the checks
        % only: the viewer is sent what was given.
        q = char(g('Quantity', ''));
        classes = {'signed', 'magnitude', 'power', 'phase', 'coupling', 'code'};
        if isempty(q)
            if isempty(hemi), hemi = hemispheres(S1); end
            shown = ismember(hemi, find(ismember({'left', 'right'}, s.views(:, 2)))) | any(strcmp(s.views(:, 2), 'both'));
            x = result.(s.Field);  xs = x(shown, fr);  xs = xs(isfinite(xs));
            if any(xs < 0) && any(xs > 0), cls = 'signed'; else, cls = 'magnitude'; end
            s.Quantity = '';
        else
            cls = pickName(q, classes, 'Quantity');  s.Quantity = cls;
        end
        qn = cls;  if strcmp(qn, 'signed'), qn = 'divergence'; end   % the tokens name each signed quantity; one map
        s.Colormap = char(g('Colormap', ''));
        cm = s.Colormap;
        if isempty(cm), cm = regexp(tok.quantity.(qn).colormap, '^[\w\-]+', 'match', 'once'); end
        tk = strrep(cm, '-', '_');   % a figure-token name stands for its colormap, as the viewer reads it: $value, 'cmo.' dropped
        if isfield(tok.colour.colormap, tk), cm = regexprep(char(tok.colour.colormap.(tk).x_value), '^cmo\.', ''); end
        if ~isKey(cmaps, cm)
            error('rheome:snapshot:colormap', '%s: the viewer has no colormap ''%s'' (%s quantity); give Colormap= or add it to the viewer.', s.Name, cm, cls);
        end
        c = g('CLim', []);
        if ~isempty(c)
            c = double(c(:)');
            if numel(c) ~= 2 || ~all(isfinite(c)) || c(1) >= c(2), error('rheome:snapshot:shot', '%s: CLim [lo hi], lo < hi.', s.Name); end
            if any(strcmp(cls, {'signed', 'coupling'})) && abs(c(1) + c(2)) > 1e-9 * max(abs(c))
                error('rheome:snapshot:clim', '%s: a diverging scale is symmetric around 0 (got [%g %g]).', s.Name, c(1), c(2));
            end
        end
        s.CLim = c;
        s.Vector = char(g('Vector', ''));
        if ~isempty(s.Vector), s.Vector = pickName(s.Vector, vec, 'Vector'); end
        s.Settings = g('Settings', struct());
        if ~isstruct(s.Settings), error('rheome:snapshot:shot', '%s: Settings is a struct of viewer settings.', s.Name); end
        s.Selection = selections(g('Selection', []), s.Name, nV);
        s.Paths = paths(g('Paths', []), s.Name, nV);
        s.Inset = inset(g('Inset', []), s, real, nV);
        s.ntiles = size(s.views, 1) * numel(fr);
        shots(k) = s;
    end
end

function in = inset(a, s, real, nV)
% a shot's inset, checked, as the viewer's Shot.inset (0-based vertices; what is not given is left to the viewer)
    in = [];
    if isempty(a), return; end
    id = 'rheome:snapshot:inset';
    if ~isstruct(a) || ~isscalar(a), error(id, '%s: Inset is a struct.', s.Name); end
    bad = setdiff(fieldnames(a), {'Vertices', 'Label', 'Field', 'Zoom', 'View', 'Shape'});
    if ~isempty(bad), error(id, '%s: %s is not an inset field.', s.Name, strjoin(bad, ', ')); end
    if numel(s.frames) > 1, error(id, '%s: a strip takes no inset.', s.Name); end
    in = struct();
    v = get(a, 'Vertices', []);
    if ~isempty(v)
        m = vertexMask(v, nV, sprintf('%s: Inset.Vertices', s.Name), id);
        if nnz(m) < 3, error(id, '%s: Inset.Vertices: a region of 3 vertices or more.', s.Name); end
        in.vertices = num2cell(find(m)' - 1);
    elseif isempty(s.Selection)
        error(id, '%s: Inset.Vertices: the region (indices in 1..%d or a logical [nV x 1]), or a Selection (the first is the region).', s.Name, nV);
    end
    shape = get(a, 'Shape', '');
    if ~isempty(shape), shape = pickName(shape, {'tile', 'disk', 'square'}, 'Inset.Shape'); end
    lab = char(get(a, 'Label', ''));
    if isempty(lab), error(id, '%s: Inset.Label: say what the inset shows (it is drawn under it).', s.Name); end
    in.label = lab;
    if isfield(a, 'Zoom') && ~isempty(a.Zoom)
        z = double(a.Zoom);
        if ~(isscalar(z) && z >= 2), error(id, '%s: Inset.Zoom >= 2 (an inset shows more: tokens annotation.inset, I1).', s.Name); end
        in.zoom = z;
    end
    if isfield(a, 'View') && ~isempty(a.View)
        vi = find(strcmpi(char(a.View), s.views(:, 1)), 1);
        if isempty(vi), error(id, '%s: Inset.View: one of %s.', s.Name, strjoin(s.views(:, 1), ', ')); end
        in.view = s.views{vi, 1};
    end
    if isfield(a, 'Field') && ~isempty(a.Field), in.field = pickName(a.Field, real, 'Inset.Field'); end
    if ~isempty(shape), in.shape = shape; end
end

function m = vertexMask(v, nV, what, id)
    if islogical(v) && numel(v) == nV, m = v(:);
    elseif isnumeric(v) && ~isempty(v) && all(v(:) >= 1 & v(:) <= nV & v(:) == round(v(:))), m = false(nV, 1);  m(v) = true;
    else, error(id, '%s: indices in 1..%d or a logical [nV x 1].', what, nV);
    end
end

function sel = selections(a, name, nV)
% the shot's selections, checked, as the viewer's Shot.selections (0-based); the disk is measured by the viewer
    sel = {};
    if isempty(a), return; end
    id = 'rheome:snapshot:selection';
    if ~isstruct(a), error(id, '%s: Selection is a struct (array).', name); end
    for j = 1:numel(a)
        b = a(j);
        bad = setdiff(fieldnames(b), {'Vertices', 'Centre', 'RadiusMM', 'Label', 'Axis', 'Dashed'});
        if ~isempty(bad), error(id, '%s: %s is not a selection field.', name, strjoin(bad, ', ')); end
        c = get(b, 'Centre', []);
        if ~isempty(c) && ~(isscalar(c) && c >= 1 && c <= nV && c == round(c)), error(id, '%s: Selection.Centre: a vertex in 1..%d.', name, nV); end
        q = struct();
        if isfield(b, 'Vertices') && ~isempty(b.Vertices)
            q.vertices = num2cell(find(vertexMask(b.Vertices, nV, sprintf('%s: Selection.Vertices', name), id))' - 1);
            if ~isempty(c), q.centre = double(c) - 1; end
        elseif ~isempty(c) && isfield(b, 'RadiusMM') && isscalar(b.RadiusMM) && b.RadiusMM > 0
            q.centre = double(c) - 1;  q.radiusMm = double(b.RadiusMM);
        else
            error(id, '%s: Selection: Vertices, or Centre with RadiusMM (a geodesic disk).', name);
        end
        lab = char(get(b, 'Label', ''));
        if isempty(lab), error(id, '%s: Selection.Label: its name for the key.', name); end
        q.label = lab;
        q.axis = pickName(get(b, 'Axis', 'space'), {'space', 'time', 'joint'}, 'Selection.Axis');
        q.dashed = logical(get(b, 'Dashed', false));
        sel{end+1} = q; %#ok<AGROW>
    end
end

function P = paths(a, name, nV)
% the shot's paths, checked, as the viewer's Shot.paths (0-based): the viewer traces the geodesic between From and To
    P = {};
    if isempty(a), return; end
    id = 'rheome:snapshot:path';
    if ~isstruct(a), error(id, '%s: Paths is a struct (array).', name); end
    for j = 1:numel(a)
        b = a(j);
        bad = setdiff(fieldnames(b), {'From', 'To', 'Points', 'How', 'Label'});
        if ~isempty(bad), error(id, '%s: %s is not a path field.', name, strjoin(bad, ', ')); end
        lab = char(get(b, 'Label', ''));
        if isempty(lab), error(id, '%s: Paths.Label: what the path is, for the key.', name); end
        if isfield(b, 'Points') && ~isempty(b.Points)
            if isempty(char(get(b, 'How', ''))), error(id, '%s: Paths.How: a given path states how it was traced on the surface (tokens cortex.paths).', name); end
            X = double(b.Points);
            if size(X, 2) ~= 3 || size(X, 1) < 2, error(id, '%s: Paths.Points: [n x 3], n >= 2.', name); end
            error(id, ['%s: Paths.Points: a given path is not drawn now: the viewer draws only the geodesic between From and To ' ...
                '(vertices in 1..%d) on the folded surface.'], name, nV);
        end
        fr = get(b, 'From', []);  to = get(b, 'To', []);
        if ~(isscalar(fr) && isscalar(to) && all([fr to] >= 1 & [fr to] <= nV & [fr to] == round([fr to])))
            error(id, '%s: Paths: From and To (vertices in 1..%d), or Points with How.', name, nV);
        end
        q = struct();  q.from = double(fr) - 1;  q.to = double(to) - 1;  q.label = lab;
        P{end+1} = q; %#ok<AGROW>
    end
end

function v = viewerShot(s, folded, o, J)
% the viewer's Shot (nxr-cortical-flow editor/src/viewer/print/panel.ts): 0-based frames and vertices, the width in
% mm; the keys in one order, so nxr.snapshot's job is this one byte for byte
    v = struct();
    v.field = s.Field;  v.unit = s.Unit;  v.label = s.Label;
    if ~isempty(s.Reference), v.reference = s.Reference; end
    if ~isempty(s.Quantity), v.quantity = s.Quantity; end
    if ~isempty(s.CLim), v.clim = s.CLim; end
    if ~isempty(s.Colormap), v.colormap = s.Colormap; end
    v.frames = num2cell(s.frames - 1);
    v.views = {s.view};
    if ~isempty(s.Hemisphere), v.hemisphere = s.Hemisphere; end
    v.surface = s.Surface;  v.folded = folded;
    if ~isempty(s.Vector), v.vector = s.Vector; end
    if ~isempty(fieldnames(s.Settings)), v.settings = s.Settings; end
    if ~isempty(s.Selection), v.selections = s.Selection; end
    if ~isempty(s.Paths), v.paths = s.Paths; end
    if ~isempty(s.Inset), v.inset = s.Inset; end
    v.scaleBarMm = double(o.ScaleBar);
    v.widthMm = panelWidth(o.Width, s.ntiles + ~isempty(s.Inset), J);
end

function r = fromReport(rep, s, base, renders, fmt)
% out: the files, and the viewer's report with its indices made 1-based
    r = struct('name', s.Name, 'svg', [base '.svg'], 'pdf', '', 'png', '', 'renders', {renders(:)'});
    if any(fmt == "pdf"), r.pdf = [base '.pdf']; end
    if any(fmt == "png"), r.png = [base '.png']; end
    r.clim = rep.clim(:)';  r.colormap = rep.colormap;  r.quantity = rep.quantity;  r.colorbar = rep.colorbar;
    r.views = cellstr(rep.views);  r.views = r.views(:)';  r.times = rep.times(:)';  r.dpi = rep.dpi;
    r.scalebar = rep.scalebar;
    if ~isempty(r.scalebar), r.scalebar.vertices = r.scalebar.vertices(:) + 1;  r.scalebar.px = reshape(r.scalebar.px, [], 2); end
    r.legend = rep.legend;
    r.inset = rep.inset;
    if ~isempty(r.inset), r.inset.vertices = r.inset.vertices(:) + 1; end
    sel = struct('label', {}, 'kind', {}, 'axis', {}, 'colour', {}, 'dashed', {}, 'vertices', {});
    for j = 1:numel(rep.selections)
        q = rep.selections(j);  if iscell(rep.selections), q = rep.selections{j}; end
        sel(end+1) = struct('label', q.label, 'kind', q.kind, 'axis', q.axis, 'colour', q.colour, 'dashed', logical(q.dashed), ...
            'vertices', q.vertices(:) + 1); %#ok<AGROW>
    end
    r.selections = sel;
    P = struct('label', {}, 'how', {}, 'points', {}, 'length_mm', {}, 'geodesic_mm', {}, 'px', {}, 'visible', {});
    for j = 1:numel(rep.paths)
        q = rep.paths(j);  if iscell(rep.paths), q = rep.paths{j}; end
        P(end+1) = struct('label', q.label, 'how', q.how, 'points', reshape(q.points, [], 3), 'length_mm', q.length_mm, ...
            'geodesic_mm', q.geodesic_mm, 'px', {perRender(q.px, 2, s.ntiles)}, 'visible', {perRender(q.visible, 1, s.ntiles)}); %#ok<AGROW>
        P(end).visible = cellfun(@logical, P(end).visible, 'UniformOutput', false);
    end
    r.paths = P;
    r.build = rep.build;
end

function c = perRender(x, w, n)
% a report's per-render list (JSON [render][point]([x y])) as one cell per render: [points x w]
    if iscell(x), c = cellfun(@(y) reshape(y, [], w), x(:)', 'UniformOutput', false); return; end
    c = cell(1, n);
    for i = 1:n
        if w == 1, c{i} = reshape(x(i, :), [], 1); else, c{i} = reshape(x(i, :, :), [], w); end
    end
end

function warnOf(name, msg)
% the viewer's warnings, under Rheome's identifiers where they are Rheome's
    id = 'rheome:snapshot:viewer';
    if contains(msg, 'scale bar'), id = 'rheome:snapshot:scalebar'; elseif contains(msg, 'dpi'), id = 'rheome:snapshot:dpi'; end
    warning(id, '%s: %s', name, msg);
end

function [w, h] = svgSize(file)
    t = fileread(file);
    m = regexp(t, '<svg[^>]*?\swidth="([\d.]+)mm"[^>]*?\sheight="([\d.]+)mm"', 'tokens', 'once');
    if isempty(m), error('rheome:snapshot:viewer', 'The panel %s has no size in mm.', file); end
    w = str2double(m{1});  h = str2double(m{2});
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

function W = panelWidth(W, ncol, J)
% mm: given, or the journal's double column for a row of views (or a view and its inset), single for one
    if isempty(W), if ncol > 1, W = 'double'; else, W = 'single'; end, end
    if ~isnumeric(W), W = J.width_mm.(pickName(W, fieldnames(J.width_mm), 'Width')); end
    W = double(W);
end

function h = hemispheres(S)
% 1 left, 2 right, 0 neither: the mesh's two largest connected components, the left the one with larger mean y (as
% the viewer finds them); for the default quantity's check only
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

%% ---------------------------------------------------------------- tokens, colormaps, tools

function [T, file, sha, shown] = tokens(file, work)
% the lab's figure tokens (Graphics): $NXR_FIGURE_TOKENS, else the library's, else the viewer's vendored copy
% (rheome.viz.tokenspath); a struct is written to the work folder so the viewer and the marks sidecar name it
    if isstruct(file)
        T = file;  file = fullfile(work, 'figure-tokens.json');
        fid = fopen(file, 'w');  fwrite(fid, jsonencode(T));  fclose(fid);
        [file, sha] = rheome.viz.tokenspath(file);
        shown = '<struct>';   % the temporary file is deleted with the work folder
        return
    end
    [file, sha] = rheome.viz.tokenspath(file);
    shown = file;   % the viewer is served this file's own bytes
    T = jsondecode(fileread(file));
end

function m = viewerColormaps(vdir)
% the viewer's own colour tables (its default library): a colormap the viewer has not is refused before it runs
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

function jf = runDriver(job, name, node, drv, work)
% the job written (UTF-8 jsonencode: nxr.snapshot writes the same bytes) and the viewer's driver run on it
    jf = fullfile(work, name);
    fid = fopen(jf, 'w');  fwrite(fid, unicode2native(jsonencode(job), 'UTF-8'));  fclose(fid);
    [st, msg] = system(sprintf('"%s" "%s" "%s" 2>&1', node, drv, jf));
    if st ~= 0, error('rheome:snapshot:viewer', 'The viewer run failed: %s', strtrim(msg)); end
end
