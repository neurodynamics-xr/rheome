function out = snapshot(surf, result, varargin)
% SNAPSHOT  Print-ready cortical panels from the Cortical Flow viewer, headless, by the lab's figure tokens.
%
%   out = rheome.snapshot(surf, result, Shots=shots, Out=folder)
%   out = rheome.snapshot(surf, result, Shots=struct('Field', "div", 'Unit', "1/s", 'Label', "divergence"), Out="figs")
%
% Every cortical panel of a figure the same way: the Cortical Flow viewer (the
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
% RAW PICTURES (a shot with StagePx): each render exactly StagePx pixels, kept as <Out>/<Name>/<label>.png with its own
% marks <Out>/<Name>/<label>.marks.json (that render's overlays, pixels in its raster: paths, selection edges, ticks,
% critical points, the bar; every vector mark left out of the PNG for the host to draw), for a host that composes them
% itself; the panel (<Name>.svg, .pdf, .png) only when 'Width' is given. <Name>.marks.json and .report.json always.
% THE PRIMARY API for printed cortex, with rheome.view (interactive); see rheome.viz.
%
% THE SHOT LIST -- Shots is a struct array or a cell array of structs, one panel each, with the fields:
%   Field      the real field painted (a field of result) (default: the first real; none when the shot paints a
%              Labelling or result has no real field). "" (or []): NO FIELD -- a labelling, arrows or the bare
%              cortex: no colour layer, no colourbar (out's clim, colormap, quantity, colorbar empty); Unit, Label,
%              Reference, Quantity, CLim and Colormap are then refused
%   Unit       its unit, as printed on the colourbar ("a.u." only when it has none)        REQUIRED with a Field
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
%              a viewer camera "left" | "right" | "top" | "bottom" | "front" | "back" with Hemisphere, or a CAMERA
%              DIRECTION [x y z] with Hemisphere (from the brain toward the camera, as the viewer's setting view; up the
%              viewer's +z; its label "<hemisphere>, [x, y, z]", its render <hemisphere>-from-<x>_<y>_<z>, the numbers
%              as given with - written m and . written p); or a
%              list of these (a cell or string array: the views in that order, each once) (default "standard"). A
%              strip of several views is a row a view
%   Hemisphere "both" | "left" | "right", with a viewer camera (default "both"; one for every camera of the list)
%   Surface    the surface shown, a name of surf's struct (default: the first)
%   Vector     a vector field drawn as arrows (default none); arrows show the same vector the colour shows
%   Cortex     the cortex surface itself: true (default), false (not drawn: the arrows, outlines and marks alone) or
%              an opacity 0-1
%   Labelling  a labelling of 'Labellings' (its name) painted in its own colours, in place of a field (default none;
%              a Field and a Labelling are refused together)
%   Outlines   TILE OUTLINES drawn by the viewer, keyed in the marks and the panel (default none). A struct with
%                Labelling     the labelling (of 'Labellings') whose tile boundaries are drawn in ink     REQUIRED
%                Coarser       a coarser level's boundaries, thicker: another labelling's name, or n >= 1 (the
%                              same labelling's codes shifted right by n: a dyadic tile's parent is its code halved)
%                Width, CoarserWidth   the lines' widths in pt, 0.05-20 (default 0.35, 0.9)
%                Select        tile codes of the labelling (codes, not indices) marked as selections "tile <code>"
%                SelectStyle   "fill" (default: the selection fill) | "extent" (a dashed outline)
%   StagePx    RAW PICTURES: [w h] integer pixels (16-8192): each render exactly this stage, reported one by one
%              (out.raw) with its marks, for a host that composes them; no inset (default none: the panel)
%   Fit        raw pictures: the camera fitted to the surface shown, a struct with Margin (a fraction of the
%              half-size to spare, 0 to 0.5, default 0.04) and Shared (true: every view fitted first, all take the
%              smallest zoom: one scale), or true for the defaults (default none: the viewer's own framing)
%   Dpi        raw pictures: the dpi they print at (72-1200), which sizes the drawn lines (outlines, isolines) in pt
%              (default the screen's, 96)
%   MarksInRaster  raw pictures: true draws the bundle's paths (BundlePaths), the selections' edges and the points'
%              badges IN each render too (still reported in its marks); default false: left out, for the host to draw
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
%                From, To  vertices: the viewer's geodesic between them on the folded surface (geometry-central's FlipOut:
%                          locally shortest; out.paths reports its length and the exact distance between them)
%                Label     what it is ("path of the event")                                            REQUIRED
%              (Points with How, a path the caller traced, is refused for now: the viewer draws only the geodesics
%              between From and To.)
%
% OPTIONS (name/value):
%   'Out'       folder written (created)                                                     REQUIRED
%   'Fs'        sampling rate of result's columns, Hz (default 1)
%   'Width'     the panel's width: "single" | "onehalf" | "double" (the token theme's widths) or mm, from the theme
%               page's least width to its double column (default: "double" for several views, a strip or an inset,
%               "single" for one view; raw pictures (StagePx): no panel unless given)
%   'ScaleBar'  the scale bar's length, mm along the folded surface; 0: none; "auto": the 1-2-5 round length nearest
%               96 px of the stage (the viewer's on-screen bar), fixed on the render that carries it (default 20)
%   'ScaleBarKind' "geodesic" (default: a path of that length along the folded surface) | "plane": a straight bar of
%               that length in the picture plane (the camera is orthographic: exact at every depth), lower left of
%               its view; out.scalebar.kind is then "plane" and .vertices empty
%   'ScaleBarTile' the render that carries the bar, 1-based in each shot's renders (out.renders) (default: the first
%               where it fits); every shot of the call must have that many renders
%   'Folded'    the name of the folded surface in surf, on which distances are measured (default "pial" when surf has
%               it, else the first); on another surface shown, the legend says the distance is on the folded one
%   'Format'    printed from the SVG (always written): "pdf", "png", both (default) or [] (the SVG only)
%   'Tokens'    the figure tokens' JSON (default $NXR_FIGURE_TOKENS, else $HOME/workspace/library/styles/tokens/
%               figure-tokens.json if it exists, else the viewer's vendored viewer/figure/figure-tokens.json). The viewer
%               is served the file used, and its sha256 and path go into every <Name>.marks.json (tokens_sha256, tokens_path).
%   'Viewer'    the viewer folder (default as rheome.view: $RHEOME_VIEWER, else prefdir/rheome-viewer); it must ship
%               driver/snapshot-driver.mjs (a Cortical Flow viewer build of 2026-11-10 or later)
%   'Chrome'    Chrome or Chromium (default $RHEOME_CHROME, else the usual install)
%   'Node'      Node.js 22+ (default $RHEOME_NODE, else the usual install)
%   'Labellings' partitions for the viewer (rheome.viz.tiles), named by a shot's Labelling and Outlines (or Settings)
%   'Points'    point sets for the viewer (rheome.viz.points: critical points), shown by Settings.points (default all)
%   'Selections' the bundle's own VERTEX SETS (rheome.io.write.bundle's 'Selections': name, vertices, style fill |
%               outline | extent, color, width), shown by a shot's Settings.selections (names; default all). Not a shot's
%               Selection (a key entry the viewer measures): these are drawn as the viewer draws a bundle's sets
%   'BundlePaths' the bundle's own PATHS (rheome.io.write.bundle's 'Paths': name, vertices or faces + weights, color,
%               width, geodesic, arrow, ticks), shown by a shot's Settings.paths (names; default all). Not a shot's Paths
%               (From/To geodesics the report measures). In a raw render (StagePx) each set's edge and each path is an
%               overlay of <Out>/<Name>/<view>.marks.json ('selection: <name>', 'path: <name>') with its pixels and its
%               ink (colour, widthPt, dash, arrow); a fill style's fill stays in the PNG
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
%         .dashed .vertices) .paths (.label .how .points [n x 3] m .length_mm the drawn path's length along the folded
%         surface (geometry-central's FlipOut geodesic, locally shortest), .geodesic_mm the EXACT geodesic distance
%         between its ends (MMP, the shortest path over the faces), so .geodesic_mm <= .length_mm (equal when the drawn
%         path is the shortest; on fsaverage5 pial up to ~8 % apart), .px and .visible: one cell per render) .build
%         (the viewer's build)
%         ONLY WHEN A SHOT OF THE CALL ASKS (else absent; [] in a shot that did not ask): .labelling (the labelling
%         painted) .outlines (.labelling .coarser ([] none) .width .coarserWidth, pt) and, for raw pictures, .stops (the
%         colourbar's stops: an n x 2 cell {offset 0-1, colour}, [] with no field) .fit (.margin .shared .zoom (the
%         shared zoom, [] unless Shared)) or [] .raw (struct array, one per render in order: .label .view .hemisphere
%         .camera .frame (1-based) .time (s) .size_px [w h] (= StagePx) .mmPerPx (in the picture plane: orthographic,
%         exact at every depth) .fit (.zoom .aspect) or [] .colorbar (.lo .hi .name .unit .quantity .colormap .stops)
%         or [] .png .marks (the files in <Out>/<Name>/)). .svg, .pdf, .png are '' for raw pictures without a Width.
%   PIXELS (.scalebar.px, .paths.px, each raw render's marks) are in the raw render's pixels (<Out>/<Name>/<view>.png): x right, y down, the
%   origin at the render's top-left CORNER (pixel (r, c) of imread covers [c-1, c] x [r-1, r]).
%
% THE SCALE BAR is a curve on the surface, not a line on the image: the viewer's exact geodesic on the folded surface,
% cut at exactly 'ScaleBar' mm of its length, placed by Rheome's rule (an anchor among the front-most vertices near
% the lower right of the silhouette, ends 1-1.2 L away, flat to the view, horizontal and leftward first) where every
% point of it is seen; one bar a row. Limits, visibility (a ray cast toward the camera) and the bar are the viewer's
% (the Cortical Flow Editor's viewer/print/limits.ts, visible.ts, scalebar.ts), by the rules Rheome used.
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
    p.addParameter('ScaleBar', 20);
    p.addParameter('ScaleBarKind', '');
    p.addParameter('ScaleBarTile', []);
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
    p.addParameter('Selections', []);
    p.addParameter('BundlePaths', []);
    p.addParameter('Timeout', 60);
    p.parse(varargin{:});
    o = p.Results;
    if isempty(o.Out), error('rheome:snapshot:input', 'Out: the folder written is expected.'); end
    o.Out = char(o.Out);  o.Format = lower(string(o.Format));
    sb = o.ScaleBar;
    if (ischar(sb) || (isstring(sb) && isscalar(sb))) && strcmpi(sb, 'auto'), o.ScaleBar = 'auto';
    elseif isnumeric(sb) && isscalar(sb) && isreal(sb) && sb >= 0, o.ScaleBar = double(sb);
    else, error('rheome:snapshot:input', 'ScaleBar: mm >= 0, or "auto", is expected.');
    end
    if ~isempty(o.ScaleBarKind), o.ScaleBarKind = pickName(o.ScaleBarKind, {'geodesic', 'plane'}, 'ScaleBarKind'); else, o.ScaleBarKind = ''; end
    t = o.ScaleBarTile;
    if ~isempty(t) && ~(isnumeric(t) && isscalar(t) && t >= 1 && t == round(t))
        error('rheome:snapshot:input', 'ScaleBarTile: a render''s index, 1-based.');
    end

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
    editor = isfile(fullfile(vdir, 'index.html')) && ~isfile(fullfile(vdir, 'viewer.html'));   % THE EDITOR's build draws the snapshots
    if ~isfile(fullfile(vdir, 'viewer.html')) && ~editor, error('rheome:snapshot:noviewer', 'No viewer at %s (Viewer= or RHEOME_VIEWER): an Editor build (index.html) or a viewer build (viewer.html).', vdir); end
    drv = fullfile(vdir, 'driver', 'snapshot-driver.mjs');
    if ~isfile(drv)
        error('rheome:snapshot:driver', ['The viewer at %s has no snapshot driver (driver/snapshot-driver.mjs): rheome.snapshot needs ' ...
            'a Cortical Flow viewer build of 2026-11-10 or later.'], vdir);
    end
    cmaps = viewerColormaps(vdir);

    folded = o.Folded;
    if isempty(folded), folded = snames{1}; if any(strcmp(snames, 'pial')), folded = 'pial'; end, end
    folded = pickName(folded, snames, 'Folded');
    J = tok.theme.(tok.theme.default).page;
    shots = normalise(o.Shots, result, surfs, T, o.Fs, tok, cmaps, labellingNames(o.Labellings), o.ScaleBarTile);

    % one bundle, one job, one viewer run for every shot
    if isfield(surf, 'Vertices'), bs = surf; else, bs = surfs; end
    labs = o.Labellings;  st = structures(surfs.(snames{1}));      % the hemispheres, as the atlas the Editor reads them from
    if ~isempty(st) && ~any(strcmp(labellingNames(labs), 'Structures'))
        if isempty(labs), labs = {}; elseif isstruct(labs), labs = num2cell(labs(:)'); end
        labs = [labs(:)' {st}];
    end
    rheome.io.write.bundle(fullfile(work, 'bundle'), bs, result, 'Fs', o.Fs, 'Name', 'rheome.snapshot', ...
        'Labellings', labs, 'Points', o.Points, 'Paths', o.BundlePaths, 'Selections', o.Selections);
    list = cell(1, numel(shots));
    for k = 1:numel(shots), list{k} = struct('name', shots(k).Name, 'shot', viewerShot(shots(k), folded, o, J)); end
    job = struct();
    job.format = 'nxr-snapshot-job@1';  if editor, job.editor = vdir; else, job.viewer = vdir; end;  job.bundle = 'bundle';  job.out = 'out';
    job.chrome = chromePath(o.Chrome);  job.tokens = tokFile;  job.tokens_sha256 = tokSha;  job.tokens_path = tokShown;
    job.size = windowSize(double(tok.cortex.viewer.size_px(:)'), shots);  job.wait = o.Wait;  job.timeout = o.Timeout;  job.shots = list;
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
        for j = 1:numel(labels)
            movefile(fullfile(from, s.Name, [labels{j} '.png']), fullfile(rdir, [labels{j} '.png']), 'f');
            mk = fullfile(from, s.Name, [labels{j} '.marks.json']);      % raw pictures: each render's own marks
            if isfile(mk), movefile(mk, fullfile(rdir, [labels{j} '.marks.json']), 'f'); end
        end
        base = fullfile(o.Out, s.Name);
        composed = isfile(fullfile(from, [s.Name '.svg']));                % raw pictures without a Width: no panel
        exts = {'.marks.json', '.report.json'};  if composed, exts = [{'.svg'} exts]; end %#ok<AGROW>
        for ext = exts, movefile(fullfile(from, [s.Name ext{1}]), [base ext{1}], 'f'); end
        r = fromReport(jsondecode(fileread([base '.report.json'])), s, base, fullfile(rdir, strcat(labels(1:s.ntiles), '.png')), o.Format, composed, rdir);
        out = joinOut(out, r);
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

function shots = normalise(in, result, surfs, T, fs, tok, cmaps, labs, barTile)
% each shot checked, as the viewer will draw it (it is refused here, by name, rather than by the viewer); labs: the names
% of the Labellings given; barTile: the ScaleBarTile option (1-based) or []
    if isempty(in), in = {struct()}; end
    if isstruct(in), in = num2cell(in); end
    if ~iscell(in), error('rheome:snapshot:input', 'Shots: a struct array or a cell array of structs is expected.'); end
    known = {'Field', 'Unit', 'Label', 'Reference', 'Quantity', 'CLim', 'Colormap', 'Time', 'Frames', 'View', ...
        'Hemisphere', 'Surface', 'Vector', 'Settings', 'Name', 'Inset', 'Selection', 'Paths', 'Cortex', 'Labelling', 'Outlines', ...
        'StagePx', 'Fit', 'Dpi', 'MarksInRaster'};
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
        % the field painted, or none (a labelling, arrows or the bare cortex): none when Field is given empty, or when it is
        % not given and the shot paints a Labelling or the result has no real field
        s.Labelling = '';
        if isfield(a, 'Labelling') && ~isempty(a.Labelling), s.Labelling = pickLabelling(a.Labelling, labs, 'Labelling', s.Name); end
        if isfield(a, 'Field'), noField = isempty(a.Field); else, noField = ~isempty(s.Labelling) || isempty(real); end
        if noField
            for f = {'Unit', 'Label', 'Reference', 'Quantity', 'CLim', 'Colormap'}
                if isfield(a, f{1}) && ~isempty(a.(f{1})), error('rheome:snapshot:shot', '%s: %s is for a field''s colours; this shot paints no Field.', s.Name, f{1}); end
            end
            s.Field = '';  s.Unit = '';  s.Label = '';  s.Reference = '';
        else
            if isempty(real), error('rheome:snapshot:shot', 'The result has no real ([nV x T]) field to paint.'); end
            if ~isempty(s.Labelling), error('rheome:snapshot:shot', '%s: a Labelling is painted in place of a field''s colours: give a Field or a Labelling.', s.Name); end
            s.Field = pickName(g('Field', real{1}), real, 'Field');
            s.Unit = char(g('Unit', ''));
            if isempty(s.Unit), error('rheome:snapshot:unit', '%s: a colourbar needs a unit (Unit="a.u." only when there is none; say why in the legend).', s.Name); end
            s.Label = char(g('Label', s.Field));
            s.Reference = char(g('Reference', ''));
        end

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

        % views: the figure system's four, one of them, or a viewer camera with a hemisphere; or a list of these (a strip
        % of several views is a row a view)
        v = g('View', 'standard');
        if isnumeric(v), v = {v};
        elseif iscell(v) || (isstring(v) && numel(v) > 1), if ~iscell(v), v = cellstr(v); end, v = v(:)'; else, v = {char(v)}; end
        views4 = {'L lateral', 'left', 'left'; 'L medial', 'left', 'right'; 'R medial', 'right', 'left'; 'R lateral', 'right', 'right'};
        s.Hemisphere = '';  s.views = cell(0, 3);  s.view = {};
        for j = 1:numel(v)
            vj = v{j};
            if isnumeric(vj)                                  % a camera direction [x y z] (from the brain toward the camera)
                if ~(numel(vj) == 3 && isreal(vj) && all(isfinite(vj(:))) && any(vj(:) ~= 0))
                    error('rheome:snapshot:shot', '%s: View: a camera direction is [x y z], three finite numbers not all 0.', s.Name);
                end
                if isempty(s.Hemisphere), s.Hemisphere = pickName(g('Hemisphere', 'both'), {'both', 'left', 'right'}, 'Hemisphere'); end
                h = s.Hemisphere;  d = double(vj(:)');
                s.views = [s.views; {sprintf('%s, [%s]', h, strjoin(arrayfun(@jsonencode, d, 'UniformOutput', false), ', ')), h, d}];  s.view{end+1} = d; %#ok<AGROW>
                continue
            end
            vj = char(vj);
            if strcmpi(vj, 'standard'), s.views = [s.views; views4];  s.view{end+1} = 'standard'; %#ok<AGROW>
            elseif any(strcmpi(vj, views4(:, 1))), w = views4(strcmpi(vj, views4(:, 1)), :);  s.views = [s.views; w];  s.view{end+1} = w{1}; %#ok<AGROW>
            else
                cam = pickName(vj, {'left', 'right', 'top', 'bottom', 'front', 'back'}, 'View');
                if isempty(s.Hemisphere), s.Hemisphere = pickName(g('Hemisphere', 'both'), {'both', 'left', 'right'}, 'Hemisphere'); end
                h = s.Hemisphere;
                s.views = [s.views; {sprintf('%s, %s', h, cam), h, cam}];  s.view{end+1} = cam; %#ok<AGROW>
            end
        end
        if numel(unique(lower(s.views(:, 1)))) < size(s.views, 1), error('rheome:snapshot:shot', '%s: View: each view once (%s).', s.Name, strjoin(s.views(:, 1), ', ')); end
        s.Surface = pickName(g('Surface', firstOf(fieldnames(surfs))), fieldnames(surfs), 'Surface');

        % the colour scale checked: the quantity's colormap is the viewer's, a diverging CLim symmetric. The class when
        % none is given is the viewer's rule (both signs over the frames and hemispheres shown: signed), for the checks
        % only: the viewer is sent what was given.
        q = char(g('Quantity', ''));
        classes = {'signed', 'magnitude', 'power', 'phase', 'coupling', 'code'};
        s.Quantity = '';  s.Colormap = '';  s.CLim = [];
        if ~noField
            if isempty(q)
                if isempty(hemi), hemi = hemispheres(S1); end
                if isempty(hemi) && any(~strcmp(s.views(:, 2), 'both'))
                    error('rheome:snapshot:hemispheres', '%s: a view of one hemisphere needs the surface''s hemispheres (its ''Structures'' atlas: Cortex L, Cortex R).', s.Name);
                end
                if isempty(hemi), hemi = ones(size(S1.Vertices, 1), 1); end
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
        end
        s.Vector = char(g('Vector', ''));
        if ~isempty(s.Vector), s.Vector = pickName(s.Vector, vec, 'Vector'); end
        s.Settings = g('Settings', struct());
        if ~isstruct(s.Settings), error('rheome:snapshot:shot', '%s: Settings is a struct of viewer settings.', s.Name); end
        s.Selection = selections(g('Selection', []), s.Name, nV);
        s.Paths = paths(g('Paths', []), s.Name, nV);
        s.Inset = inset(g('Inset', []), s, real, nV);
        s.Cortex = [];
        if isfield(a, 'Cortex') && ~isempty(a.Cortex)
            c = a.Cortex;
            if islogical(c) && isscalar(c), s.Cortex = c;
            elseif isnumeric(c) && isscalar(c) && isreal(c) && c >= 0 && c <= 1, s.Cortex = double(c);
            else, error('rheome:snapshot:shot', '%s: Cortex: true, false or an opacity 0-1.', s.Name);
            end
        end
        s.Outlines = outlines(g('Outlines', []), s.Name, labs);
        [s.StagePx, s.Fit, s.Dpi, s.MarksInRaster] = raw(a, s);
        s.ntiles = size(s.views, 1) * numel(fr);
        if ~isempty(barTile) && barTile > s.ntiles
            error('rheome:snapshot:shot', '%s: ScaleBarTile %d: there are %d renders.', s.Name, barTile, s.ntiles);
        end
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

function o = outlines(a, name, labs)
% a shot's tile outlines, checked, as the viewer's Shot.outlines (codes are the labelling's own: no index to convert)
    o = [];
    if isempty(a), return; end
    id = 'rheome:snapshot:outlines';
    if ~isstruct(a) || ~isscalar(a), error(id, '%s: Outlines is a struct.', name); end
    bad = setdiff(fieldnames(a), {'Labelling', 'Coarser', 'Width', 'CoarserWidth', 'Select', 'SelectStyle'});
    if ~isempty(bad), error(id, '%s: %s is not an outlines field.', name, strjoin(bad, ', ')); end
    lab = get(a, 'Labelling', '');
    if isempty(lab), error(id, '%s: Outlines.Labelling: the labelling whose tiles are outlined.', name); end
    o = struct();
    o.labelling = pickLabelling(lab, labs, 'Outlines.Labelling', name);
    c = get(a, 'Coarser', []);
    if ~isempty(c)
        if isnumeric(c) && isscalar(c) && c >= 1 && c == round(c), o.coarser = double(c);
        elseif ischar(c) || (isstring(c) && isscalar(c)), o.coarser = pickLabelling(c, labs, 'Outlines.Coarser', name);
        else, error(id, '%s: Outlines.Coarser: a labelling''s name, or n >= 1 (the codes shifted right by n).', name);
        end
    end
    for f = {'Width', 'width'; 'CoarserWidth', 'coarserWidth'}'
        w = get(a, f{1}, []);
        if isempty(w), continue; end
        if ~(isnumeric(w) && isscalar(w) && w >= 0.05 && w <= 20), error(id, '%s: Outlines.%s: a line width in pt, 0.05-20.', name, f{1}); end
        o.(f{2}) = double(w);
    end
    sel = get(a, 'Select', []);
    if ~isempty(sel)
        if ~(isnumeric(sel) && all(sel(:) == round(sel(:)))), error(id, '%s: Outlines.Select: tile codes (integers) of Outlines.Labelling.', name); end
        o.select = num2cell(double(sel(:)'));
    end
    st = get(a, 'SelectStyle', '');
    if ~isempty(st), o.selectStyle = pickName(st, {'fill', 'extent'}, 'Outlines.SelectStyle'); end
end

function [px, fit, dpi, mir] = raw(a, s)
% RAW PICTURES: StagePx [w h] (each render exactly that), Fit (Margin, Shared), Dpi, MarksInRaster; as the viewer's Shot fields
    id = 'rheome:snapshot:raw';
    px = get(a, 'StagePx', []);  fit = [];  dpi = [];  mir = [];
    if ~isempty(px)
        if ~(isnumeric(px) && numel(px) == 2 && all(px(:) == round(px(:)) & px(:) >= 16 & px(:) <= 8192))
            error(id, '%s: StagePx: [width height], integer pixels 16-8192.', s.Name);
        end
        px = double(px(:)');
        if ~isempty(s.Inset), error(id, '%s: raw pictures (StagePx) take no inset (the host composes them).', s.Name); end
    end
    if isfield(a, 'Fit') && ~isempty(a.Fit) && ~(islogical(a.Fit) && isscalar(a.Fit) && ~a.Fit)
        if isempty(px), error(id, '%s: Fit is for raw pictures: give StagePx.', s.Name); end
        f = a.Fit;  fit = struct();
        if islogical(f) && isscalar(f), f = struct(); end
        if ~isstruct(f) || ~isscalar(f), error(id, '%s: Fit is a struct (Margin, Shared) or true.', s.Name); end
        bad = setdiff(fieldnames(f), {'Margin', 'Shared'});
        if ~isempty(bad), error(id, '%s: %s is not a fit field.', s.Name, strjoin(bad, ', ')); end
        m = get(f, 'Margin', []);
        if ~isempty(m)
            if ~(isnumeric(m) && isscalar(m) && m >= 0 && m < 0.5), error(id, '%s: Fit.Margin: 0 to 0.5 (a fraction of the half-size).', s.Name); end
            fit.margin = double(m);
        end
        if isfield(f, 'Shared') && ~isempty(f.Shared)
            if ~(isscalar(f.Shared) && (islogical(f.Shared) || (isnumeric(f.Shared) && any(f.Shared == [0 1])))), error(id, '%s: Fit.Shared: true or false.', s.Name); end
            fit.shared = logical(f.Shared);
        end
    end
    d = get(a, 'Dpi', []);
    if ~isempty(d)
        if isempty(px) || ~(isnumeric(d) && isscalar(d) && d >= 72 && d <= 1200), error(id, '%s: Dpi (72-1200) is for raw pictures (StagePx).', s.Name); end
        dpi = double(d);
    end
    m = get(a, 'MarksInRaster', []);
    if ~isempty(m)
        if isempty(px) || ~(isscalar(m) && (islogical(m) || (isnumeric(m) && any(m == [0 1]))))
            error(id, '%s: MarksInRaster (true or false) is for raw pictures (StagePx).', s.Name);
        end
        mir = logical(m);
    end
end

function x = pickLabelling(x, labs, what, name)
% a labelling's name among the Labellings given
    if isempty(labs)
        error('rheome:snapshot:option', '%s: %s: ''%s'': no labelling was given (Labellings=).', name, what, char(x));
    end
    x = pickName(x, labs, what);
end

function n = labellingNames(L)
% the names of the Labellings option (struct array or cell of structs)
    n = {};
    if isempty(L), return; end
    if isstruct(L), L = num2cell(L); end
    for k = 1:numel(L), if isstruct(L{k}) && isfield(L{k}, 'name'), n{end+1} = char(L{k}.name); end, end %#ok<AGROW>
end

function sz = windowSize(sz, shots)
% the window: the tokens' stage, or larger to hold the largest raw stage (StagePx) asked
    for k = 1:numel(shots), if ~isempty(shots(k).StagePx), sz = max(sz, shots(k).StagePx); end, end
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
% the viewer's Shot (the Cortical Flow Editor's viewer/print/panel.ts): 0-based frames and vertices, the width in
% mm; the keys in one order, so nxr.snapshot's job is this one byte for byte
    v = struct();
    if ~isempty(s.Field), v.field = s.Field;  v.unit = s.Unit;  v.label = s.Label; end
    if ~isempty(s.Reference), v.reference = s.Reference; end
    if ~isempty(s.Quantity), v.quantity = s.Quantity; end
    if ~isempty(s.CLim), v.clim = s.CLim; end
    if ~isempty(s.Colormap), v.colormap = s.Colormap; end
    v.frames = num2cell(s.frames - 1);
    v.views = s.view;
    if ~isempty(s.Hemisphere), v.hemisphere = s.Hemisphere; end
    v.surface = s.Surface;  v.folded = folded;
    if ~isempty(s.Vector), v.vector = s.Vector; end
    if ~isempty(s.Cortex), v.cortex = s.Cortex; end
    if ~isempty(s.Labelling), v.labelling = s.Labelling; end
    if ~isempty(s.Outlines), v.outlines = s.Outlines; end
    if ~isempty(fieldnames(s.Settings)), v.settings = s.Settings; end
    if ~isempty(s.Selection), v.selections = s.Selection; end
    if ~isempty(s.Paths), v.paths = s.Paths; end
    if ~isempty(s.Inset), v.inset = s.Inset; end
    v.scaleBarMm = o.ScaleBar;
    if ~isempty(o.ScaleBarKind), v.scaleBarKind = o.ScaleBarKind; end
    if ~isempty(o.ScaleBarTile), v.scaleBarTile = double(o.ScaleBarTile) - 1; end
    if ~isempty(s.StagePx), v.stagePx = s.StagePx; end
    if ~isempty(s.Fit), v.fit = s.Fit; end
    if ~isempty(s.Dpi), v.dpi = s.Dpi; end
    if ~isempty(s.MarksInRaster), v.marksInRaster = s.MarksInRaster; end
    if isempty(s.StagePx) || ~isempty(o.Width), v.widthMm = panelWidth(o.Width, s.ntiles + ~isempty(s.Inset), J); end   % raw pictures: a panel only with a Width
end

function r = fromReport(rep, s, base, renders, fmt, composed, rdir)
% out: the files, and the viewer's report with its indices made 1-based
    r = struct('name', s.Name, 'svg', '', 'pdf', '', 'png', '', 'renders', {renders(:)'});
    if composed, r.svg = [base '.svg']; end
    if composed && any(fmt == "pdf"), r.pdf = [base '.pdf']; end
    if composed && any(fmt == "png"), r.png = [base '.png']; end
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
    % only when asked: the labelling painted, the outlines drawn, and raw pictures' stops, fit and renders
    if isfield(rep, 'labelling'), r.labelling = rep.labelling; end
    if isfield(rep, 'outlines'), r.outlines = rep.outlines; end
    if isfield(rep, 'stops'), r.stops = stopsOf(rep.stops); end
    if isfield(rep, 'fit'), r.fit = rep.fit;  if ~isempty(r.fit), r.fit.shared = logical(r.fit.shared); end, end
    if isfield(rep, 'renders')
        R = struct('label', {}, 'view', {}, 'hemisphere', {}, 'camera', {}, 'frame', {}, 'time', {}, 'size_px', {}, 'mmPerPx', {}, ...
            'fit', {}, 'colorbar', {}, 'png', {}, 'marks', {});
        for j = 1:numel(rep.renders)
            q = rep.renders(j);  if iscell(rep.renders), q = rep.renders{j}; end
            cb = q.colorbar;  if ~isempty(cb), cb.stops = stopsOf(cb.stops); end
            mk = fullfile(rdir, [q.label '.marks.json']);  if ~isfile(mk), mk = ''; end
            R(end+1) = struct('label', q.label, 'view', q.view, 'hemisphere', q.hemisphere, 'camera', q.camera, 'frame', q.frame + 1, ...
                'time', q.time, 'size_px', q.size_px(:)', 'mmPerPx', q.mmPerPx, 'fit', q.fit, 'colorbar', cb, ...
                'png', fullfile(rdir, [q.label '.png']), 'marks', mk); %#ok<AGROW>
        end
        r.raw = R;
    end
end

function c = stopsOf(x)
% a colourbar's stops (JSON [[offset, colour], ...]) as an n x 2 cell {offset, colour}; [] for none
    c = [];
    if isempty(x), return; end
    if ~iscell(x), x = num2cell(x, 2); end
    c = cell(numel(x), 2);
    for i = 1:numel(x), q = x{i};  if ~iscell(q), q = num2cell(q); end, c(i, :) = {double(q{1}), char(q{2})}; end
end

function out = joinOut(out, r)
% shots asking for different extras (labelling, outlines, raw pictures) give one struct array: each lacks the others' as []
    if isempty(out), out = r; return; end
    for f = setdiff(fieldnames(r), fieldnames(out))', [out.(f{1})] = deal([]); end
    for f = setdiff(fieldnames(out), fieldnames(r))', r.(f{1}) = []; end
    out = [out orderfields(r, out(1))];
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
% 1 left, 2 right, 0 neither: the surface's ATLAS (its 'Structures' hemispheres, rheome.io.read.surface's Hemi / HemiLabel:
% Cortex L, Cortex R) -- a hemisphere is the set of vertices the atlas names, never found from the mesh (2026-11-10). [] without one
    h = [];
    if ~isfield(S, 'Hemi') || isempty(S.Hemi), return; end
    h = zeros(size(S.Vertices, 1), 1);
    for k = 1:numel(S.Hemi)
        lab = char(S.HemiLabel{k});
        if lab(end) == 'L', h(S.Hemi{k}) = 1; else, h(S.Hemi{k}) = 2; end
    end
end

function st = structures(S)
% the surface's hemispheres as the bundle's 'Structures' labelling (codes 1 Cortex L, 2 Cortex R, 0 neither) -- what the
% Editor shows a hemisphere from; [] without them
    st = [];
    h = hemispheres(S);
    if isempty(h), return; end
    st = struct('name', 'Structures', 'codes', int32(h), 'entries', struct('code', {0, 1, 2}, 'name', {'unassigned', 'Cortex L', 'Cortex R'}));
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
    drv = char(java.io.File(drv).getCanonicalPath());   % the real path: a viewer release reached through a symlink otherwise does nothing (driver before 0cfda31)
    [st, msg] = system(sprintf('"%s" "%s" "%s" 2>&1', node, drv, jf));
    if st ~= 0, error('rheome:snapshot:viewer', 'The viewer run failed: %s', strtrim(msg)); end
end
