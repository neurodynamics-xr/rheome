%% How to make figures of the cortex
% *How-to guide.* A cortical pattern is read from its figures, and figures
% compare only when they are drawn the same way. Rheome draws every figure
% of the cortex with one renderer, the Cortical Flow viewer: |rheome.view|
% to look at a result, |rheome.snapshot| to print it by the lab's figure
% tokens. This guide puts a flow map, the current, a tile level and the
% critical points on one hemisphere and hands them to the viewer. It uses
% the demo study (see <howto_ingest_brainstorm.html How to import a
% Brainstorm study>). MATLAB graphics (|rheome.show|) are kept for quick
% looks, tests and machines without WebGPU; no figure is made with them.

name = 'cfdemo';
if ~rheome.load.has(name)                     % the demo study, imported as in the ingest guide
    demo = cfdemo_study();
    rheome.import.surface(name, demo.cortexFile);
    rheome.import.study(name, demo.studyDir, demo.dataName);
    rheome.import.bases(name, 100, 100, struct('overwrite', true));
end
ctx = rheome.flow.context(name, [], 100);
K   = rheome.flow.build(ctx);
[~, t0] = max(vecnorm(ctx.F));
b = ctx.F(:, t0);

%% The result: a map and the current on one hemisphere
% Each field is a column of values a vertex (real) or three a vertex
% (vector, interleaved x, y, z). The viewer paints a real field and any
% form of a vector one (|v|, its divergence, its rotation), and draws the
% vector as arrows.

B = rheome.load.bases(name, [], [], 100);
H = B.L;
Psi = K.stream.vertexOperator * b;
rows = reshape((H.gv(:)' - 1)*3 + (1:3)', [], 1);
J = ctx.currentKernel(rows, :) * b;                % [3nVh x 1]
result = struct('psi', Psi(H.gv), 'J', J);

%% Tiles by level and critical points
% |rheome.viz.tiles| turns a level of the cortical pyramid
% (|rheome.geom.tree|) into a labelling the viewer paints;
% |rheome.viz.points| turns critical points into point sets, one a type, so
% each type has its key.

T  = rheome.geom.tree(H.S, MaxDepth=4);
L  = rheome.viz.tiles(T, 4);
cp = rheome.detect.criticalPoints(J, H.S);
P  = rheome.viz.points(cp, H.S);

%% Look at it
% |rheome.view| writes the bundle and opens it in the viewer (inside MATLAB
% R2026a+, or the Cortical Flow app). Layers name a field, a form and how
% it is drawn. Without a viewer it falls back to a quick look. Here
% |Viewer="none"| writes the bundle only, so the guide runs anywhere.

out = rheome.view(H.S, result, Viewer="none", Labellings=L, Points=P, ...
    Layers={{"psi", "f", "colour"}, {"J", "xt", "arrows"}});
assert(isfile(out.manifest))

%%
% With the viewer open (|Viewer="auto"|), |out.screenshot(file)| writes the
% view as a PNG and, beside it, its marks sidecar |<file>.marks.json|: every
% overlay drawn with its key, which the figure check reads.

%% Print it
% |rheome.snapshot| prints panels by the figure tokens: the four standard
% views, the colormap and limits by quantity, the colourbar with its unit,
% the geodesic scale bar, strips over time, multiscale insets (see
% <howto_snapshot.html How to make print-ready cortical panels>). It needs
% the viewer, Chrome and Node.js.

if ~isempty(getenv('RHEOME_VIEWER')) && isfile(fullfile(getenv('RHEOME_VIEWER'), 'viewer.html'))
    shot = struct('Field', "psi", 'Label', "stream function", 'Unit', "a.u.", 'View', "left", ...
        'Hemisphere', "left", 'Settings', struct('labelling', L.name));
    rheome.snapshot(H.S, result, Shots=shot, Out=rheome.load.outpath('figures', 'stream', name), ...
        Labellings=L, Points=P);
end

%% A quick look without the viewer
% |rheome.show.surface| and |rheome.show.vectors| draw in MATLAB graphics.
% Use them to check a result, in tests, or where there is no WebGPU; a
% figure for a paper or a talk comes from |rheome.snapshot|.

ax = rheome.show.surface(H.S, Psi(H.gv), 'Title', 'Stream function (quick look)', 'Visible', 'off');
close(ancestor(ax, 'figure'))

%% See also
% <howto_view.html How to view a result in 3D>,
% <howto_snapshot.html How to make print-ready cortical panels>,
% <reference/rheome.viz.tiles.html rheome.viz.tiles>,
% <reference/rheome.viz.points.html rheome.viz.points>,
% <reference/rheome.viz.marks.html rheome.viz.marks>,
% <reference/rheome.show.surface.html rheome.show.surface>.
%
% _Written for Rheome @COMMIT@._
