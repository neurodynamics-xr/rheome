%% How to make figures
% *How-to guide.* This guide draws flow maps on one hemisphere, draws the
% current as arrows, and saves the figures where the toolbox keeps its
% outputs. It uses the demo study (see <howto_ingest_brainstorm.html How to
% import a Brainstorm study>).

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

%% Draw a map on one hemisphere
% Take the hemisphere's surface and its vertices' values. |rheome.show.surface|
% chooses a diverging colormap, symmetric about zero, for signed maps.

B = rheome.load.bases(name, [], [], 100);
H = B.L;
Psi = K.stream.vertexOperator * b;
ax = rheome.show.surface(H.S, Psi(H.gv), 'Title', 'Stream function, left hemisphere');

%% Draw the current as arrows
% |rheome.show.vectors| colours the surface by the current's magnitude and draws
% arrows at the strongest vertices.

rows = reshape((H.gv(:)' - 1)*3 + (1:3)', [], 1);
J = ctx.currentKernel(rows, :) * b;                % [3nVh x 1], interleaved x, y, z
rheome.show.vectors(H.S, J, 'NumArrows', 200, 'Title', 'Current');

%% Put several maps in one figure
% Pass the axes with |'Parent'|:

fig = figure;
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
rheome.show.surface(H.S, K.potential.vertexOperator(H.gv, :) * b, 'Parent', nexttile(tl), 'Title', 'Scalar potential');
rheome.show.surface(H.S, Psi(H.gv), 'Parent', nexttile(tl), 'Title', 'Stream function');

%% Save it
% |rheome.load.outpath| returns the file's place under |rheome.load.outroot()|, creating
% the folder; set |RHEOME_OUT| to move the whole output tree.

file = rheome.load.outpath('figures', 'potentials.png', name);
exportgraphics(fig, file, 'Resolution', 150);
assert(isfile(file))

%%
% For scripts without a display, pass |'Visible', 'off'| to |rheome.show.surface|
% and |rheome.show.vectors| and export as above.
%
%% See also
% <howto_flow_maps.html How to compute flow maps>,
% <reference/rheome.show.surface.html rheome.show.surface>,
% <reference/rheome.show.vectors.html rheome.show.vectors>,
% <reference/rheome.load.outpath.html rheome.load.outpath>.
%
% _Written for Rheome @COMMIT@._
