%% How to make print-ready cortical panels
% *How-to guide.* A cortical pattern in a paper is read from still panels,
% and a reader compares panels only when they are drawn the same way: the
% same views, the same colour scale for the same quantity, a colourbar that
% says what the colour measures and against what, and a scale bar that
% measures distance on the cortex. |rheome.snapshot| draws every cortical
% panel that way. The Cortical Flow viewer (the renderer of |rheome.view|)
% draws the views in headless Chrome. The layout follows the lab's figure
% tokens (|library/styles/tokens/figure-tokens.json|): colormap and limits
% by quantity, colourbar, view labels, a geodesic scale bar, and time
% strips with their interval stated once.
%
% It needs the viewer (|RHEOME_VIEWER|), Chrome and Node.js 22+. Nothing
% else in Rheome does.

%% A result
% A cortex from the data cache, a signed field (anterior position) and a
% bump that moves over 100 ms at 100 Hz. Any [nV x T] field works the same way.

S = rheome.load.surface('omega_sub0002');
V = S.Vertices;  c = mean(V);
x = (V(:, 1) - c(1)) * 1e3;                                   % mm, signed
T = 10;  m = zeros(size(V, 1), T);
for t = 1:T
    m(:, t) = exp(-sum((V - (c + [0.06 - 0.012*t, 0.05, 0.03])).^2, 2) / (2*0.015^2));
end
result = struct('position', repmat(x, 1, T), 'bump', m);

%% The shot list
% One struct per panel. |Unit| is required, because a colourbar without a
% unit is not allowed. |Quantity| picks the colormap and the limits from the
% tokens. A signed field gets PuOr_r with limits symmetric around 0;
% magnitudes get cividis from 0. |View| defaults to the four standard
% views: L lateral, L medial, R medial, R lateral. A vector of times makes a
% strip of one view.

shots = { ...
    struct('Field', "position", 'Label', "anterior position", 'Unit', "mm", 'Name', "four-views"), ...
    struct('Field', "bump", 'Label', "amplitude", 'Unit', "a.u.", 'Time', (1:2:9)/100, ...
           'View', "L lateral", 'Name', "strip")};

%% Draw them
vdir = getenv('RHEOME_VIEWER');
if isempty(vdir) || ~isfile(fullfile(vdir, 'viewer.html'))
    disp('No viewer here (set RHEOME_VIEWER to the dist-viewer folder); nothing drawn.')
    return
end
out = rheome.snapshot(S, result, Shots=shots, Out=tempname, Fs=100);

%% What comes out
% Each shot is written as an SVG (editable text and bars, the views at
% their rendered resolution), a PDF and a PNG at the journal's width and
% dpi. The viewer's raw renders are kept too. |legend| is the sentence the
% figure legend needs: views, surface, the scale bar's distance and the
% interval.

for k = 1:numel(out)
    fprintf('%s: %s\n  %s\n  printed at %.0f dpi\n', out(k).name, out(k).pdf, out(k).legend, out(k).dpi);
end
imshow(imread(out(1).png))

%% The scale bar
% The bar is a curve on the cortex, not a straight line on the image. It is
% the shortest path along the folded surface's mesh edges, cut at exactly
% 20 mm of its length, and drawn where it lies in the view. A mesh-edge path
% is 3-13 % longer than the geodesic between its ends, so the ends are a
% little closer than 20 mm, never farther. On an inflated surface the legend
% says the distance is measured on the folded one.

disp(out(1).scalebar)
