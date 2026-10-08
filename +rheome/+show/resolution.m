function hFig = resolution(marks, varargin)
% SHOW.RESOLUTION  The spatial scales of the problem on one log ruler, with the cortex above it.
%
%   rheome.show.resolution(marks)
%   rheome.show.resolution(marks, 'Surface', S, 'Patches', P, 'PatchTitles', t)
%   rheome.show.resolution(marks, ..., 'Fields', X, 'FieldTitles', t2, 'Bands', bands)
%
% Renders what `rheome.inverse.resolution` measures next to what the analysis wants to measure, so the
% gap is visible rather than tabulated: a row of nested pyramid patches, a row of blobs drawn on
% the SAME mesh with the SAME camera (a real point-spread function beside a real atlas scout),
% and underneath them a logarithmic ruler carrying every length with the floors marked.
%
% ⭐ THE POINT OF THE SHARED CAMERA. A point-spread function and a Desikan-Killiany scout are
% usually shown in different figures at different zooms, which is exactly how a 52 mm PSF and a
% 53 mm scout come to look like different orders of magnitude. Drawn on one mesh at one camera
% they are the same size, which is the finding.
%
% ⚠ THIS FUNCTION COMPUTES NOTHING. It takes lengths and fields that were measured elsewhere
% (rheome.inverse.resolution, rheome.geom.tree, rheome.load.atlas) exactly as +show always has. If a number here
% disagrees with a docstring, the docstring is the record.
%
% INPUTS:
%   marks  table with .name (string), .mm (double), .kind (string). `kind` selects the marker
%          and the row the label sits on; the conventional set is "instrument", "pyramid",
%          "atlas", "floor". Order does not matter, it is sorted by .mm.
%
% OPTIONS (name/value):
%   'Surface'      surface struct (.Vertices, .Faces) for the cortex rows; [] draws the ruler alone
%   'Patches'      {1 x nP} vertex index vectors, drawn as filled patches in row 1
%   'PatchTitles'  {1 x nP} titles
%   'Fields'       [nV x nF] per-vertex scalars, drawn in row 2 with a SHARED colour scale
%                  after each column is normalised to its own peak
%   'FieldTitles'  {1 x nF} titles
%   'FieldColormaps' {1 x nF} colormap per field panel (default hot for all). ⚠ A BINARY
%                  INDICATOR UNDER 'hot' IS INVISIBLE -- zero is black and the panel reads as
%                  empty. Pass a sequential map for set-valued panels (an atlas scout) and keep
%                  hot for graded ones (a point-spread function).
%   'Bands'        {lo hi label} per row: shaded regions of the ruler (mm)
%   'View'         camera [az el] for every cortex panel (default [-90 10], left lateral)
%   'XLim'         ruler limits in mm (default from the marks, padded)
%   'Title'        figure title
%   'Visible'      'on' (default) | 'off' for headless export
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.inverse.resolution, rheome.geom.tree, rheome.show.surface, rheome.show.filmstrip
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Surface', []);
    p.addParameter('Patches', {});
    p.addParameter('PatchTitles', {});
    p.addParameter('Fields', []);
    p.addParameter('FieldTitles', {});
    p.addParameter('FieldColormaps', {});
    p.addParameter('Bands', {});
    p.addParameter('View', [-90 10]);
    p.addParameter('XLim', []);
    p.addParameter('Title', 'Spatial scales of the cortical problem');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;

    if ~istable(marks) || ~all(ismember({'name','mm','kind'}, marks.Properties.VariableNames))
        error('show:resolution:marks', 'marks must be a table with name, mm and kind.');
    end
    marks = sortrows(marks, 'mm');
    S = o.Surface;
    nP = numel(o.Patches);
    nF = 0;  if ~isempty(o.Fields), nF = size(o.Fields, 2); end
    hasCortex = ~isempty(S) && (nP > 0 || nF > 0);

    nCol = max([nP, nF, 1]);
    nRow = double(nP > 0) + double(nF > 0) + 1;            % + the ruler
    hFig = figure('Visible', o.Visible, 'Color', 'w', ...
                  'Position', [60 60 max(300*nCol, 1200) 260*(nRow-1) + 360]);
    tl = tiledlayout(hFig, nRow, nCol, 'TileSpacing', 'none', 'Padding', 'compact');
    title(tl, o.Title, 'FontWeight', 'bold');

    r = 0;
    if nP > 0
        r = r + 1;
        for i = 1:nP
            ax = nexttile(tl, (r-1)*nCol + i);
            d = zeros(size(S.Vertices,1), 1);  d(o.Patches{i}) = 1;
            rheome.show.surface(S, d, 'Parent', ax, 'View', o.View, 'CLim', [0 1], ...
                'Colormap', i_seq(), 'Colorbar', false, 'Title', i_title(o.PatchTitles, i));
        end
        for i = nP+1:nCol, ax = nexttile(tl, (r-1)*nCol + i); axis(ax, 'off'); end
    end
    if nF > 0
        r = r + 1;
        X = double(o.Fields);
        X = X ./ max(max(abs(X), [], 1), eps);             % each to its own peak: shapes compare
        for i = 1:nF
            ax = nexttile(tl, (r-1)*nCol + i);
            cm = hot(256);
            if numel(o.FieldColormaps) >= i && ~isempty(o.FieldColormaps{i})
                cm = o.FieldColormaps{i};
            end
            rheome.show.surface(S, X(:,i), 'Parent', ax, 'View', o.View, 'CLim', [0 1], ...
                'Colormap', cm, 'Colorbar', false, 'Title', i_title(o.FieldTitles, i));
        end
        for i = nF+1:nCol, ax = nexttile(tl, (r-1)*nCol + i); axis(ax, 'off'); end
    end

    % ---- the ruler ----
    axR = nexttile(tl, r*nCol + 1, [1 nCol]);
    hold(axR, 'on');
    if isempty(o.XLim)
        xl = [0.85*min(marks.mm), 1.2*max(marks.mm)];
    else
        xl = o.XLim;
    end
    kinds = unique(marks.kind, 'stable');
    yOf = containers.Map(cellstr(string(kinds)), num2cell(1:numel(kinds)));
    for b = 1:size(o.Bands, 1)
        lo = o.Bands{b,1};  hi = min(o.Bands{b,2}, xl(2));
        patch(axR, [lo hi hi lo], [0 0 numel(kinds)+2.7 numel(kinds)+2.7], i_bandcolor(b), ...
              'EdgeColor', 'none', 'FaceAlpha', 0.16);
        text(axR, exp(mean(log([max(lo,xl(1)) hi]))), 0.16, o.Bands{b,3}, ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'FontAngle', 'italic', ...
             'Color', [0.2 0.2 0.2], 'Interpreter', 'none');
    end
    mk = {'o','s','^','d','v','p'};
    % ⚠ LABELS COLLIDE ON A LOG AXIS. Several floors land within a few mm of each other
    % (99, 103, 110 on this cortex), so the text is staggered in height by rank within each
    % kind rather than drawn at one offset, and the leader line says which mark it belongs to.
    % marks are already sorted by mm, so the row index IS the rank in x: staggering on it
    % guarantees that two marks close in x are never drawn at the same label height.
    for kk = 1:numel(kinds)
        sel = find(strcmp(string(marks.kind), string(kinds(kk))));
        for q = 1:numel(sel)
            i = sel(q);  y = kk;
            plot(axR, marks.mm(i), y, mk{min(kk, numel(mk))}, 'MarkerSize', 8, ...
                 'MarkerFaceColor', i_kindcolor(kk), 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
            if strcmpi(char(string(marks.kind(i))), 'floor')
                plot(axR, [marks.mm(i) marks.mm(i)], [0.35 numel(kinds)+0.6], '--', ...
                     'Color', [0.15 0.15 0.15 0.5], 'LineWidth', 1.1);  %#ok<*NBRAK2>
            end
            yt = y + 0.12 + 0.30*mod(i-1, 5);
            plot(axR, [marks.mm(i) marks.mm(i)], [y yt], '-', 'Color', [0.5 0.5 0.5], ...
                 'LineWidth', 0.5);
            text(axR, marks.mm(i), yt, sprintf(' %s %.0f', string(marks.name(i)), marks.mm(i)), ...
                 'Rotation', 30, 'FontSize', 8.5, 'Interpreter', 'none', ...
                 'VerticalAlignment', 'bottom', 'Clipping', 'on');
        end
    end
    tk = [3 5 10 20 30 50 75 100 150 200 300 400 600 1000];
    tk = tk(tk >= xl(1) & tk <= xl(2));
    set(axR, 'XScale', 'log', 'XLim', xl, 'YLim', [0 numel(kinds)+2.7], ...
             'YTick', 1:numel(kinds), 'YTickLabel', cellstr(string(kinds)), 'Box', 'on', ...
             'XTick', tk, 'XTickLabel', compose('%g', tk), 'XGrid', 'on', 'GridAlpha', 0.3, ...
             'XMinorGrid', 'off', 'FontSize', 9);
    xlabel(axR, 'length on the cortex (mm, log scale) -- diameter for a patch, wavelength for a pattern');
    hold(axR, 'off');
end

function t = i_title(C, i)
    if numel(C) >= i, t = C{i}; else, t = ''; end
end

function c = i_seq()
    c = [linspace(0.82,0.10,256)', linspace(0.85,0.30,256)', linspace(0.90,0.65,256)'];
end

function c = i_kindcolor(i)
    P = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56; 0.93 0.69 0.13];
    c = P(mod(i-1, size(P,1)) + 1, :);
end

function c = i_bandcolor(i)
    P = [0.80 0.20 0.20; 0.95 0.75 0.20; 0.20 0.65 0.30; 0.45 0.45 0.45];
    c = P(mod(i-1, size(P,1)) + 1, :);
end
% Author: Diellor Basha, 2026
