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
%   'FieldQuantity' what the fields are (tokens quantity, default 'magnitude'): their colormap.
%   'FieldColormaps' {1 x nF} colormap per field panel, in place of the quantity's.
%                  ⚠ under 'hot' a BINARY INDICATOR is invisible -- zero is black and the panel
%                  reads as empty; pass a two-colour map for set-valued panels (an atlas scout).
%   'Bands'        {lo hi label} per row: shaded regions of the ruler (mm)
%   'View'         camera [az el] for every cortex panel (default [-90 10], left lateral)
%   'XLim'         ruler limits in mm (default from the marks, padded)
%   'Title'        the ruler's panel title
%   'Export'       file stem: write <stem>.pdf/.svg/.png at print size (nxr_save)
%   'Visible'      'on' (default) | 'off' for headless export
%   'Theme'        figure theme ('' = the tokens' default, Pergammon; 'neuroimage' for MS1)
% Sizes and colours are the figure tokens; double column, one panel per row; the fields' shared colourbar
% in their panel's key band.
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
    p.addParameter('FieldQuantity', 'magnitude');
    p.addParameter('Bands', {});
    p.addParameter('View', [-90 10]);
    p.addParameter('XLim', []);
    p.addParameter('Title', 'Spatial scales of the cortical problem');
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.addParameter('Theme', '');
    p.parse(varargin{:});
    o = p.Results;
    T = tokens(o.Theme);
    o.FieldQuantity = quantity(o.FieldQuantity, []);

    if ~istable(marks) || ~all(ismember({'name','mm','kind'}, marks.Properties.VariableNames))
        error('show:resolution:marks', 'marks must be a table with name, mm and kind.');
    end
    marks = sortrows(marks, 'mm');
    S = o.Surface;
    nP = numel(o.Patches);
    nF = 0;  if ~isempty(o.Fields), nF = size(o.Fields, 2); end
    hasCortex = ~isempty(S) && (nP > 0 || nF > 0);

    nCol = max([nP, nF, 1]);
    rows = {};  kinds = struct();  ttl = struct();  hr = [];  keyb = struct();
    fw = (T.page.width_mm.double - 2*T.grid.margin_mm) / nCol;
    viewH = T.grid.header_mm + 3 + 0.76*fw;          % header, view titles, the render aspect (1000 x 760)
    if hasCortex && nP > 0, rows{end+1} = 'patch'; hr(end+1) = viewH; end
    KB = T.legend.keyBand;
    if hasCortex && nF > 0, rows{end+1} = 'field'; hr(end+1) = viewH + KB.scaleRow_mm + KB.gap_below_mm; end
    rows{end+1} = 'ruler';  hr(end+1) = 45;
    keys = char('A' + (0:numel(rows)-1));
    for i = 1:numel(rows)
        kinds.(keys(i)) = ternary(strcmp(rows{i}, 'ruler'), 'plot', 'image');
        if strcmp(rows{i}, 'field'), keyb.(keys(i)) = 's'; end
        ttl.(keys(i)) = ternary(strcmp(rows{i}, 'ruler'), o.Title, ...
                        ternary(strcmp(rows{i}, 'patch'), 'Nested pyramid patches', 'Fields, each to its own peak'));
    end
    H = min(sum(hr) + (numel(hr)-1)*T.grid.gutter_mm + 2*T.grid.margin_mm, T.page.max_height_mm);
    F = nxr_figure('double', H, cellstr(keys'), 'Theme', o.Theme, 'Heights', hr, 'Kinds', kinds, 'Titles', ttl, 'Keys', keyb);
    hFig = F.fig;  hFig.Visible = o.Visible;

    for r = 1:numel(rows) - 1
        hp = F.ax.(keys(r));  P = hp.Position;  axis(hp, 'off');
        if strcmp(rows{r}, 'field')                  % one shared scale, in the panel's key band
            nicebar(nxr_colorbar(hp, 'field', 'ratio to its peak', o.FieldQuantity, [0 1], 'F', F), hp);
            hp.Position = P;
        end
        g = T.spacing.viewGap_mm / F.size_mm(1);  w = (P(3) - g*(nCol-1)) / nCol;
        hv = P(4) - 3/F.size_mm(2);
        for i = 1:nCol
            ax = axes(hFig, 'Position', [P(1) + (i-1)*(w+g), P(2), w, hv]);
            if strcmp(rows{r}, 'patch') && i <= nP
                d = zeros(size(S.Vertices,1), 1);  d(o.Patches{i}) = 1;
                rheome.show.surface(S, d, 'Parent', ax, 'View', o.View, 'CLim', [0 1], ...
                    'Colormap', [T.rgb.neutral_mist; T.rgb.atlas_space], 'Colorbar', false, 'Title', i_title(o.PatchTitles, i));
            elseif strcmp(rows{r}, 'field') && i <= nF
                X = double(o.Fields(:, i));  X = X / max(max(abs(X)), eps);   % each to its own peak: shapes compare
                cm = nxr_cmap(T.quantity.(o.FieldQuantity).colormap);
                if numel(o.FieldColormaps) >= i && ~isempty(o.FieldColormaps{i}), cm = o.FieldColormaps{i}; end
                rheome.show.surface(S, X, 'Parent', ax, 'View', o.View, 'CLim', [0 1], ...
                    'Colormap', cm, 'Colorbar', false, 'Title', i_title(o.FieldTitles, i));
            else
                axis(ax, 'off');
            end
        end
    end

    % ---- the ruler ----
    axR = F.ax.(keys(end));
    hold(axR, 'on');
    if isempty(o.XLim)
        xl = [0.85*min(marks.mm), 1.2*max(marks.mm)];
    else
        xl = o.XLim;
    end
    kinds = unique(marks.kind, 'stable');
    top = numel(kinds) + 2.7;
    for b = 1:size(o.Bands, 1)
        lo = o.Bands{b,1};  hi = min(o.Bands{b,2}, xl(2));
        patch(axR, [lo hi hi lo], [0 0 top top], ternary(mod(b,2) == 1, T.rgb.neutral_mist, T.rgb.neutral_snow), ...
              'EdgeColor', 'none', 'FaceAlpha', 0.6);
        text(axR, exp(mean(log([max(lo,xl(1)) hi]))), 0.16, o.Bands{b,3}, ...
             'HorizontalAlignment', 'center', 'Interpreter', 'none');
    end
    % up to 4 kinds told apart by shape AND fill (tokens: marker); more kinds repeat the cycle
    mk = T.marker.shapes.x_value;  mk = strrep(mk, 'D', 'd');
    fill = {T.rgb.neutral_ink, T.rgb.neutral_slate, T.rgb.neutral_paper, T.rgb.neutral_mist};
    for kk = 1:numel(kinds)
        sel = find(strcmp(string(marks.kind), string(kinds(kk))));
        j = mod(kk-1, numel(mk)) + 1;
        for q = 1:numel(sel)
            i = sel(q);  y = kk;
            plot(axR, marks.mm(i), y, mk{j}, 'MarkerSize', 2*T.marker.size_pt, ...
                 'MarkerFaceColor', fill{j}, 'MarkerEdgeColor', T.rgb.neutral_ink, 'LineWidth', T.marker.edge_pt);
            if strcmpi(char(string(marks.kind(i))), 'floor')
                plot(axR, [marks.mm(i) marks.mm(i)], [0.35 numel(kinds)+0.6], '--', ...
                     'Color', T.rgb.neutral_slate, 'LineWidth', T.line.leader_pt);  %#ok<*NBRAK2>
            end
            yt = y + 0.12 + 0.30*mod(i-1, 5);
            plot(axR, [marks.mm(i) marks.mm(i)], [y yt], '-', 'Color', T.rgb.neutral_slate, ...
                 'LineWidth', T.line.leader_pt);
            text(axR, marks.mm(i), yt, sprintf(' %s %.0f', string(marks.name(i)), marks.mm(i)), ...
                 'Interpreter', 'none', 'VerticalAlignment', 'bottom', 'Clipping', 'on');   % staggered rows, upright
        end
    end
    tk = [3 5 10 20 30 50 75 100 150 200 300 400 600 1000];
    tk = tk(tk >= xl(1) & tk <= xl(2));
    set(axR, 'XScale', 'log', 'XLim', xl, 'YLim', [0 top], ...
             'YTick', 1:numel(kinds), 'YTickLabel', cellstr(string(kinds)), ...
             'XTick', tk, 'XTickLabel', compose('%g', tk), 'XMinorTick', 'off');
    nxr_axes_style(axR, 'length on the cortex (mm, log scale)');
    P = axR.Position;  inset = 8 / F.size_mm(1);       % room for the kind names on the y axis
    axR.Position = [P(1) + inset, P(2), P(3) - inset, P(4)];
    hold(axR, 'off');
    finish(F, o.Title, o.Export);
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

function t = i_title(C, i)
    if numel(C) >= i, t = C{i}; else, t = ''; end
end

% Author: Diellor Basha, 2026
