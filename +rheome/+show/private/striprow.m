function [F, ax, lim] = striprow(parent, grid, row, n, visible, label, bar, theme)
% STRIPROW  One row of a snapshot strip on the house grid: a panel (letter A, B ... per row) holding
% n cortex views side by side, one view gap apart, with ONE colourbar in the panel's key band.
%   parent  [] (a new column-width figure with grid(1) rows) or a figure made by a strip before
%   bar     {label, unit, quantity, lim, colormap} for nxr_colorbar: quantity is a token key (divergence,
%           magnitude, amplitude ...: its own colormap); colormap [] = the quantity's, else kept as given
%   theme   figure theme ('' = the tokens' default; 'neuroimage' for MS1)
%   F       the house figure (nxr_figure);  ax [1 x n] the view axes, left to right
%   lim     the colour limits used (bar{4} widened to round values): give them to the views
%
% Author: Diellor Basha, 2026
    if nargin < 8, theme = ''; end
    T = tokens(theme);
    over = 3;                                        % mm over the views: their t = ... labels
    if isempty(parent)
        KB = T.legend.keyBand;
        fw = (T.page.width_mm.double - 2*T.grid.margin_mm) / grid(2);
        rowH = T.grid.header_mm + KB.scaleRow_mm + KB.gap_below_mm + over + 0.76*fw;   % viewer aspect 1000 x 760
        h = min(grid(1)*rowH + (grid(1)-1)*T.grid.gutter_mm + 2*T.grid.margin_mm, T.page.max_height_mm);
        keys = cellstr(char('A' + (0:grid(1)-1))');
        F = nxr_figure('double', h, keys, 'Theme', theme, 'Kinds', cell2struct(repmat({'image'}, grid(1), 1), keys, 1), ...
            'Keys', cell2struct(repmat({'s'}, grid(1), 1), keys, 1));
        F.fig.Visible = visible;
    else
        F = getappdata(parent, 'nxr');
        if isempty(F), error('rheome:show:parent', '''Parent'' must be a figure made by a rheome.show strip.'); end
    end
    k = char('A' + row - 1);
    hp = F.ax.(k);  P = hp.Position;  W = F.size_mm(1);  H = F.size_mm(2);
    if ~isempty(label)
        r = F.rect_mm.(k);
        annotation(F.fig, 'textbox', [(r(1) + 4.5)/W, 1 - (r(2) + T.grid.header_mm)/H, (r(3) - 4.5)/W, ...
            (T.grid.header_mm - 0.55)/H], 'String', label, 'EdgeColor', 'none', 'Interpreter', 'none', ...
            'FontSize', T.type.panelTitle.size_pt, 'FontWeight', 'bold', 'Color', T.rgb.neutral_ink, ...
            'VerticalAlignment', 'top', 'Margin', 0, 'Tag', ['nxr-title-' k]);
    end
    q = bar{3};  keep = ~isempty(bar{5});
    if keep, colormap(hp, bar{5}); end
    lim = bar{4};  m = max(abs(lim));
    diverging = strcmp(F.tokens.colour.colormap.(matlab.lang.makeValidName(F.tokens.quantity.(q).colormap)).scale, 'diverging');
    if diverging && m > 0                            % -L, 0, +L with L at two significant digits
        e = 10^(floor(log10(m)) - 1);  lim = [-1 1] * ceil(m/e) * e;
    end
    cb = nxr_colorbar(hp, bar{1}, bar{2}, q, lim, 'F', F, 'Keep', keep);
    if ~diverging, nicebar(cb, hp);  lim = hp.CLim; end
    if ~keep, bar{5} = colormap(hp); end
    hp.Position = P;                                 % the colourbar must not resize the holder
    g = T.spacing.viewGap_mm / W;  w = (P(3) - g*(n-1)) / n;  hv = P(4) - over/H;
    ax = gobjects(1, n);
    for i = 1:n, ax(i) = axes(F.fig, 'Position', [P(1) + (i-1)*(w+g), P(2), w, hv]); colormap(ax(i), bar{5}); end
end
