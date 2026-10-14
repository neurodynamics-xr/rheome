function [F, ax, lim] = striprow(parent, grid, row, n, visible, label, bar)
% STRIPROW  One row of a snapshot strip on the house grid: a panel (letter A, B ... per row) holding
% n cortex views side by side, one view gap apart, over ONE horizontal colourbar centred under the row.
%   parent  [] (a new journal-width figure with grid(1) rows) or a figure made by a strip before
%   bar     {quantity, unit, kind, lim, colormap} for nxr_colorbar (kind 'signed' | 'magnitude' ...)
%   F       the house figure (nxr_figure);  ax [1 x n] the view axes, left to right
%   lim     the colour limits used (bar{4} widened to round values): give them to the views
%
% Author: Diellor Basha, 2026
    T = tokens();
    under = 10;                                      % mm under the views: bar, its ticks and its label
    over = 3;                                        % mm over the views: their t = ... labels
    if isempty(parent)
        fw = (T.journalSel.width_mm.double - 2*T.grid.margin_mm) / grid(2);
        rowH = T.grid.header_mm + over + 0.76*fw + under;   % the viewer's render aspect (1000 x 760)
        h = min(grid(1)*rowH + (grid(1)-1)*T.grid.gutter_mm + 2*T.grid.margin_mm, T.journalSel.max_height_mm);
        keys = cellstr(char('A' + (0:grid(1)-1))');
        F = nxr_figure('double', h, keys, 'Kinds', cell2struct(repmat({'image'}, grid(1), 1), keys, 1));
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
            'VerticalAlignment', 'top', 'Margin', 0);
    end
    colormap(hp, bar{5});
    lim = bar{4};  m = max(abs(lim));
    if strcmp(bar{3}, 'signed') && m > 0          % -L, 0, +L with L at two significant digits
        e = 10^(floor(log10(m)) - 1);  lim = [-1 1] * ceil(m/e) * e;
    end
    cb = nxr_colorbar(hp, bar{1}, bar{2}, bar{3}, lim, '', 'keep');
    if ~strcmp(bar{3}, 'signed'), nicebar(cb, hp);  lim = hp.CLim; end
    L = T.spacing.colorbarLength_frac * P(3);  th = T.spacing.colorbarThickness_mm / H;
    cb.Position = [P(1) + (P(3) - L)/2, P(2) + 6/H, L, th];
    hp.Position = P;                                 % the colourbar must not resize the holder
    g = T.spacing.viewGap_mm / W;  w = (P(3) - g*(n-1)) / n;  y = P(2) + under/H;  hv = P(4) - (under + over)/H;
    ax = gobjects(1, n);
    for i = 1:n, ax(i) = axes(F.fig, 'Position', [P(1) + (i-1)*(w+g), y, w, hv]); end
end
