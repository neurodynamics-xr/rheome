function [hFig, G] = levelstack(duration, levels, varargin)
% SHOW.LEVELSTACK  The dyadic levels as what they are: the grid over time, one row of tiles per level.
%
%   hFig = rheome.show.levelstack(duration, levels)
%   hFig = rheome.show.levelstack(duration, levels, 'Base', 1/64, 'Slice', 3, 'Tile', [3 5])
%   hFig = rheome.show.levelstack(duration, levels, 'Field', vals, 'Quantity', 'power', 'Label', 'power', 'Unit', 'dB')
%
% A level is not a measurement: it is the spacing of the grid on which measurements are made, as the grid spacing is
% in Euclidean space (the step of the differential operators, the support of a tile). So levels get no colormap. Level
% j tiles time into intervals of Base * 2^j s; each tile of level j is two tiles of level j - 1, and the rows stacked
% form the constant-Q volume (time x level; over the cortex each level also has its spatial tile, 'SpaceMM'). One
% level is a SLICE through that stack ('Slice', drawn in ink and marked on the z-stack indicator, rheome.show.zstack),
% one tile a cell of it ('Tile'). A FIELD measured on the tiles (a scalogram: a value per tile per level) is painted in
% the tiles with its quantity's colormap and a colourbar: the field is coloured, the grid is not.
%
% INPUTS:
%   duration  s, the time shown (from 0)
%   levels    vector of integer levels (e.g. 0:6), drawn bottom (finest) to top (coarsest)
%
% OPTIONS (name/value):
%   'Base'      s, the tile of level 0 (default: duration / 2^max(levels), one tile at the coarsest level)
%   'Slice'     a level drawn as the current slice (default none)
%   'Tile'      [level index] (index 1-based along time) a tile outlined (default none)
%   'SpaceMM'   [1 x nL] the spatial tile size at each level, mm (the cortical pyramid's), on the right axis
%   'Field'     {1 x nL} a value per tile of each level (a scalogram on the grid); with Quantity, Label and Unit
%   'Quantity'  the field's quantity (figure tokens: power, magnitude, divergence ...; default 'power')
%   'Label','Unit'  the field's name and unit on the colourbar (required with Field)
%   'Width'     'single' | 'onehalf' | 'double' (default 'double')
%   'Export','Visible','Title'  as every rheome.show figure
%
% OUTPUT:
%   hFig the figure;  G struct per level: .level .t0 .t1 (s, the tiles' edges)
%
% Example:
%   rheome.show.levelstack(4, 0:6, 'Slice', 3, 'Tile', [3 5], 'SpaceMM', 122 ./ 2.^(6:-1:0));
%
% See also: rheome.show.zstack, rheome.show.timetracks
%
% Author: Diellor Basha, 2026
    T = tokens();
    p = inputParser;
    p.addParameter('Base', []);
    p.addParameter('Slice', []);
    p.addParameter('Tile', []);
    p.addParameter('SpaceMM', []);
    p.addParameter('Field', {});
    p.addParameter('Quantity', 'power');
    p.addParameter('Label', '');
    p.addParameter('Unit', '');
    p.addParameter('Width', 'double');
    p.addParameter('Title', 'The dyadic level stack');
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;
    levels = double(levels(:)');
    if ~(isscalar(duration) && duration > 0), error('rheome:show:levelstack', 'duration: s, > 0.'); end
    if isempty(levels) || any(levels ~= round(levels)) || any(diff(levels) <= 0), error('rheome:show:levelstack', 'levels: increasing integers.'); end
    base = o.Base;  if isempty(base), base = duration / 2^max(levels); end
    nL = numel(levels);
    hasField = ~isempty(o.Field);
    if hasField && (numel(o.Field) ~= nL || isempty(o.Label) || isempty(o.Unit))
        error('rheome:show:levelstack', 'Field: one vector per level, with Label and Unit.');
    end
    if ~isempty(o.SpaceMM) && numel(o.SpaceMM) ~= nL, error('rheome:show:levelstack', 'SpaceMM: one size per level.'); end
    if ~isempty(o.Slice) && ~ismember(o.Slice, levels), error('rheome:show:levelstack', 'Slice: one of the levels.'); end

    keys = 'g';  if hasField, keys = 'gs'; end   % the colour row last: its ticks hang below it
    F = nxr_figure(o.Width, 34 + 8 * nL, {'A'}, 'Kinds', struct('A', 'plot'), 'Keys', struct('A', keys));
    hFig = F.fig;  hFig.Visible = o.Visible;
    ax = F.ax.A;  hold(ax, 'on');
    P = ax.Position;  ax.Position = [P(1) + 6 / F.size_mm(1), P(2), P(3) - 36 / F.size_mm(1), P(4)];   % room: level names, z-stack
    G = struct('level', {}, 't0', {}, 't1', {});
    if hasField
        lim = [min(cellfun(@min, o.Field)) max(cellfun(@max, o.Field))];
        cm = nxr_cmap(T.quantity.(o.Quantity).colormap);
    end
    for i = 1:nL
        d = base * 2^levels(i);  n = ceil(duration / d - 1e-9);
        t0 = (0:n-1) * d;  t1 = min(t0 + d, duration);
        G(i) = struct('level', levels(i), 't0', t0, 't1', t1);
        y = i - 0.42;  hgt = 0.84;
        if hasField
            v = double(o.Field{i}(:)');
            if numel(v) ~= n, error('rheome:show:levelstack', 'Field{%d}: %d tiles at level %d, got %d values.', i, n, levels(i), numel(v)); end
            c = cm(1 + round((min(max(v, lim(1)), lim(2)) - lim(1)) / max(diff(lim), eps) * (size(cm, 1) - 1)), :);
            for k = 1:n, patch(ax, [t0(k) t1(k) t1(k) t0(k)], y + [0 0 hgt hgt], c(k, :), 'EdgeColor', 'none'); end
        end
        on = ~isempty(o.Slice) && levels(i) == o.Slice;
        ec = T.rgb.neutral_slate;  lw = T.line.outline_pt;
        if on, ec = T.rgb.neutral_ink;  lw = 2.5 * T.line.outline_pt; end
        for k = 1:n
            plot(ax, [t0(k) t1(k) t1(k) t0(k) t0(k)], y + [0 0 hgt hgt 0], '-', 'Color', ec, 'LineWidth', lw);
        end
    end
    if ~isempty(o.Tile)
        i = find(levels == o.Tile(1), 1);
        if isempty(i) || o.Tile(2) < 1 || o.Tile(2) > numel(G(i).t0), error('rheome:show:levelstack', 'Tile: [level index] of a tile shown.'); end
        a = G(i).t0(o.Tile(2));  b = G(i).t1(o.Tile(2));
        patch(ax, [a b b a], i - 0.42 + [0 0 0.84 0.84], T.rgb.atlas_time, 'FaceAlpha', 0.35 * ~hasField, 'EdgeColor', T.rgb.atlas_time, 'LineWidth', 2.5 * T.line.outline_pt);
    end
    ax.XLim = [0 duration];  ax.YLim = [0.4 nL + 0.6];
    ax.YTick = 1:nL;  ax.YTickLabel = compose('%d', levels);
    nxr_axes_style(ax, 'time (s)', 'level j (tile = 2^j base)');
    ax.YLabel.Interpreter = 'tex';
    if ~isempty(o.SpaceMM)
        yyaxis(ax, 'right');  ax.YLim = [0.4 nL + 0.6];  ax.YTick = 1:nL;
        ax.YTickLabel = compose('%.3g', o.SpaceMM);  ax.YColor = T.rgb.neutral_ink;
        ylabel(ax, 'spatial tile (mm)', 'FontSize', T.type.axisTitle.size_pt);
        yyaxis(ax, 'left');
    end
    hold(ax, 'off');
    if ~isempty(o.Slice)   % the z-stack indicator, right of the grid
        Q = ax.Position;  w = 12 / F.size_mm(1);
        zi = axes(hFig, 'Position', [Q(1) + Q(3) + 6 / F.size_mm(1) + (~isempty(o.SpaceMM)) * 12 / F.size_mm(1), Q(2) + Q(4) / 2 - w, w, 2 * w]);
        rheome.show.zstack(zi, nL, find(levels == o.Slice));
    end
    e = {struct('kind', 'patch', 'text', sprintf('tile of a level (%g s %s 2^j)', base, char(215)), 'color', T.rgb.neutral_paper)};
    if ~isempty(o.Slice), e{end+1} = struct('kind', 'line', 'text', sprintf('slice: level %d', o.Slice), 'color', T.rgb.neutral_ink); end
    if ~isempty(o.Tile), e{end+1} = struct('kind', 'patch', 'text', sprintf('tile %d of level %d', o.Tile(2), o.Tile(1)), 'color', T.rgb.atlas_time); end
    nxr_key(F, 'A', e);
    if hasField   % three ticks: the ends and the middle (the kit's default adds MATLAB's, which crowd a short bar)
        pos = ax.Position;   % the kit's southoutside bar shrinks the axes before it moves the bar to the key band
        cb = nxr_colorbar(ax, o.Label, o.Unit, o.Quantity, lim, 'F', F);  ax.Position = pos;
        tk = linspace(lim(1), lim(2), 3);  cb.Ticks = tk;  cb.TickLabels = compose('%.3g', tk);
    end
    finish(F, o.Title, o.Export);
end
