function [hFig, G] = levelstack(duration, levels, varargin)
% SHOW.LEVELSTACK  The volume of the dyadic levels: time x level, each level a row of tiles of halving width.
%
%   hFig = rheome.show.levelstack(duration, levels)
%   hFig = rheome.show.levelstack(duration, 0:5, 'Slice', 3, 'Tile', [3 5], 'SpaceMM', s)
%   hFig = rheome.show.levelstack(duration, 0:5, 'Field', vals, 'Quantity', 'power', 'Label', 'power', 'Unit', 'dB')
%
% A level is not a measurement: it is the spacing of the grid on which measurements are made, as the grid spacing is
% in Euclidean space (the step of the differential operators, the support of a tile). So a level never gets a colour
% (tokens levels.rule). Level l tiles the time shown into 2^l tiles of duration / 2^l; each tile of level l is two of
% level l + 1, and the rows stacked (coarse at the bottom, finer upward) form the constant-Q volume. Each level is
% stated by its spacing on the axis ("3 (8 x 0.5 s)"); over the cortex a level also has its spatial tile ('SpaceMM',
% right axis). One level is a SLICE through the stack (an ink frame; the z-stack indicator beside the volume,
% rheome.show.zstack), one tile a cell of it ('Tile': it and the tiles holding it, the roll-up, outlined in ink). A
% FIELD measured on the tiles (a scalogram: a value per tile per level) fills the cells on its quantity's own map with
% a colourbar: the field is coloured, the grid is not. Same drawing as style.level_scalogram (Python); an idea for
% review with Science (tokens levels.volume).
%
% INPUTS:
%   duration  s, the time shown (from 0)
%   levels    increasing integer levels >= 0 (e.g. 0:6)
%
% OPTIONS (name/value):
%   'Slice'     a level drawn as the current slice (default none)
%   'Tile'      [level index] (index 1-based along time): the tile and its ancestors outlined (default none)
%   'SpaceMM'   [1 x nL] the spatial tile at each level, mm (the cortical pyramid's), on the right axis
%   'Field'     {1 x nL} 2^l values for level l (a scalogram on the grid); with Quantity, Label and Unit
%   'Quantity'  the field's quantity (figure tokens quantity: power, magnitude, divergence ...; default 'power')
%   'Label','Unit'  the field's name and unit on the colourbar (required with Field)
%   'Width'     'single' | 'onehalf' | 'double' (default 'double')
%   'Export','Visible','Title'  as every rheome.show figure
%
% OUTPUT:
%   hFig the figure;  G struct per level: .level .t0 .t1 (s, the tiles' edges)
%
% Example:
%   rheome.show.levelstack(4, 0:5, 'Slice', 3, 'Tile', [5 11], 'SpaceMM', 122 ./ 2.^(0:5));
%
% See also: rheome.show.zstack, rheome.show.timetracks
%
% Author: Diellor Basha, 2026
    T = tokens();
    p = inputParser;
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
    if isempty(levels) || any(levels < 0 | levels ~= round(levels)) || any(diff(levels) <= 0)
        error('rheome:show:levelstack', 'levels: increasing integers >= 0.');
    end
    nL = numel(levels);
    hasField = ~isempty(o.Field);
    if hasField && (numel(o.Field) ~= nL || isempty(o.Label) || isempty(o.Unit))
        error('rheome:show:levelstack', 'Field: one vector per level, with Label and Unit.');
    end
    if ~isempty(o.SpaceMM) && numel(o.SpaceMM) ~= nL, error('rheome:show:levelstack', 'SpaceMM: one size per level.'); end
    if ~isempty(o.Slice) && ~ismember(o.Slice, levels), error('rheome:show:levelstack', 'Slice: one of the levels.'); end
    if ~isempty(o.Tile) && (~ismember(o.Tile(1), levels) || o.Tile(2) < 1 || o.Tile(2) > 2^o.Tile(1) || o.Tile(2) ~= round(o.Tile(2)))
        error('rheome:show:levelstack', 'Tile: [level index], index in 1..2^level.');
    end

    keys = 'g';  if hasField, keys = 'gs'; end   % the colour row last: its ticks hang below it
    F = nxr_figure(o.Width, 34 + 8 * nL, {'A'}, 'Kinds', struct('A', 'plot'), 'Keys', struct('A', keys));
    hFig = F.fig;  hFig.Visible = o.Visible;
    ax = F.ax.A;  hold(ax, 'on');
    P = ax.Position;  ax.Position = [P(1) + 10 / F.size_mm(1), P(2), P(3) - 40 / F.size_mm(1), P(4)];   % level names; z-stack
    if hasField
        lim = [min(cellfun(@min, o.Field)) max(cellfun(@max, o.Field))];
        cm = nxr_cmap(T.quantity.(o.Quantity).colormap);
    end
    G = struct('level', {}, 't0', {}, 't1', {});
    for i = 1:nL
        n = 2^levels(i);  d = duration / n;
        t0 = (0:n-1) * d;  t1 = t0 + d;
        G(i) = struct('level', levels(i), 't0', t0, 't1', t1);
        y = i - 1;
        if hasField
            v = double(o.Field{i}(:)');
            if numel(v) ~= n, error('rheome:show:levelstack', 'Field{%d}: %d tiles at level %d, got %d values.', i, n, levels(i), numel(v)); end
            c = cm(1 + round((min(max(v, lim(1)), lim(2)) - lim(1)) / max(diff(lim), eps) * (size(cm, 1) - 1)), :);
            ec = T.rgb.neutral_paper;
        else
            c = repmat(T.rgb.neutral_paper, n, 1);  ec = T.rgb.neutral_mist;
        end
        for k = 1:n
            patch(ax, [t0(k) t1(k) t1(k) t0(k)], y + [0 0 1 1], c(k, :), 'EdgeColor', ec, 'LineWidth', T.line.outline_pt);
        end
    end
    if ~isempty(o.Slice)
        i = find(levels == o.Slice);
        plot(ax, [0 duration duration 0 0], i - 1 + [0 0 1 1 0], '-', 'Color', T.rgb.neutral_ink, 'LineWidth', 3 * T.line.outline_pt);
    end
    if ~isempty(o.Tile)   % the tile and the tiles holding it (its roll-up)
        l = o.Tile(1);  k = o.Tile(2) - 1;
        for j = levels(levels <= l)
            i = find(levels == j);  d = duration / 2^j;  kk = floor(k / 2^(l - j));
            plot(ax, kk * d + [0 d d 0 0], i - 1 + [0 0 1 1 0], '-', 'Color', T.rgb.neutral_ink, 'LineWidth', 2 * T.line.outline_pt);
        end
    end
    ax.XLim = [0 duration];  ax.YLim = [0 nL];
    ax.YTick = (1:nL) - 0.5;
    ax.YTickLabel = arrayfun(@(l) sprintf('%d (%d %s %.3g s)', l, 2^l, char(215), duration / 2^l), levels, 'UniformOutput', false);
    nxr_axes_style(ax, 'time (s)', 'level (tiles)');
    if ~isempty(o.SpaceMM)
        yyaxis(ax, 'right');  ax.YLim = [0 nL];  ax.YTick = (1:nL) - 0.5;
        ax.YTickLabel = compose('%.3g', o.SpaceMM);  ax.YColor = T.rgb.neutral_ink;
        ylabel(ax, 'spatial tile (mm)', 'FontSize', T.type.axisTitle.size_pt);
        yyaxis(ax, 'left');
    end
    hold(ax, 'off');
    if ~isempty(o.Slice)   % the z-stack indicator, right of the volume
        Q = ax.Position;  w = 14 / F.size_mm(1);
        zi = axes(hFig, 'Position', [Q(1) + Q(3) + (6 + 12 * ~isempty(o.SpaceMM)) / F.size_mm(1), Q(2) + Q(4) / 4, w, Q(4) / 2]);
        rheome.show.zstack(zi, levels, o.Slice);
    end
    e = {struct('kind', 'patch', 'text', sprintf('tile: level l has 2^l tiles of %.3g s / 2^l', duration), 'color', T.rgb.neutral_paper)};
    if ~isempty(o.Slice), e{end+1} = struct('kind', 'line', 'text', sprintf('slice: level %d', o.Slice), 'color', T.rgb.neutral_ink); end
    if ~isempty(o.Tile), e{end+1} = struct('kind', 'isoline', 'text', sprintf('tile %d of level %d and the tiles holding it', o.Tile(2), o.Tile(1)), 'color', T.rgb.neutral_ink); end
    nxr_key(F, 'A', e);
    if hasField   % three ticks: the ends and the middle
        pos = ax.Position;   % the kit's southoutside bar shrinks the axes before it moves the bar to the key band
        cb = nxr_colorbar(ax, o.Label, o.Unit, o.Quantity, lim, 'F', F);  ax.Position = pos;
        tk = linspace(lim(1), lim(2), 3);  cb.Ticks = tk;  cb.TickLabels = compose('%.3g', tk);
    end
    finish(F, o.Title, o.Export);
end
