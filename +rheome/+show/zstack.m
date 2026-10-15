function h = zstack(ax, levels, shown)
% SHOW.ZSTACK  The level indicator: the dyadic stack as rows of 2^l cells, the shown level the ink row, numbered.
%
%   h = rheome.show.zstack(ax, levels, shown)        e.g. rheome.show.zstack(ax, 0:5, 3)
%
% Beside a cortical map or a scalogram slice it says WHICH level of the constant-Q stack is shown, without a colour
% (a level is the grid spacing, not a field: tokens levels.slice). Coarse on top; each row is 2^l cells outlined in
% mist; the shown row is outlined in ink on snow and labelled "l = <level>". The same glyph as style.level_indicator
% (Python). ax's axes are hidden.
%
% See also: rheome.show.levelstack
%
% Author: Diellor Basha, 2026
    T = tokens();
    levels = double(levels(:)');
    if isempty(levels) || any(levels < 0 | levels ~= round(levels)) || ~ismember(shown, levels)
        error('rheome:show:zstack', 'levels: integers >= 0; shown: one of them.');
    end
    hold(ax, 'on');
    h = gobjects(0);
    lo = 0.12;  hi = 0.88;   % a cell's bottom and top within its row
    for i = 1:numel(levels)
        nc = 2^levels(i);  on = levels(i) == shown;
        for c = 0:nc-1
            if on
                h(end+1) = patch(ax, [c c+1 c+1 c] / nc, i - 1 + [lo lo hi hi], T.rgb.neutral_snow, 'EdgeColor', T.rgb.neutral_ink, 'LineWidth', T.line.outline_pt); %#ok<AGROW>
            else
                h(end+1) = patch(ax, [c c+1 c+1 c] / nc, i - 1 + [lo lo hi hi], T.rgb.neutral_paper, 'EdgeColor', T.rgb.neutral_mist, 'LineWidth', 0.75 * T.line.outline_pt); %#ok<AGROW>
            end
        end
        if on
            h(end+1) = text(ax, 1.06, i - 0.5, sprintf('%s = %d', char(8467), shown), 'FontSize', T.type.annotation.size_pt, ...
                'Color', T.rgb.neutral_ink, 'VerticalAlignment', 'middle', 'Interpreter', 'none'); %#ok<AGROW>
        end
    end
    ax.YDir = 'reverse';  ax.XLim = [0 1];  ax.YLim = [0 numel(levels)];  axis(ax, 'off');
    hold(ax, 'off');
end
