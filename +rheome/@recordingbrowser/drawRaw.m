function drawRaw(app)
% DRAWRAW  Phase 2: the actual samples, once the window is short enough to justify reading them.
%
% ⚠ THE NOTE IS PART OF THE ANSWER. When the axis is empty it says why -- nothing loaded
% yet, the window is wider than loadRaw will read, group scope, or no recording store on
% this disk -- because an empty trace otherwise reads as "no signal here", which is the
% opposite of what it means.
%
% ⭐ SAMPLES ARRIVE ONLY THROUGH loadRaw. Panning and zooming never read the recording, so
% what is drawn here is whatever was last asked for, and the note says which span that is
% once the window has moved off it.
%
% Author: Diellor Basha, 2026

    ax = app.RawAx;  cla(ax, 'reset');
    if isempty(app.Raw.x)
        text(ax, 0.5, 0.5, app.Raw.note, 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45], 'Interpreter', 'none');
        set(ax, 'XTick', [], 'YTick', []);
        xlim(ax, app.Window);
        return
    end
    plot(ax, app.Raw.t, app.Raw.x, 'Color', [0.25 0.25 0.25]);
    hold(ax, 'on');
    R = app.Rows;
    for k = 1:height(R)
        xline(ax, R.t_lo(k), 'Color', [0.75 0.75 0.75]);        % the tile edges, for scale
    end
    grid(ax, 'on');
    xlim(ax, app.Window);
    ylabel(ax, app.Db.meta.units);  xlabel(ax, 'time (s)');
    title(ax, sprintf('raw samples  --  %s', app.Raw.note), 'Interpreter', 'none', 'FontSize', 9);
end
% Author: Diellor Basha, 2026
