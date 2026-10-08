function shadeMargin(app, ax)
% SHADEMARGIN  Grey the frames this page does NOT own.
% Author: Diellor Basha, 2026

    t = app.Page.Time;  c = app.Page.Core;
    if all(c), return; end
    yl = ylim(ax);
    hold(ax, 'on');
    i1 = find(c, 1);  i2 = find(c, 1, 'last');
    col = [0.86 0.86 0.86];
    if i1 > 1
        patch(ax, [t(1) t(i1) t(i1) t(1)], [yl(1) yl(1) yl(2) yl(2)], col, ...
            'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
    end
    if i2 < numel(t)
        patch(ax, [t(i2) t(end) t(end) t(i2)], [yl(1) yl(1) yl(2) yl(2)], col, ...
            'EdgeColor','none','FaceAlpha',0.45,'HandleVisibility','off');
    end
    hold(ax, 'off');
end

% Author: Diellor Basha, 2026
