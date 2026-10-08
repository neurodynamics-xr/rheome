function updateCursors(app)
% UPDATECURSORS  Move the timeline cursors without redrawing the traces.
% Author: Diellor Basha, 2026

    t = app.Page.Time(app.Frame);
    hs = [app.TimeCursors, app.SensorCursor];
    for h = hs
        if ~isempty(h) && isvalid(h), h.Value = t; end
    end
end

% Author: Diellor Basha, 2026
