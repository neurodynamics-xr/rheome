function refreshMarks(app)
% REFRESHMARKS  Re-read the measurements from the sidecar and redraw.
%
% ⚠ READ ONCE, NOT PER FRAME. The measurements live in a file beside the store, so reading
% them on every redraw would put a disk read in the middle of a pan. They are read here --
% at construction and after a write -- and filtered to the window when drawn.
%
% Author: Diellor Basha, 2026

    try
        app.Marks = rheome.select.measures(app.Db);
    catch
        app.Marks = table();
    end
    if ~isempty(app.Fig) && isvalid(app.Fig) && ~app.Busy
        app.refresh();
    end
end
% Author: Diellor Basha, 2026
