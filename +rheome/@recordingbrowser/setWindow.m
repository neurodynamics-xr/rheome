function setWindow(app, t0, t1)
% SETWINDOW  Move the visible span, clamped to the record, and re-pick the level.
%
%   app.setWindow(120, 180)
%   app.setWindow([120 180])
%
% ⚠ CLAMPED, NEVER THROWN. A drag, a key repeat or a zoom at the edge all overshoot; a
% browser that errors in a callback for that is unusable. The span keeps its width when it
% is pushed against an edge, and it never gets shorter than one level-0 tile.
%
% Author: Diellor Basha, 2026

    if nargin == 2, v = double(t0(:)'); t0 = v(1); t1 = v(2); end
    D = app.Db.meta.duration;  fmin = app.Db.grid.tExtent(1);
    t0 = double(t0);  t1 = double(t1);
    if t1 < t0, [t0, t1] = deal(t1, t0); end
    w = max(t1 - t0, fmin);
    if w >= D
        t0 = 0;  w = D;
    else
        t0 = min(max(t0, 0), D - w);
    end
    app.Window = [t0 t0 + w];
    app.refresh();
end
% Author: Diellor Basha, 2026
