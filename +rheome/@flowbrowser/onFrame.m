function onFrame(app, i)
% ONFRAME  Scrub time. The fast path -- a GEMV per scale, no recompute.
% Author: Diellor Basha, 2026
    i = min(max(round(i), 1), app.Page.NumFrames);
    if i == app.Frame, return; end
    app.Frame = i;
    app.refresh();
end

% Author: Diellor Basha, 2026
