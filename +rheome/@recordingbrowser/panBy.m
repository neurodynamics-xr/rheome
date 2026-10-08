function panBy(app, frac)
% PANBY  Slide the window by a fraction of its width (negative = earlier).
% Author: Diellor Basha, 2026
    d = double(frac) * diff(app.Window);
    app.setWindow(app.Window(1) + d, app.Window(2) + d);
end
% Author: Diellor Basha, 2026
