function onMode(app, m)
% ONMODE  Signed vorticity (CCW +/CW -) vs magnitude (the envelope). See @flowpage/map.
% Author: Diellor Basha, 2026
    app.MapMode = m;
    app.drawCortex();
end

% Author: Diellor Basha, 2026
