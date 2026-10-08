function setAuto(app, on)
% SETAUTO  Hand the level back to the level-of-detail rule (or take it away).
% Author: Diellor Basha, 2026
    app.Auto = logical(on);
    app.refresh();
end
% Author: Diellor Basha, 2026
