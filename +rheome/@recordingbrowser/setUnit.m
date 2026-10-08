function setUnit(app, u)
% SETUNIT  Choose the channel (or sensor-tree node) the view reads. Accepts a channel name.
% Author: Diellor Basha, 2026
    app.Unit = rheome.recordingbrowser.unitid(app.Db, app.Scope, u);
    app.refresh();     % the loaded samples stay: they carry their own channel (loadRaw)
end
% Author: Diellor Basha, 2026
