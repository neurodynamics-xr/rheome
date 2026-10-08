function setLevel(app, L)
% SETLEVEL  Pin the level, leaving Auto off until setAuto(true) turns the rule back on.
%
% A pinned level is still clamped to the diagonal's floor and to Lmax -- the store holds
% nothing below the floor, so pinning there would ask for rows that do not exist.
%
% Author: Diellor Basha, 2026

    g = app.Db.grid;
    [fl, why] = rheome.recordingbrowser.floorlevel(app.Db, app.Scope, app.Band, app.Stat);
    app.Floor = fl;  app.FloorReason = why;
    Lw = min(max(round(double(L)), fl), g.Lmax);
    app.Auto = false;
    app.Level = Lw;
    app.refresh();
    if Lw ~= round(double(L))
        app.setStatus(sprintf('level %d is below the floor: %s', round(double(L)), why));
    end
end
% Author: Diellor Basha, 2026
