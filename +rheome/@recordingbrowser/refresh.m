function refresh(app)
% REFRESH  Pick the level for the window, derive its statistics, redraw.
%
% The whole level-of-detail loop is here and nowhere else: every setter calls this, so the
% level, the table, the strip and the three axes cannot disagree about what is on screen.
%
% ⭐ ONE LEVEL PER FRAME. rheome.select.derive reads only the moments the chosen statistic needs,
% and rheome.select.level caches each level's arrays in the handle, so panning inside a level costs
% arithmetic on a few hundred rows and no disk at all.
%
% ⚠ REENTRANCY. Drawing sets axis limits, and an axis-limit listener calls setWindow; the
% Busy flag breaks that loop. Without it a zoom recurses until the stack gives out.
%
% Author: Diellor Basha, 2026

    if app.Busy, return; end
    app.Busy = true;
    c = onCleanup(@() i_clear(app));

    g = app.Db.grid;
    env = app.Envelope || strcmp(app.View, 'channels');        % the stack is min/max only
    [fl, why] = rheome.recordingbrowser.floorlevel(app.Db, app.Scope, app.Band, app.Stat, env);
    app.Floor = fl;  app.FloorReason = why;
    w = diff(app.Window);
    bud = app.budget();
    if app.Auto
        app.Level = rheome.recordingbrowser.levelfor(app.Db, w, bud, fl);
    else
        app.Level = min(max(app.Level, fl), g.Lmax);
    end
    L = app.Level;

    if strcmp(app.View, 'channels')
        app.Stack = i_stack(app, L);
        app.Rows = table();  app.Strip = struct('mode', 'none', 'level', L, 'bands', zeros(1,0), ...
                                                'levels', zeros(1,0), 'rows', struct([]), 'hasBand', false);
        app.Raw = i_raw(app);
        if isempty(app.Fig) || ~isvalid(app.Fig), return; end
        app.drawChannels();  app.drawRaw();  app.updateLabels();
        return
    end

    % ---- the rows on screen ----
    % In envelope mode the level can sit below the channel floor, where only min and max
    % exist; asking for anything else there would fail, so the top axis asks for exactly
    % what it draws and the strip fetches its own rows at a level that carries bands.
    if app.Envelope
        stats = ["min", "max"];
    elseif isempty(app.Db.bands)
        stats = string(app.Stat);                          % a preview store: no band columns
    else
        stats = unique([string(app.Stat), "bandPower"], 'stable');
    end
    R = i_derive(app, L, stats, []);
    app.Rows = R;
    if isempty(app.Ladder) && ~isempty(app.Db.bands)
        app.Ladder = rheome.select.ladder(app.Db, Scope=string(app.Scope));
    end

    % ---- the strip: one row per band, at the current level or at each band's own ----
    app.Strip = i_strip(app, R);

    % ---- the raw axis: phase 2, and only when the window is short enough ----
    app.Raw = i_raw(app);

    if isempty(app.Fig) || ~isvalid(app.Fig), return; end
    app.drawStat();
    app.drawTiles();
    app.drawRaw();
    app.updateLabels();
end

function i_clear(app)
    if isvalid(app), app.Busy = false; end
end

% ⭐ ONE DERIVE CALL FOR THE WHOLE STACK. min and max for every selected channel come out of
% the same level array, so drawing thirty-two channels costs what drawing one costs plus the
% arithmetic; the level is read once and cached in the handle.
function S = i_stack(app, L)
    u = app.Units(:)';
    u = u(1:min(numel(u), app.MaxRows));
    if strcmp(app.Scope, 'group')
        R = rheome.select.derive(app.Db, Level=L, Stats=["min","max"], Scope="group", Groups=u, Window=app.Window);
    else
        R = rheome.select.derive(app.Db, Level=L, Stats=["min","max"], Channels=u, Window=app.Window);
    end
    k = unique(R.k, 'stable');
    nK = numel(k);
    S = struct('units', u, 'k', k, 'tLo', R.t_lo(1:nK), 'tHi', R.t_hi(1:nK), ...
               'min', zeros(nK, numel(u)), 'max', zeros(nK, numel(u)));
    for i = 1:numel(u)
        r = R(R.unit_id == u(i), :);
        S.min(:, i) = r.min;  S.max(:, i) = r.max;
    end
end

% One derive call for the current unit, whatever the scope. Bands empty = all carried.
function R = i_derive(app, L, stats, bandSel)
    if strcmp(app.Scope, 'group')
        R = rheome.select.derive(app.Db, Level=L, Stats=stats, Scope="group", Groups=app.Unit, ...
                          Bands=bandSel, Window=app.Window);
    else
        R = rheome.select.derive(app.Db, Level=L, Stats=stats, Channels=app.Unit, ...
                          Bands=bandSel, Window=app.Window);
    end
end

% The strip always needs a level that carries bands, which in channel scope is the channel
% floor even when the envelope above it is drawn finer. Its own level is reported in the
% caption, so a strip coarser than the axis above never reads as the same tiling.
%
% ⭐ THE STRIP IS WHERE THE CONSTANT-Q TILING BECOMES VISIBLE. In 'natural' mode each band
% is drawn at its OWN tile length -- the first level at least its support long, or as fine
% as the budget and the channel floor allow -- so a row at 64-128 Hz is cut into 0.25 s
% tiles while a row at 1-2 Hz is cut into 16 s tiles, in the same picture. That staircase IS
% the store's diagonal. In 'level' mode every row is cut at the level the statistic uses, so
% all the columns line up with the axis above and nothing mixes.
function S = i_strip(app, R)
    g = app.Db.grid;  b = app.Db.bands;
    base = 0;
    if ~strcmp(app.Scope, 'group'), base = app.Db.channelLevel; end
    S = struct('mode', app.StripMode, 'level', app.Level, 'hasBand', false, ...
               'bands', [], 'levels', [], 'rows', struct('band', {}, 'level', {}, 'tLo', {}, 'tHi', {}, 'power', {}), ...
               'power', [], 'tLo', [], 'tHi', []);
    if isempty(b)
        S.level = app.Level;  S.bands = zeros(1, 0);  S.levels = zeros(1, 0);
        return                                             % a preview store has no strip
    end
    if strcmp(app.StripMode, 'level')
        Ls = max(app.Level, max(base, min(b.naturalLevel)));
        if Ls == app.Level && ismember('bandPower', R.Properties.VariableNames)
            Rs = R;                                            % already in hand
        else
            Rs = i_derive(app, Ls, "bandPower", []);
        end
        S.level = Ls;
        S.bands = Rs.Properties.CustomProperties.bands;
        S.levels = repmat(Ls, 1, numel(S.bands));
        S.power = Rs.bandPower';  S.tLo = Rs.t_lo;  S.tHi = Rs.t_hi;
        for i = 1:numel(S.bands)
            S.rows(i) = struct('band', S.bands(i), 'level', Ls, 'tLo', Rs.t_lo, 'tHi', Rs.t_hi, 'power', Rs.bandPower(:, i)');
        end
    else
        w = diff(app.Window);
        js = b.j(:)';
        for i = 1:numel(js)
            k = find(b.j == js(i), 1);
            Lj = rheome.recordingbrowser.levelfor(app.Db, w, app.budget(), max(b.naturalLevel(k), base));
            Rj = i_derive(app, Lj, "bandPower", js(i));
            S.rows(i) = struct('band', js(i), 'level', Lj, 'tLo', Rj.t_lo, 'tHi', Rj.t_hi, 'power', Rj.bandPower(:)');
        end
        S.bands = js;  S.levels = [S.rows.level];
    end
    S.hasBand = ismember(app.Band, S.bands);
end

% ⚠ NOTHING IS READ HERE. Raw samples arrive only when asked for (recordingbrowser/loadRaw,
% the load button, key l): a read pulls every channel's hyperslab, so 8 s of one sensor
% moves 10 MB and 2 minutes moves 156 MB, and a browser that did that on every zoom would
% stall on the way past. What refresh does is decide what the panel SAYS about the samples
% already in hand, and whether they still cover the window.
function raw = i_raw(app)
    raw = app.Raw;
    if ~isfield(raw, 'span'), raw.span = [NaN NaN]; end
    if ~isfield(raw, 'unit'), raw.unit = app.Unit; end
    w = app.Window;
    if ~app.ShowRaw
        raw.note = 'raw panel off (r)';  return
    end
    if strcmp(app.Scope, 'group')
        raw.note = 'raw samples are per channel: switch Scope to channel';  return
    end
    if ~isempty(raw.x) && raw.span(1) <= w(1) + 1e-9 && raw.span(2) >= w(2) - 1e-9
        raw.note = sprintf('%s: %d samples loaded for %.1f-%.1f s', ...
            i_uname(app, raw.unit), numel(raw.x), raw.span(1), raw.span(2));
        if raw.unit ~= app.Unit
            raw.note = sprintf('%s  (focus is %s: press load for it)', raw.note, i_uname(app, app.Unit));
        end
        return
    end
    n = round(diff(w) * app.Db.meta.fs);
    if diff(w) > app.RawSeconds
        raw.note = sprintf('%.0f s window: narrow it below %g s, then load (%d samples, %.0f MB off disk)', ...
            diff(w), app.RawSeconds, n, n * app.Db.meta.C * 8 / 1e6);
    elseif isempty(raw.x)
        raw.note = sprintf('press load (l) to read %d samples of %s for this window', n, i_uname(app, app.Unit));
    else
        raw.note = sprintf('%s loaded for %.1f-%.1f s; press load (l) for this window', ...
            i_uname(app, raw.unit), raw.span(1), raw.span(2));
    end
end

function nm = i_uname(app, u)
    if strcmp(app.Scope, 'group') || isempty(u), nm = sprintf('node %d', u);
    else, nm = char(string(app.Db.meta.ChannelName{u})); end
end

% Author: Diellor Basha, 2026
