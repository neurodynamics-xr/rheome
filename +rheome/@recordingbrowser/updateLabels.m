function updateLabels(app)
% UPDATELABELS  The readout and the controls, kept in step with the state.
%
% The readout names the window, the level, the tile length, how many tiles are on screen and
% what set the floor -- everything needed to know whether what is drawn is the finest the
% store can answer or only the finest the budget allows.
%
% Author: Diellor Basha, 2026

    if isempty(app.Fig) || ~isvalid(app.Fig), return; end
    g = app.Db.grid;
    i_fileline(app);
    if app.Auto, mode = 'auto'; else, mode = 'pinned'; end
    line1 = sprintf('%s store   %.1f-%.1f s (%.1f s)   L%d %s, %g s tiles, %d of %d in view   floor L%d: %s', ...
        app.StoreKind, app.Window(1), app.Window(2), diff(app.Window), app.Level, mode, ...
        g.tExtent(app.Level+1), i_ntiles(app), app.budget(), app.Floor, app.FloorReason);
    set(app.LevelLbl, 'String', {line1, i_bandline(app)});
    if ~isempty(app.LevelSpin) && isvalid(app.LevelSpin)
        app.syncLevels();                                  % the floor moves with stat and scope
        set(app.LevelSpin, 'Enable', i_onoff(~app.Auto));
    end
    if ~isempty(app.AutoBox)  && isvalid(app.AutoBox),  set(app.AutoBox,  'Value', app.Auto); end
    if ~isempty(app.EnvBox)   && isvalid(app.EnvBox),   set(app.EnvBox,   'Value', app.Envelope); end
    if ~isempty(app.ScopeDrop) && isvalid(app.ScopeDrop), set(app.ScopeDrop, 'Value', 1 + strcmp(app.Scope, 'group')); end
    if ~isempty(app.StatDrop) && isvalid(app.StatDrop)
        names = rheome.recordingbrowser.statsfor(app.Db);
        i = find(names == string(app.Stat), 1);
        if ~isempty(i), set(app.StatDrop, 'Value', i); end
    end
    if ~isempty(app.BandDrop) && isvalid(app.BandDrop) && ~isempty(app.Db.bands)
        set(app.BandDrop, 'Value', find(app.Db.bands.j == app.Band, 1));
    end
    if ~isempty(app.UnitDrop) && isvalid(app.UnitDrop)
        i = find(app.UnitIds_ == app.Unit, 1);
        if ~isempty(i), set(app.UnitDrop, 'Value', i); end
    end
    if ~isempty(app.ViewDrop) && isvalid(app.ViewDrop)
        set(app.ViewDrop, 'Value', 1 + strcmp(app.View, 'detail'));
    end
    if ~isempty(app.ChanList) && isvalid(app.ChanList) && ~isempty(app.UnitIds_)
        sel = find(ismember(app.UnitIds_, app.Units));
        if ~isequal(sort(get(app.ChanList, 'Value')), sort(sel)), set(app.ChanList, 'Value', sel); end
    end
end

% The band's own row of the ladder: what this band costs in time and in rate. This is the
% answer to "which band am I working with, and at what tile and rate" without leaving the
% window.
function n = i_ntiles(app)
    if strcmp(app.View, 'channels'), n = numel(app.Stack.k); else, n = height(app.Rows); end
end

% ⭐ THE HEADER NAMES THE FILE. A view of a stored matrix that does not say which matrix is
% a screenshot waiting to be misattributed: the store, its size, the recording it summarises,
% and the shape of that recording.
function i_fileline(app)
    if isempty(app.FileLbl) || ~isvalid(app.FileLbl), return; end
    m = app.Db.meta;
    d = dir(app.Db.file);  mb = 0;  if ~isempty(d), mb = d.bytes/1e6; end
    [~, stem, ext] = fileparts(app.Db.file);
    src = '';  srcmb = 0;
    if isfield(m, 'source') && ~isempty(m.source)
        [~, s2, e2] = fileparts(char(m.source));  src = [s2 e2];
        d2 = dir(char(m.source));  if ~isempty(d2), srcmb = d2.bytes/1e6; end
    end
    set(app.FileLbl, 'String', sprintf('%s%s  (%s store, %.0f MB)   <-  %s (%.0f MB)   %d ch x %d samples @ %g Hz = %.0f s', ...
        stem, ext, app.StoreKind, mb, src, srcmb, m.C, m.nT, m.fs, m.duration));
end

function s = i_bandline(app)
    if isempty(app.Db.bands)                               % a preview store: the window's support
        T = rheome.ingest.support(diff(app.Window));
        s = sprintf('preview store, %g s tiles: this %.3g s window supports %.3g Hz and above (lowest full octave %.3g-%.3g Hz, %.0f cycles)', ...
            app.Db.grid.tExtent(app.Level+1), T.window_s, T.f_lowest, T.band_lo, T.band_hi, T.cycles_in_window);
        return
    end
    T = app.Ladder;
    if isempty(T), s = ''; return; end
    r = T(T.band_id == app.Band, :);
    if isempty(r), s = ''; return; end
    s = sprintf('band %d  %.4g-%.4g Hz  Q %.1f  support %.3g s  ->  %g s tile (L%d), %.1f cycles, %.3g Hz rate, %.0f samples per tile (%.0f at full rate)', ...
        r.band_id, r.f_lo, r.f_hi, r.q, r.support_s, r.tile_s, r.natural_level, ...
        r.cycles_per_tile, r.rate_hz, r.samples_per_tile, r.samples_full);
    if ismember('floor_level', T.Properties.VariableNames) && r.floor_level > r.natural_level
        s = sprintf('%s   [this sensor: %g s tiles at L%d]', s, r.floor_tile_s, r.floor_level);
    end
end

function v = i_onoff(tf)
    if tf, v = 'on'; else, v = 'off'; end
end
% Author: Diellor Basha, 2026
