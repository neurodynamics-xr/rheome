classdef recordingbrowser < handle
% RECORDINGBROWSER  Look at a stored recording through its pyramid, before loading any of it.
%
%   app = rheome.recordingbrowser('sub01')                        % the best pyramid it has
%   app = rheome.recordingbrowser('sub01', 'Store', 'preview')    % the min/max pyramid
%   app = rheome.recordingbrowser(storeFile, 'Units', 1:16)               % those channels, stacked
%   app = rheome.recordingbrowser(db, 'View', 'detail', 'Unit', 'MLO32')  % one channel, in full
%
% ⭐ WHAT IT IS. A stored recording is a matrix too big to draw. The pyramid holds, for every
% channel and every tile, the minimum and the maximum of the samples under it -- and those
% two merge by extremum, so the pair at any level BOUNDS every sample it covers. Drawing them
% as a band is therefore an exact picture of the matrix at any zoom: too wide at worst, never
% too narrow, and no transient can hide between the drawn columns. That is the whole app.
% The channels view stacks those bands, one row per channel, so you can see where in the
% recording and in which sensors something happens, and only then read samples.
%
% ⭐ IT BROWSES A RECORDING, NOT A FILE. A recording usually carries more than one pyramid:
% the tile store with bands and moments from 0.25 s up, and a preview store (rheome.ingest.preview)
% with min and max alone but reaching far finer. 'Store' picks; the default takes the richer
% one; the header names the file that is open and the recording it came from.
%
% ⭐ A TILE IS A PLACE TO PUT THINGS. Click a row to focus a channel and select the tile
% under the pointer, then app.measure("amplitude", v) writes that number onto that tile,
% keyed by (recording, scope, unit, level, k) -- the pyramid's own key. It lands in the
% second bookkeeping matrix (rheome.select.measure), which is sparse and does NOT merge, so it can
% hold anything a moment cannot: a vortex count, a rotation sense, a hand-read amplitude.
% Marks for what is stored are drawn on the rows, and rheome.select.context joins them back to the
% statistics of their tile and of any ancestor.
%
% ⚠ RAW SAMPLES ARE LOADED ON REQUEST, NEVER BY ZOOMING. The load button (key l) reads the
% current window for the focused channel and the samples stay until you ask again; nothing
% about panning or zooming touches the recording. A read pulls every channel's hyperslab, so
% it is the one expensive thing here and it is the one thing that is explicit.
%
% ⚠ THE DIAGONAL IS A FLOOR. A store keeps per-channel min and max down to
% grid.envelopeLevel and its other columns only from grid.channelLevel; below those the view
% stops refining and the readout says which constraint bit.
%
% VIEWS
%   'channels'  the stack: one min/max band per selected channel, on a shared scale, with
%               the channel list beside it. This is the explorer.
%   'detail'    one channel in full: a statistic or its envelope, the constant-Q band strip,
%               and the raw trace. This is the analysis view, and it needs a tile store.
%
% CONTROLS  the channel list, View, Store-driven statistic and band menus, Level (Auto or
%           pinned), zoom / pan / full extent, load, Snapshot.
%           Keyboard: left/right pan, +/- zoom, 0 full extent, l load, r raw panel,
%           v switch view, e envelope, s strip mode, [ ] step band, a auto.
%
% THE CALLBACKS ARE THE API. setWindow/setLevel/setUnits/setView/loadRaw are ordinary methods
% the controls call, so a script or a test drives the browser without synthesising UI events.
%
% See also: rheome.select.derive, rheome.select.open, rheome.select.frames, rheome.ingest.build, rheome.flowbrowser
%
% Author: Diellor Basha, 2026

    properties
        Stat    (1,:) char = 'rms'              % any rheome.select.derive scalar statistic
        Band    (1,1) double = 1                % band_id: the strip cursor, and the band of a band stat
        Unit    (1,1) double = 1                % channel_id, or group_id when Scope is group
        Scope   (1,:) char = 'channel'          % 'channel' | 'group'
        Budget = "auto"                         % tiles on screen: "auto" = the axis width in
                                                % pixels (one tile per pixel), or a number
        StripMode (1,:) char = 'natural'        % 'natural': each band at its own tile length
                                                % 'level':   every band at the current level
        FollowBand (1,1) logical = true         % choosing a band sets the level and the window
        Envelope (1,1) logical = false          % detail view: the min/max band, not a statistic
        View    (1,:) char = 'channels'         % 'channels': stacked min/max envelopes
                                                % 'detail':   one channel, bands and statistics
        Units   double = []                     % the channels stacked in the channels view
        MaxRows (1,1) double = 32               % how many envelopes are drawn at once
        Gain    (1,1) double = 1                % the stack's vertical gain (keys , and .)
        TilesPerView (1,1) double = 24          % how many of a band's tiles a band jump shows
        RawSeconds (1,1) double = 20            % the widest window loadRaw will read
        ShowRaw (1,1) logical = true            % the raw panel is drawn at all
    end

    properties (SetAccess = private)
        Db                                      % the select handle (open store)
        Name        (1,:) char = ''
        StoreKind   (1,:) char = 'tile'         % which pyramid is open: tile | preview | native
        Window      (1,2) double = [0 1]        % visible span, seconds
        Level       (1,1) double = 0            % the level in view
        Auto        (1,1) logical = true        % level chosen by the level-of-detail rule
        Floor       (1,1) double = 0            % the finest level the diagonal allows
        FloorReason (1,:) char = ''
        Rows        = table()                   % rheome.select.derive for the visible tiles
        Strip       = struct()                  % the strip as drawn: one entry per band row
        Ladder      = table()                   % rheome.select.ladder: tile, cycles and rate per band
        Raw         = struct('t', [], 'x', [], 'span', [NaN NaN], 'note', '')   % raw, once loaded
        Stack       = struct('units', [], 'k', [], 'tLo', [], 'tHi', [], 'min', [], 'max', [])
        Marks       = table()                   % rheome.select.measures on this store, cached
        Selected    = struct()                  % the clicked tile: level, k, parent, children
        Fig                                     % the uifigure
        StatAx; TileAx; RawAx
        StatusLbl; LevelLbl; FileLbl
    end

    properties (Access = private)
        StatDrop; BandDrop; UnitDrop; ScopeDrop; LevelSpin; AutoBox; EnvBox; LoadBtn
        ChanList; ViewDrop
        Cursors = []
        Busy (1,1) logical = false
    end

    % The raw reader is kept open across redraws: constructing a pagedrecording re-reads the
    % store header, which is wasted work when only the window moved.
    properties (Hidden)
        Pr_ = []
        PrUnit_ = NaN
        UnitIds_ = []                           % the ids behind the unit menu's strings
        UnitNames_ = {}
        LevelIds_ = []                          % the levels behind the level menu's strings
        CB_ = []                                % the strip's colour bar, kept across redraws
    end

    methods
        function app = recordingbrowser(src, varargin)
            p = inputParser;
            p.addParameter('Stat',    'rms');
            p.addParameter('Band',    []);
            p.addParameter('Unit',    []);
            p.addParameter('Scope',   'channel');
            p.addParameter('Window',  []);
            p.addParameter('Budget',  "auto");
            p.addParameter('Visible', true);
            p.addParameter('RawSeconds', 20);
            p.addParameter('Store', 'auto');
            p.addParameter('StripMode', 'natural');
            p.addParameter('FollowBand', true);
            p.addParameter('Envelope', false);
            p.addParameter('View',  'channels');
            p.addParameter('Units', []);
            p.addParameter('MaxRows', 32);
            p.addParameter('Gain', 1);
            p.parse(varargin{:});
            o = p.Results;

            [app.Db, app.StoreKind] = rheome.recordingbrowser.resolve(src, o.Store);
            app.Name   = char(app.Db.recording_id);
            app.Scope  = char(o.Scope);
            app.Budget = o.Budget;
            app.RawSeconds = o.RawSeconds;
            app.Stat   = char(o.Stat);
            app.StripMode  = char(o.StripMode);
            app.FollowBand = logical(o.FollowBand);
            app.Envelope   = logical(o.Envelope);
            app.View       = char(o.View);
            app.MaxRows    = o.MaxRows;
            app.Gain       = o.Gain;

            b = app.Db.bands;
            if isempty(b)
                % a preview store (rheome.ingest.preview): min and max, no bank, no bands. The
                % envelope is the only view it can serve, so it opens on it.
                app.Band = 0;  app.Envelope = true;  app.StripMode = 'level';
                if ~ismember(app.Stat, {'min','max','ptp'}), app.Stat = 'ptp'; end
            elseif isempty(o.Band)
                oc = find(strcmp(string(b.kind), "octave") & b.fLo >= 8, 1, 'last');
                if isempty(oc), oc = 1; end
                app.Band = b.j(oc);                       % an alpha-ish octave by default
            else
                app.Band = o.Band;
            end
            app.Unit = rheome.recordingbrowser.unitid(app.Db, app.Scope, o.Unit);
            if isempty(o.Units)
                app.Units = 1:min(app.MaxRows, app.Db.meta.C);    % the first screenful
            else
                if isstring(o.Units) || iscellstr(o.Units)
                    u = cellstr(o.Units);
                    app.Units = arrayfun(@(i) rheome.recordingbrowser.unitid(app.Db, app.Scope, u{i}), 1:numel(u));
                else
                    app.Units = arrayfun(@(u) rheome.recordingbrowser.unitid(app.Db, app.Scope, u), double(o.Units(:)'));
                end
                % ⚠ named channels win the focus. Selecting a set and being focused on a
                % channel outside it would load raw samples from a row that is not drawn.
                if isempty(o.Unit), app.Unit = app.Units(1); end
            end

            if isempty(o.Window), app.Window = [0 app.Db.meta.duration]; else, app.Window = sort(o.Window(:)'); end

            app.buildUI(o.Visible);
            app.refreshMarks();                                  % reads the sidecar once
            app.refresh();
        end

        function delete(app)
            if ~isempty(app.Fig) && isvalid(app.Fig), delete(app.Fig); end
        end
    end

    % ---- the programmatic API: every control goes through one of these ----
    methods
        setWindow(app, t0, t1)
        setLevel(app, L)
        setAuto(app, on)
        setStat(app, name)
        setBand(app, b)
        stepBand(app, d)
        setStripMode(app, m)
        setEnvelope(app, on)
        setView(app, v)
        setUnits(app, ids)
        n = loadRaw(app)
        t = measure(app, kind, value, opts)
        M = measurements(app, opts)
        refreshMarks(app)
        b = budget(app)
        setUnit(app, u)
        setScope(app, s)
        zoomBy(app, f)
        panBy(app, frac)
        fullExtent(app)
        selectTime(app, t)
        refresh(app)
        f = snapshot(app, file)
        onKey(app, e)
    end

    methods (Static)
        s = statsfor(db)
        u  = unitid(db, scope, want)
        L  = levelfor(db, width, budget, floorLevel)
        [L, why] = floorlevel(db, scope, band, stat, envelope)
        [db, choice] = resolve(src, kind)
    end

    methods (Access = private)
        buildUI(app, visible)
        drawStat(app)
        drawChannels(app)
        layout(app)
        drawTiles(app)
        drawRaw(app)
        updateLabels(app)
        syncUnits(app)
        syncLevels(app)
        setStatus(app, s)
    end
end

% Author: Diellor Basha, 2026
