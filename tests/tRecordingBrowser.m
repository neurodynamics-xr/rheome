classdef tRecordingBrowser < matlab.unittest.TestCase
% @recordingbrowser -- the level-of-detail rule and the programmatic API.
%
% These do not test graphics. They test the two things that would make the browser lie: the
% level it chooses for a window, and whether what it draws is what rheome.select.derive returns for
% that level. The setters ARE the API, so every check drives them directly on an invisible
% figure -- which is also how the browser gets scripted and how snapshots are made headless.
%
% Author: Diellor Basha, 2026

    properties
        fx; db; App
    end

    methods (TestClassSetup)
        function open(tc)
            tc.fx = selFixture();
            tc.db = tc.fx.db;
        end
    end

    methods (TestMethodTeardown)
        function shut(tc)
            if ~isempty(tc.App) && isvalid(tc.App), delete(tc.App); end
            tc.App = [];
        end
    end

    methods
        function app = openApp(tc, varargin)
            % Most of these test the DETAIL view (one channel, bands, statistics), which was
            % the whole app before the channel stack became the default. They say so here
            % once rather than in every case; the stack has its own cases below.
            if ~any(strcmpi(varargin(1:2:end), 'View'))
                varargin = [varargin, {'View', 'detail'}];
            end
            app = rheome.recordingbrowser(tc.db, 'Visible', false, varargin{:});
            tc.App = app;
        end
    end

    methods (Test)

        function itOpensOnTheWholeRecordWithThreeAxes(tc)
            app = tc.openApp();
            tc.verifyTrue(isvalid(app.Fig));
            tc.verifyEqual(app.Window, [0 tc.db.meta.duration]);
            tc.verifyTrue(all(isvalid([app.StatAx app.TileAx app.RawAx])));
            tc.verifyEqual(app.Scope, 'channel');
            tc.verifyTrue(app.Auto);
            tc.verifyGreaterThan(height(app.Rows), 0);
        end

        function theLevelIsTheFinestUnderTheBudget(tc)
            % ⭐ THE RULE. For each width the chosen level must fit the budget, and the level
            % below it must not -- otherwise the view is coarser than it needs to be.
            g = tc.db.grid;
            for budget = [10 50 400]
                for w = [g.tExtent(end) 20 5 1 0.3]
                    L = rheome.recordingbrowser.levelfor(tc.db, w, budget, 0);
                    tc.verifyLessThanOrEqual(ceil(w / g.tExtent(L+1)), max(budget, 1), ...
                        sprintf('width %g, budget %d, level %d', w, budget, L));
                    if L > 0
                        tc.verifyGreaterThan(ceil(w / g.tExtent(L)), budget, ...
                            sprintf('level %d is not the finest for width %g', L, w));
                    end
                end
            end
        end

        function theFloorOverridesTheRule(tc)
            g = tc.db.grid;
            tc.verifyEqual(rheome.recordingbrowser.levelfor(tc.db, 0.1, 400, 4), 4);        % floor wins
            tc.verifyEqual(rheome.recordingbrowser.levelfor(tc.db, 1e9, 1, 0), g.Lmax);     % clamped at the top
        end

        function zoomingInRefinesTheLevelAndStopsAtTheFloor(tc)
            app = tc.openApp('Budget', 40, 'Stat', 'rms');
            levels = app.Level;
            for i = 1:8
                app.zoomBy(0.5);
                levels(end+1) = app.Level;                          %#ok<AGROW>
                tc.verifyLessThanOrEqual(height(app.Rows), 40 + 1); % the budget, plus the edge tile
            end
            tc.verifyTrue(all(diff(levels) <= 0), 'zooming in must never coarsen the level');
            tc.verifyEqual(app.Level, app.Floor);
            tc.verifyLessThan(levels(end), levels(1));
        end

        function zoomingOutCoarsensAndTheWidestViewIsTheWholeRecord(tc)
            app = tc.openApp('Budget', 40, 'Window', [9 10]);
            L0 = app.Level;
            for i = 1:10, app.zoomBy(2); end
            tc.verifyEqual(app.Window, [0 tc.db.meta.duration]);
            tc.verifyGreaterThan(app.Level, L0);
        end

        function whatIsDrawnIsExactlyWhatSelectDeriveReturns(tc)
            % ⭐ The browser must own no arithmetic of its own: the table it plots has to be
            % rheome.select.derive's, for the same level, unit and window.
            app = tc.openApp('Budget', 60, 'Stat', 'crest', 'Unit', 2, 'StripMode', 'level');
            app.setWindow(4, 12);
            R = rheome.select.derive(tc.db, Level=app.Level, Stats=["crest","bandPower"], ...
                              Channels=2, Window=app.Window);
            tc.verifyEqual(app.Rows.k, R.k);
            tc.verifyEqual(app.Rows.crest, R.crest);
            tc.verifyEqual(app.Rows.bandPower, R.bandPower);
            tc.verifyEqual(app.Strip.power, R.bandPower');
            tc.verifyEqual(app.Strip.bands, tc.db.grid.bandsAt{app.Level+1});
        end

        function theStatisticIsRecomputedAtTheNewLevelNotRescaled(tc)
            % ⚠ THE FAILURE THIS GUARDS. A level-of-detail plot that carried values from one
            % level to another would show the fine level's numbers on coarse tiles. rms at a
            % coarser level is the rms OF that tile -- from merged moments, not an average.
            app = tc.openApp('Budget', 1000, 'Stat', 'rms', 'Unit', 1);
            app.setLevel(2);  fine = app.Rows;
            app.setLevel(4);  coarse = app.Rows;
            tc.verifyNotEqual(height(fine), height(coarse));
            ratio = tc.db.grid.tExtent(5) / tc.db.grid.tExtent(3);
            for k = 1:height(coarse)
                i = (k-1)*ratio + (1:ratio);  i = i(i <= height(fine));
                p = sum(fine.n(i) .* fine.rms(i).^2) / sum(fine.n(i));      % power merges
                tc.verifyEqual(coarse.rms(k)^2, p, 'RelTol', 1e-10);
            end
        end

        function aPerBandStatisticRaisesTheFloorToItsOwnBand(tc)
            % The diagonal: share of a low band does not exist at a fine level, so the floor
            % moves and the reason says which band moved it.
            b = tc.db.bands;
            [~, i] = max(b.naturalLevel);
            app = tc.openApp('Stat', 'rms', 'Band', b.j(i));
            fl0 = app.Floor;
            app.setStat('share');
            tc.verifyEqual(app.Floor, b.naturalLevel(i));
            tc.verifyGreaterThan(app.Floor, fl0);
            tc.verifySubstring(app.FloorReason, 'starts at level');
            tc.verifyGreaterThanOrEqual(app.Level, app.Floor);
            app.setStat('rms');
            tc.verifyEqual(app.Floor, fl0);                      % a time stat needs no band
        end

        function aBandTheLevelDoesNotCarryIsReportedNotFaked(tc)
            % Below its natural level a band is simply not in the store. The strip must show
            % the bands the level holds and say the selected one is absent -- never borrow a
            % neighbouring band's numbers to fill the row.
            b = tc.db.bands;
            [nl, i] = max(b.naturalLevel);
            tc.assumeGreaterThan(nl, 0);
            app = tc.openApp('Stat', 'share', 'Band', b.j(i), 'Budget', 1000, ...
                             'StripMode', 'level', 'FollowBand', false);
            app.setLevel(app.Floor);                             % the finest level carrying it
            tc.verifyEqual(app.Level, nl);
            tc.verifyTrue(app.Strip.hasBand);
            tc.verifyTrue(ismember(b.j(i), app.Strip.bands));
            app.setStat('rms');                                  % a time stat frees the floor
            app.setLevel(app.Floor);
            tc.verifyLessThan(app.Level, nl);
            tc.verifyFalse(app.Strip.hasBand);
            tc.verifyFalse(ismember(b.j(i), app.Strip.bands));
            tc.verifyEqual(app.Strip.bands, tc.db.grid.bandsAt{app.Level+1});
        end

        function pinningALevelBelowTheFloorIsClampedAndSaidSo(tc)
            app = tc.openApp('Stat', 'share');
            fl = app.Floor;
            if fl == 0, app.setBand(tc.db.bands.j(end)); fl = app.Floor; end
            app.setLevel(-3);
            tc.verifyEqual(app.Level, fl);
            tc.verifyFalse(app.Auto);
            tc.verifySubstring(app.StatusLbl.String, 'below the floor');
            app.setLevel(1e6);
            tc.verifyEqual(app.Level, tc.db.grid.Lmax);
            app.setAuto(true);
            tc.verifyTrue(app.Auto);
        end

        function theWindowIsClampedRatherThanThrowing(tc)
            app = tc.openApp();
            D = tc.db.meta.duration;
            app.setWindow(-50, 1e6);
            tc.verifyEqual(app.Window, [0 D]);
            app.setWindow(12, 4);                                 % reversed
            tc.verifyEqual(app.Window, [4 12]);
            app.setWindow(D - 1, D + 100);                        % over the end
            tc.verifyLessThanOrEqual(app.Window(2), D);
            app.setWindow(5, 5);                                  % degenerate: one level-0 tile
            tc.verifyEqual(diff(app.Window), tc.db.grid.tExtent(1));
            app.panBy(-100);
            tc.verifyGreaterThanOrEqual(app.Window(1), 0);
            app.panBy(+100);
            tc.verifyLessThanOrEqual(app.Window(2), D);
            tc.verifyError(@() app.zoomBy(0), 'recordingbrowser:zoom');
        end

        function keyboardAndButtonsGoThroughTheSameSetters(tc)
            app = tc.openApp('Window', [8 12], 'Budget', 40);
            w0 = app.Window;
            app.onKey(struct('Key', 'rightarrow', 'Modifier', {{}}));
            tc.verifyEqual(app.Window, w0 + 0.25*diff(w0), 'AbsTol', 1e-9);
            app.onKey(struct('Key', 'leftarrow', 'Modifier', {{}}));
            tc.verifyEqual(app.Window, w0, 'AbsTol', 1e-9);
            app.onKey(struct('Key', 'equal', 'Modifier', {{}}));
            tc.verifyEqual(diff(app.Window), diff(w0)/2, 'AbsTol', 1e-9);
            app.onKey(struct('Key', '0', 'Modifier', {{}}));
            tc.verifyEqual(app.Window, [0 tc.db.meta.duration]);
            L = app.Level;
            app.onKey(struct('Key', 'uparrow', 'Modifier', {{}}));
            tc.verifyEqual(app.Level, min(L + 1, tc.db.grid.Lmax));
            tc.verifyFalse(app.Auto);
            app.onKey(struct('Key', 'a', 'Modifier', {{}}));
            tc.verifyTrue(app.Auto);
            tc.verifyEqual(app.Level, L);
        end

        function selectingATileReportsItsKeyParentAndChildren(tc)
            % The hierarchy is integer arithmetic on the key, which is what lets a database do
            % the same join without an index.
            app = tc.openApp('Budget', 1000, 'Stat', 'rms');
            app.setLevel(3);
            app.selectTime(9.7);
            s = app.Selected;
            g = tc.db.grid;
            tc.verifyEqual(s.level, 3);
            tc.verifyEqual(s.k, floor(9.7 / g.tExtent(4)) + 1);
            tc.verifyEqual(s.parent, [4, floor((s.k-1)/2) + 1]);
            tc.verifyEqual(s.children, [2, 2*s.k-1; 2, 2*s.k]);
            tc.verifyTrue(s.t_lo <= 9.7 && 9.7 < s.t_hi);
            tc.verifySubstring(app.StatusLbl.String, 'parent');
            app.selectTime(1e6);                                  % clamped, not thrown
            tc.verifyEqual(app.Selected.k, g.K(4));
        end

        function rawSamplesArriveOnlyWhenAskedFor(tc)
            % ⭐ PHASE 2 IS IMPERATIVE. Everything else is answered from the pyramid; a raw
            % read pulls every channel's hyperslab, so it happens when the load button (or
            % loadRaw, or key l) says so and never because a zoom got short.
            app = tc.openApp('RawSeconds', 8, 'Unit', 3);
            app.setWindow(6, 7);
            tc.verifyEmpty(app.Raw.x);                            % a short window reads nothing
            tc.verifySubstring(app.Raw.note, 'press load');
            n = app.loadRaw();
            fs = tc.db.meta.fs;
            a = floor(6*fs) + 1;  b = ceil(7*fs);
            tc.verifyEqual(n, b - a + 1);
            tc.verifyEqual(app.Raw.x, tc.fx.X(a:b, 3), 'RelTol', 1e-12);
            tc.verifyEqual(app.Raw.t, ((a:b)' - 1)/fs, 'AbsTol', 1e-12);
            tc.verifyEqual(app.Raw.span, [(a-1)/fs, b/fs], 'AbsTol', 1e-12);
            tc.verifySubstring(app.Raw.note, 'loaded');
        end

        function whatIsLoadedStaysLoadedUntilAskedAgain(tc)
            app = tc.openApp('RawSeconds', 8, 'Unit', 1);
            app.setWindow(6, 7);  app.loadRaw();
            x = app.Raw.x;
            app.panBy(2);                                         % off the loaded span
            tc.verifyEqual(app.Raw.x, x);                         % kept, not discarded
            tc.verifySubstring(app.Raw.note, '6.0-7.0 s');
            tc.verifySubstring(app.Raw.note, 'press load');
            app.setWindow(6, 7);                                  % back onto the loaded span
            app.setUnit(2);                                       % the samples stay with unit 1
            tc.verifyEqual(app.Raw.x, x);
            tc.verifyEqual(app.Raw.unit, 1);
            tc.verifySubstring(app.Raw.note, 'focus is');         % and say whose they are
            app.setWindow(6, 7);
            app.onKey(struct('Key', 'l', 'Modifier', {{}}));      % now load unit 2's
            tc.verifyEqual(app.Raw.unit, 2);
            tc.verifyEqual(app.Raw.x, tc.fx.X(round(6*tc.fx.fs)+1:round(7*tc.fx.fs), 2), 'RelTol', 1e-12);
        end

        function aWindowWiderThanRawSecondsIsRefusedWithItsCost(tc)
            % ⚠ A read moves every channel, so the refusal names the megabytes rather than
            % stalling on them.
            app = tc.openApp('RawSeconds', 2);
            app.setWindow(0, tc.db.meta.duration);
            tc.verifySubstring(app.Raw.note, 'narrow it below');
            n = app.loadRaw();
            tc.verifyEqual(n, 0);
            tc.verifyEmpty(app.Raw.x);
            tc.verifySubstring(get(app.StatusLbl, 'String'), 'refused');
            tc.verifySubstring(get(app.StatusLbl, 'String'), 'MB');
            app.ShowRaw = false;  app.refresh();
            tc.verifySubstring(app.Raw.note, 'off');
        end

        function theLoadButtonGoesThroughTheSameMethod(tc)
            app = tc.openApp('RawSeconds', 8, 'Unit', 1);
            app.setWindow(6, 7);
            c = findobj(app.Fig, 'Type', 'uicontrol');
            c = c(strcmp({c.Style}, 'pushbutton'));
            bt = c(strcmp({c.String}, 'load'));
            tc.verifyNumElements(bt, 1);
            feval(get(bt, 'Callback'), bt, []);
            tc.verifyNotEmpty(app.Raw.x);
        end

        function groupScopeReadsTheSensorTreeAndHasNoRawSignal(tc)
            app = tc.openApp('Unit', 2);
            app.setScope('group');
            tc.verifyEqual(app.Scope, 'group');
            tc.verifyTrue(ismember(app.Unit, tc.db.groupNodes));
            R = rheome.select.derive(tc.db, Level=app.Level, Stats=["rms","bandPower"], ...
                              Scope="group", Groups=app.Unit, Window=app.Window);
            tc.verifyEqual(app.Rows.rms, R.rms);
            app.setWindow(6, 7);
            tc.verifyEmpty(app.Raw.x);
            tc.verifySubstring(app.Raw.note, 'per channel');
            app.setScope('channel');
            tc.verifyEqual(app.Unit, 1);                          % a node id is not a channel id
        end

        function unitsAreNamedOrNumberedAndBadOnesAreRefused(tc)
            app = tc.openApp();
            app.setUnit('C');
            tc.verifyEqual(app.Unit, 3);
            app.setUnit(2);
            tc.verifyEqual(app.Unit, 2);
            tc.verifyError(@() app.setUnit('nosuchsensor'), 'recordingbrowser:unit');
            tc.verifyError(@() app.setUnit(99), 'recordingbrowser:unit');
            tc.verifyError(@() app.setBand(999), 'recordingbrowser:band');
            tc.verifyError(@() app.setStat('nonsense'), 'recordingbrowser:stat');
            tc.verifyError(@() app.setScope('cortex'), 'recordingbrowser:scope');
        end

        function everyEntryInTheLevelMenuChangesTheLevel(tc)
            % ⚠ THE BUG THIS PINS. The menu used to list 0..Lmax whatever the store held, so
            % on a diagonal store more than half its entries clamped back to the floor and
            % the window did not move -- a control that looks broken. It now lists the levels
            % that exist, and the test drives the control's OWN callback, not setLevel, since
            % that is the wiring that failed.
            app = tc.openApp('Stat', 'rms', 'Budget', 1000);
            lv = i_levelmenu(app);
            tc.verifyEqual(get(lv, 'Enable'), 'off');            % auto owns the level
            ab = i_checkbox(app, 'auto');
            set(ab, 'Value', 0);  feval(get(ab, 'Callback'), ab, []);
            tc.verifyEqual(get(lv, 'Enable'), 'on');
            items = get(lv, 'String');
            tc.verifyEqual(numel(items), tc.db.grid.Lmax - app.Floor + 1);
            tc.verifySubstring(items{1}, sprintf('L%d', app.Floor));
            seen = [];
            for e = 1:numel(items)
                set(lv, 'Value', e);  feval(get(lv, 'Callback'), lv, []);
                tc.verifyEqual(app.Level, app.Floor + e - 1, ...
                    sprintf('menu entry %d (%s)', e, items{e}));
                tc.verifyEqual(height(app.Rows), tc.db.grid.K(app.Level+1));
                tc.verifyEqual(get(lv, 'Value'), e);             % the menu keeps the choice
                seen(end+1) = height(app.Rows);                  %#ok<AGROW>
            end
            tc.verifyEqual(numel(unique(seen)), numel(seen));    % every entry moved the view
        end

        function theLevelMenuFollowsTheFloorWhenTheStatisticChanges(tc)
            b = tc.db.bands;
            [nl, i] = max(b.naturalLevel);
            tc.assumeGreaterThan(nl, 0);
            app = tc.openApp('Stat', 'rms', 'Band', b.j(i));
            lv = i_levelmenu(app);
            n0 = numel(get(lv, 'String'));
            app.setStat('share');                                % per band: the floor rises
            tc.verifyEqual(app.Floor, nl);
            items = get(lv, 'String');
            tc.verifyEqual(numel(items), tc.db.grid.Lmax - nl + 1);
            tc.verifyLessThan(numel(items), n0);
            tc.verifySubstring(items{1}, sprintf('L%d', nl));
            tc.verifyEqual(tc.db.grid.Lmax, app.Floor + numel(items) - 1);
            app.setStat('rms');
            tc.verifyEqual(numel(get(lv, 'String')), n0);        % and back down again
        end

        function choosingABandSetsTheLevelAndTheTimeWindow(tc)
            % ⭐ THE LADDER READ BACKWARDS. A band fixes its support, the support fixes the
            % shortest tile that can hold it, and the tile fixes the level -- so asking for a
            % band is asking for a time window. Stepping down an octave must double the tile
            % and double the window, leaving the number of tiles on screen unchanged.
            app = tc.openApp('Scope', 'group', 'Stat', 'rms');      % no channel floor
            b = tc.db.bands;
            widths = [];  levels = [];
            for j = b.j(:)'
                app.setBand(j);
                k = find(b.j == j, 1);
                tc.verifyEqual(app.Level, min(b.naturalLevel(k), tc.db.grid.Lmax), ...
                    sprintf('band %d', j));
                w = app.TilesPerView * tc.db.grid.tExtent(app.Level+1);
                tc.verifyEqual(diff(app.Window), min(w, tc.db.meta.duration), 'AbsTol', 1e-9);
                widths(end+1) = diff(app.Window);  levels(end+1) = app.Level;   %#ok<AGROW>
            end
            tc.verifyTrue(all(diff(levels) >= 0));                  % the ladder is monotone
            tc.verifyTrue(all(diff(widths) >= 0));
        end

        function theChannelFloorStillOverridesABandChoice(tc)
            % In channel scope a band below the diagonal cannot get its own tile; the level
            % clamps and the readout says the sensor's tile instead.
            db2 = tc.db;
            fl = db2.channelLevel;
            b = tc.db.bands;
            j = b.j(find(b.naturalLevel < max(fl, 1), 1));
            app = tc.openApp('Scope', 'channel', 'Stat', 'rms');
            app.setBand(j);
            k = find(b.j == j, 1);
            tc.verifyEqual(app.Level, max(b.naturalLevel(k), fl));
            tc.verifyGreaterThanOrEqual(app.Level, b.naturalLevel(k));
        end

        function theNaturalStripDrawsEveryBandAtItsOwnTileLength(tc)
            % Nothing is shaded in this mode: every band is drawn at a level that carries it,
            % and each row's level is the finest one at or above the band's natural level
            % that fits the budget.
            app = tc.openApp('Scope', 'group', 'StripMode', 'natural', 'Budget', 64);
            app.setWindow(0, tc.db.meta.duration);
            S = app.Strip;
            b = tc.db.bands;
            tc.verifyEqual(S.mode, 'natural');
            tc.verifyEqual(sort(S.bands), sort(b.j(:)'));            % every band has a row
            tc.verifyTrue(S.hasBand);
            for i = 1:numel(S.rows)
                r = S.rows(i);
                k = find(b.j == r.band, 1);
                tc.verifyGreaterThanOrEqual(r.level, b.naturalLevel(k));
                tc.verifyEqual(r.level, rheome.recordingbrowser.levelfor(tc.db, diff(app.Window), 64, b.naturalLevel(k)));
                tc.verifyEqual(numel(r.power), numel(r.tLo));
                tc.verifyLessThanOrEqual(numel(r.power), 64 + 1);
                R = rheome.select.derive(tc.db, Level=r.level, Scope="group", Groups=app.Unit, ...
                                  Stats="bandPower", Bands=r.band, Window=app.Window);
                tc.verifyEqual(r.power(:), R.bandPower(:));          % the strip owns no arithmetic
            end
            tc.verifyGreaterThan(numel(unique(S.levels)), 1);        % it really is a staircase
        end

        function theTwoStripModesDisagreeOnlyWhereTheDiagonalDoes(tc)
            app = tc.openApp('Scope', 'group', 'StripMode', 'level', 'Budget', 400);
            app.setWindow(0, tc.db.meta.duration);
            lv = app.Strip;
            tc.verifyEqual(lv.bands, tc.db.grid.bandsAt{app.Level+1});
            tc.verifyTrue(all(lv.levels == app.Level));
            app.setStripMode('natural');
            nat = app.Strip;
            tc.verifyEqual(sort(nat.bands), sort(tc.db.bands.j(:)'));
            tc.verifyGreaterThanOrEqual(numel(nat.bands), numel(lv.bands));
            tc.verifyError(@() app.setStripMode('spectrogram'), 'recordingbrowser:stripmode');
            app.onKey(struct('Key', 's', 'Modifier', {{}}));
            tc.verifyEqual(app.StripMode, 'level');
        end

        function steppingBandsWalksTheLadder(tc)
            app = tc.openApp('Scope', 'group');
            js = tc.db.bands.j(:)';
            app.setBand(js(1));
            app.onKey(struct('Key', 'rightbracket', 'Modifier', {{}}));
            tc.verifyEqual(app.Band, js(2));
            app.onKey(struct('Key', 'leftbracket', 'Modifier', {{}}));
            tc.verifyEqual(app.Band, js(1));
            app.stepBand(-5);
            tc.verifyEqual(app.Band, js(1));                        % clamped at the ends
            app.setBand(js(end));  app.stepBand(+5);
            tc.verifyEqual(app.Band, js(end));
        end

        function theReadoutCarriesTheBandsRowOfTheLadder(tc)
            app = tc.openApp('Scope', 'group', 'Stat', 'rms');
            app.setBand(tc.db.bands.j(2));
            s = get(app.LevelLbl, 'String');
            tc.verifyClass(s, 'cell');
            tc.verifyNumElements(s, 2);
            tc.verifySubstring(s{2}, 'cycles');
            tc.verifySubstring(s{2}, 'support');
            tc.verifySubstring(s{2}, sprintf('band %d', tc.db.bands.j(2)));
            tc.verifyEqual(height(app.Ladder), height(tc.db.bands));
        end

        function theEnvelopeBandContainsEveryRawSampleItCovers(tc)
            % ⭐ THE PROOF. min and max merge by extremum and are stored rounded outward, so
            % the band drawn at ANY level contains every sample of the span it covers. A
            % viewer that decimated by sampling could not promise this.
            app = tc.openApp('Envelope', true, 'Unit', 1, 'Budget', 1000);
            tc.verifyTrue(app.Envelope);
            tc.verifyTrue(all(ismember({'min','max'}, app.Rows.Properties.VariableNames)));
            for L = [0 2 4]
                app.setLevel(L);
                R = app.Rows;
                for k = 1:height(R)
                    a = round(R.t_lo(k) * tc.fx.fs) + 1;
                    b = min(round(R.t_hi(k) * tc.fx.fs), tc.fx.nT);
                    x = tc.fx.X(a:b, 1);
                    tc.verifyLessThanOrEqual(R.min(k), min(x), sprintf('level %d tile %d', L, k));
                    tc.verifyGreaterThanOrEqual(R.max(k), max(x), sprintf('level %d tile %d', L, k));
                end
            end
        end

        function theRawTraceStaysInsideTheDrawnBand(tc)
            % The same statement where a person can see it: a short window loads raw samples,
            % and every one of them must sit inside the tile it belongs to.
            app = tc.openApp('Envelope', true, 'Unit', 3, 'RawSeconds', 8);
            app.setWindow(5, 11);
            app.loadRaw();                                        % phase 2 is asked for
            tc.verifyNotEmpty(app.Raw.x);
            R = app.Rows;
            for k = 1:height(R)
                in = app.Raw.t >= R.t_lo(k) & app.Raw.t < R.t_hi(k);
                if ~any(in), continue; end
                tc.verifyGreaterThanOrEqual(min(app.Raw.x(in)), R.min(k));
                tc.verifyLessThanOrEqual(max(app.Raw.x(in)), R.max(k));
            end
        end

        function aSingleTileInViewStillDraws(tc)
            % ⚠ A window inside one tile is the normal end of a zoom, and repelem of a scalar
            % returns a row: the fill threw inside a callback until the shapes were forced.
            app = tc.openApp('Envelope', true, 'Budget', 1000);
            app.setLevel(tc.db.grid.Lmax);
            app.setWindow(5, 6);
            tc.verifyEqual(height(app.Rows), 1);
            tc.verifyTrue(isvalid(app.StatAx));
            tc.verifyNotEmpty(findobj(app.StatAx, 'Type', 'patch'));
            app.setEnvelope(false);
            app.setEnvelope(true);
            tc.verifyEqual(height(app.Rows), 1);
        end

        function theEnvelopeHasItsOwnLowerFloor(tc)
            % On a store with a channel floor, a statistic stops there and the envelope keeps
            % going, because min and max are kept below it.
            d = tc.db;  d.channelLevel = 6;  d.envelopeLevel = 0;
            [Ls, whyS] = rheome.recordingbrowser.floorlevel(d, 'channel', d.bands.j(1), 'rms', false);
            [Le, whyE] = rheome.recordingbrowser.floorlevel(d, 'channel', d.bands.j(1), 'rms', true);
            tc.verifyEqual(Ls, 6);
            tc.verifyEqual(Le, 0);
            tc.verifySubstring(whyS, 'rows start at level 6');
            tc.verifySubstring(whyE, 'the record');
            % group scope has no floor either way
            tc.verifyEqual(rheome.recordingbrowser.floorlevel(d, 'group', d.bands.j(1), 'rms', true), 0);
        end

        function theEnvelopeTogglesFromTheKeyboardAndTheCheckbox(tc)
            app = tc.openApp('Stat', 'rms');
            tc.verifyFalse(app.Envelope);
            app.onKey(struct('Key', 'e', 'Modifier', {{}}));
            tc.verifyTrue(app.Envelope);
            tc.verifyFalse(ismember('rms', app.Rows.Properties.VariableNames));
            app.onKey(struct('Key', 'e', 'Modifier', {{}}));
            tc.verifyFalse(app.Envelope);
            tc.verifyTrue(ismember('rms', app.Rows.Properties.VariableNames));
            env = i_checkbox(app, 'env');
            set(env, 'Value', 1);  feval(get(env, 'Callback'), env, []);
            tc.verifyTrue(app.Envelope);
        end

        function theStripKeepsALevelThatCarriesBandsWhileTheEnvelopeGoesFiner(tc)
            % In envelope mode the axis can sit below the channel floor, where no band exists.
            % The strip must not follow it there, and must say which level it used.
            app = tc.openApp('Envelope', true, 'StripMode', 'level', 'Budget', 1000);
            app.setLevel(0);
            tc.verifyEqual(app.Level, 0);
            tc.verifyGreaterThanOrEqual(app.Strip.level, min(tc.db.bands.naturalLevel));
            tc.verifyNotEmpty(app.Strip.bands);
        end

        function theBudgetIsTheAxisWidthInPixelsUnlessItIsGiven(tc)
            % ⭐ A tile narrower than a pixel is invisible work; a tile several pixels wide is
            % a visible box where the store has detail. The default asks the axis how wide it
            % is, so the view is right on any screen and follows a resize.
            app = tc.openApp();
            p = getpixelposition(app.StatAx);
            tc.verifyEqual(app.budget(), round(p(3)));
            tc.verifyGreaterThan(app.budget(), 100);
            set(app.Fig, 'Position', [60 60 2560 1000]);
            drawnow;
            p2 = getpixelposition(app.StatAx);
            tc.verifyGreaterThan(p2(3), p(3));
            tc.verifyEqual(app.budget(), round(p2(3)));
            app.Budget = 250;                                    % an explicit number wins
            tc.verifyEqual(app.budget(), 250);
            app.Budget = "auto";
            tc.verifyEqual(app.budget(), round(p2(3)));
        end

        function theAutoBudgetKeepsTilesAtAboutOnePixel(tc)
            % Dyadic levels mean the count lands between half the budget and the budget, so a
            % tile is between one and two pixels wide -- never the 4 px boxes a fixed 400 gave
            % on a 1100 px axis. ⚠ UNLESS THE STORE RUNS OUT FIRST: on a 20 s fixture the
            % finest level is 80 tiles, so the pixels per tile are bounded by the data and the
            % level sits on its floor. Both cases are the rule working.
            app = tc.openApp();
            app.setWindow(0, tc.db.meta.duration);
            px = getpixelposition(app.StatAx);
            perTile = px(3) / height(app.Rows);
            tc.verifyTrue(perTile <= 2.5 || app.Level == app.Floor, ...
                sprintf('%.1f px per tile at level %d (floor %d)', perTile, app.Level, app.Floor));
            tc.verifyGreaterThanOrEqual(height(app.Rows), 1);
            app.Budget = 40;                                     % deliberately coarse
            app.refresh();
            tc.verifyGreaterThan(px(3) / height(app.Rows), perTile);
        end

        function aPreviewStoreOpensAsAnExplorerWithNoBands(tc)
            % ⭐ The explorer case: a recording with only a min/max pyramid beside it. The
            % browser must open on it, browse channels, zoom, and say what it cannot do --
            % not fail on a missing array.
            f = rheome.ingest.preview(tc.fx.rec, Floor=0.05, Verbose=false, Overwrite=true);
            db2 = rheome.select.open(f);
            app = rheome.recordingbrowser(db2, 'Visible', false, 'View', 'detail');
            tc.addTeardown(@() delete(app));
            tc.verifyTrue(app.Envelope);                       % the only view it can serve
            tc.verifyEqual(app.Band, 0);
            tc.verifyTrue(ismember(app.Stat, {'min','max','ptp'}));
            tc.verifyEmpty(app.Strip.bands);
            av = rheome.recordingbrowser.statsfor(db2);
            tc.verifyEqual(sort(av(:)), ["max"; "min"; "ptp"]);
            tc.verifyError(@() app.setStat('rms'), 'recordingbrowser:stat');
            tc.verifyError(@() app.setBand(1), 'recordingbrowser:band');
            app.setUnit(2);
            tc.verifyEqual(app.Unit, 2);
            app.setWindow(4, 6);
            tc.verifyGreaterThan(height(app.Rows), 0);
            tc.verifyTrue(all(ismember({'min','max'}, app.Rows.Properties.VariableNames)));
            s = get(app.LevelLbl, 'String');
            tc.verifySubstring(s{2}, 'supports');              % the window's band support
            tc.verifySubstring(s{2}, 'preview store');
        end

        function thePreviewGoesFinerThanTheTileStoreDoes(tc)
            % The point of a finer floor: the explorer keeps pixel-exact detail down to where
            % raw reads take over, which the 0.25 s tile store cannot.
            f = rheome.ingest.preview(tc.fx.rec, Floor=0.01, Verbose=false, Overwrite=true);
            db2 = rheome.select.open(f);
            tc.verifyLessThan(db2.grid.tExtent(1), tc.db.grid.tExtent(1));
            app = rheome.recordingbrowser(db2, 'Visible', false, 'Budget', 400);
            tc.addTeardown(@() delete(app));
            app.setWindow(4, 6);
            tc.verifyEqual(app.Level, 0);
            tc.verifyEqual(tc.db.grid.tExtent(1) / db2.grid.tExtent(1), 25, 'RelTol', 1e-9);
            R = rheome.select.derive(db2, Level=0, Channels=1, Stats=["min","max"], Window=[4 6]);
            for k = 1:height(R)
                a = round(R.t_lo(k) * tc.fx.fs) + 1;
                b = min(round(R.t_hi(k) * tc.fx.fs), tc.fx.nT);
                x = tc.fx.X(a:b, 1);
                tc.verifyLessThanOrEqual(R.min(k), min(x));
                tc.verifyGreaterThanOrEqual(R.max(k), max(x));
            end
        end

        function itOpensOnTheChannelStack(tc)
            % ⭐ THE DEFAULT VIEW IS THE POINT OF THE APP: a stored matrix, drawn as one
            % min/max band per channel, at whatever level the zoom asks for.
            app = rheome.recordingbrowser(tc.db, 'Visible', false);
            tc.addTeardown(@() delete(app));
            tc.verifyEqual(app.View, 'channels');
            tc.verifyEqual(app.Units, 1:min(app.MaxRows, tc.db.meta.C));
            tc.verifyEqual(app.Unit, 1);
            tc.verifyEqual(numel(app.Stack.units), tc.db.meta.C);
            tc.verifySize(app.Stack.min, [numel(app.Stack.k), tc.db.meta.C]);
            tc.verifyEqual(app.Stack.max >= app.Stack.min, true(size(app.Stack.min)));
        end

        function theStackIsExactlyWhatDeriveGivesForThoseChannels(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', [1 3]);
            tc.addTeardown(@() delete(app));
            app.setWindow(4, 12);
            R = rheome.select.derive(tc.db, Level=app.Level, Stats=["min","max"], ...
                              Channels=[1 3], Window=app.Window);
            for i = 1:2
                r = R(R.unit_id == app.Stack.units(i), :);
                tc.verifyEqual(app.Stack.min(:, i), r.min);
                tc.verifyEqual(app.Stack.max(:, i), r.max);
            end
            tc.verifyEqual(app.Stack.tLo, R.t_lo(1:numel(app.Stack.k)));
        end

        function everyRowOfTheStackContainsItsSamples(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', 1:tc.fx.C);
            tc.addTeardown(@() delete(app));
            for L = [0 3]
                app.setLevel(L);
                S = app.Stack;
                for i = 1:numel(S.units)
                    for k = 1:numel(S.k)
                        a = round(S.tLo(k) * tc.fx.fs) + 1;
                        b = min(round(S.tHi(k) * tc.fx.fs), tc.fx.nT);
                        x = tc.fx.X(a:b, S.units(i));
                        tc.verifyLessThanOrEqual(S.min(k, i), min(x));
                        tc.verifyGreaterThanOrEqual(S.max(k, i), max(x));
                    end
                end
            end
        end

        function theChannelListDrivesTheStack(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false);
            tc.addTeardown(@() delete(app));
            c = findobj(app.Fig, 'Type', 'uicontrol');
            lb = c(strcmp({c.Style}, 'listbox'));
            tc.verifyNumElements(lb, 1);
            tc.verifyEqual(numel(get(lb, 'String')), tc.db.meta.C);
            set(lb, 'Value', [2 4]);
            feval(get(lb, 'Callback'), lb, []);
            tc.verifyEqual(app.Units, [2 4]);
            tc.verifyEqual(app.Unit, 2);                       % the first is the focus
            tc.verifyEqual(app.Stack.units, [2 4]);
            app.setUnits("C");
            tc.verifyEqual(app.Units, 3);
            tc.verifyError(@() app.setUnits([]), 'recordingbrowser:units');
        end

        function onlyMaxRowsAreDrawnAndTheTitleSaysSo(tc)
            % ⚠ Each row is a filled band of a few thousand vertices; a few dozen draw
            % instantly and a few hundred do not. The selection keeps them all.
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', 1:tc.fx.C, 'MaxRows', 2);
            tc.addTeardown(@() delete(app));
            tc.verifyEqual(numel(app.Units), tc.fx.C);
            tc.verifyEqual(numel(app.Stack.units), 2);
            t = get(get(app.StatAx, 'Title'), 'String');
            tc.verifySubstring(t, sprintf('of %d selected', tc.fx.C));
        end

        function theHeaderNamesTheStoreAndTheRecording(tc)
            % A picture of a stored matrix that does not say which matrix is a screenshot
            % waiting to be misattributed.
            app = rheome.recordingbrowser(tc.db, 'Visible', false);
            tc.addTeardown(@() delete(app));
            h = get(app.FileLbl, 'String');
            [~, stem, ext] = fileparts(tc.db.file);
            tc.verifySubstring(h, [stem ext]);
            tc.verifySubstring(h, 'tile store');
            tc.verifySubstring(h, 'recording.mat');
            tc.verifySubstring(h, sprintf('%d ch', tc.db.meta.C));
            tc.verifySubstring(h, sprintf('%g Hz', tc.db.meta.fs));
        end

        function switchingViewsKeepsTheWindowAndTheFocus(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', [2 3]);
            tc.addTeardown(@() delete(app));
            app.setWindow(5, 9);
            w = app.Window;
            app.onKey(struct('Key', 'v', 'Modifier', {{}}));
            tc.verifyEqual(app.View, 'detail');
            tc.verifyEqual(app.Window, w);
            tc.verifyEqual(app.Unit, 2);
            tc.verifyGreaterThan(height(app.Rows), 0);
            app.setView('channels');
            tc.verifyEqual(app.View, 'channels');
            tc.verifyEqual(app.Stack.units, [2 3]);
            tc.verifyError(@() app.setView('waterfall'), 'recordingbrowser:view');
        end

        function theGainScalesTheStackWithoutTouchingTheData(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', 1:tc.fx.C);
            tc.addTeardown(@() delete(app));
            before = app.Stack.max;
            app.onKey(struct('Key', 'period', 'Modifier', {{}}));
            tc.verifyGreaterThan(app.Gain, 1);
            tc.verifyEqual(app.Stack.max, before);             % a drawing gain, not a rescale
            app.onKey(struct('Key', 'comma', 'Modifier', {{}}));
            tc.verifyEqual(app.Gain, 1, 'RelTol', 1e-12);
        end

        function clickingATileGivesTheKeyToWriteAMeasurementOn(tc)
            % ⭐ THE LOOP THE APP EXISTS FOR: look, pick a tile, put a number on it. The key
            % is the pyramid's own -- (scope, unit, level, k) -- so the measurement lands on
            % that channel's tile and joins back to its statistics.
            if exist(tc.db.labelFile, 'file') == 2, delete(tc.db.labelFile); end
            tc.addTeardown(@() i_rm(tc.db.labelFile));
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', [1 2]);
            tc.addTeardown(@() delete(app));
            app.setLevel(2);
            app.selectTime(9.7);
            s = app.Selected;
            tc.verifyEqual(s.unit, 1);
            tc.verifyEqual(s.scope, 'channel');
            tc.verifyEqual(s.k, floor(9.7 / tc.db.grid.tExtent(3)) + 1);
            t = app.measure("amplitude", 4.2e-13);
            tc.verifyFalse(t.existed);
            M = app.measurements();
            tc.verifyEqual(height(M), 1);
            tc.verifyEqual(M.unit_id, s.unit);
            tc.verifyEqual(M.level, s.level);
            tc.verifyEqual(M.k, s.k);
            tc.verifyEqual(M.kind, "amplitude");
            tc.verifyEqual(M.value, 4.2e-13);
            tc.verifyEqual(M.scope, "channel");
            t2 = app.measure("amplitude", 4.2e-13);
            tc.verifyTrue(t2.existed);                            % idempotent, like any write
            tc.verifyEqual(height(app.measurements()), 1);
        end

        function aMeasurementJoinsBackToTheTilesItWasMadeOn(tc)
            if exist(tc.db.labelFile, 'file') == 2, delete(tc.db.labelFile); end
            tc.addTeardown(@() i_rm(tc.db.labelFile));
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', [1 2]);
            tc.addTeardown(@() delete(app));
            app.setLevel(2);  app.selectTime(9.7);  app.measure("amplitude", 1);
            app.setUnit(2);   app.selectTime(3.2);  app.measure("amplitude", 2);
            V = app.measurements(Kind="amplitude", Context=["min","max"], Ancestors=[0 2]);
            tc.verifyEqual(height(V), 2);
            for i = 1:2
                R = rheome.select.derive(tc.db, Level=V.level(i), Channels=V.unit_id(i), Stats=["min","max"]);
                j = find(R.k == V.k(i), 1);
                tc.verifyEqual(V.min(i), R.min(j));
                tc.verifyEqual(V.max(i), R.max(j));
            end
            tc.verifyEqual(V.k_up2, floor((V.k - 1) / 4) + 1);
            tc.verifyTrue(all(ismember(["min_up2","max_up2"], string(V.Properties.VariableNames))));
        end

        function measurementsShowAsMarksOnTheirOwnRow(tc)
            if exist(tc.db.labelFile, 'file') == 2, delete(tc.db.labelFile); end
            tc.addTeardown(@() i_rm(tc.db.labelFile));
            app = rheome.recordingbrowser(tc.db, 'Visible', false, 'Units', [1 2 3]);
            tc.addTeardown(@() delete(app));
            app.setLevel(2);  app.selectTime(9.7);  app.measure("amplitude", 1);
            tc.verifyEqual(height(app.Marks), 1);
            t = get(get(app.StatAx, 'Title'), 'String');
            tc.verifySubstring(t, '1 measured');
            tc.verifyEqual(height(app.measurements(InView=true)), 1);
            app.setWindow(0, 2);                                  % away from the marked tile
            tc.verifyEqual(height(app.measurements(InView=true)), 0);
            tc.verifyEqual(height(app.Marks), 1);                 % still cached, just not in view
            rheome.select.rollback(tc.db, app.Marks.txn_id(1));
            app.refreshMarks();
            tc.verifyEmpty(app.Marks);
        end

        function measuringWithNothingSelectedIsRefused(tc)
            app = rheome.recordingbrowser(tc.db, 'Visible', false);
            tc.addTeardown(@() delete(app));
            tc.verifyError(@() app.measure("x", 1), 'recordingbrowser:measure');
        end

        function aSnapshotOfAnInvisibleWindowStillWritesAPng(tc)
            % ⚠ THE REASON THE WINDOW IS A CLASSIC FIGURE. exportgraphics refuses an
            % unrendered uifigure, which is every headless case there is.
            app = tc.openApp('Stat', 'spectralCentroid');
            d = tempname;  mkdir(d);  tc.addTeardown(@() rmdir(d, 's'));
            f = app.snapshot(fullfile(d, 'view.png'));
            i = dir(f);
            tc.verifyEqual(exist(f, 'file'), 2);
            tc.verifyGreaterThan(i.bytes, 1e4);
            tc.verifySubstring(app.StatusLbl.String, 'saved');
        end

        function aDatasetNameOrAPathOrAHandleAllResolve(tc)
            [db1, k1] = rheome.recordingbrowser.resolve(tc.db);
            tc.verifyEqual(db1.file, tc.db.file);
            tc.verifyEqual(k1, 'tile');
            [db2, k2] = rheome.recordingbrowser.resolve(tc.fx.store);
            tc.verifyEqual(db2.file, tc.fx.store);
            tc.verifyEqual(k2, 'tile');
            tc.verifyError(@() rheome.recordingbrowser.resolve('no_such_dataset'), 'recordingbrowser:resolve');
            tc.verifyError(@() rheome.recordingbrowser.resolve(tc.db, 'sideways'), 'recordingbrowser:resolve');
        end

        function itBrowsesARecordingAndPicksBetweenItsPyramids(tc)
            % ⭐ The app is about a RECORDING, not a file: a recording can carry a tile store
            % and a preview store, and the caller says which rather than naming a path.
            [dbp, kp] = rheome.recordingbrowser.resolve(rheome.select.open(rheome.ingest.preview(tc.fx.rec, ...
                             Floor=0.05, Verbose=false, Overwrite=true)));
            tc.verifyEqual(kp, 'preview');
            tc.verifyTrue(dbp.preview);
            app = rheome.recordingbrowser(dbp, 'Visible', false);
            tc.addTeardown(@() delete(app));
            tc.verifyEqual(app.StoreKind, 'preview');
            s1 = get(app.LevelLbl, 'String');
            tc.verifySubstring(s1{1}, 'preview store');
            app2 = rheome.recordingbrowser(tc.db, 'Visible', false);
            tc.addTeardown(@() delete(app2));
            tc.verifyEqual(app2.StoreKind, 'tile');
            s2 = get(app2.LevelLbl, 'String');
            tc.verifySubstring(s2{1}, 'tile store');
        end

    end
end

% The controls are found by style, the way a person finds them by eye: these tests drive the
% callbacks the UI installs, which is the layer a test calling methods directly cannot reach.
function h = i_checkbox(app, name)
    c = findobj(app.Fig, 'Type', 'uicontrol');
    c = c(strcmp({c.Style}, 'checkbox'));
    h = c(strcmp({c.String}, name));      % by label: there is more than one checkbox now
    h = h(1);
end

function i_rm(f)
    if exist(f, 'file') == 2, delete(f); end
end

function h = i_levelmenu(app)
    c = findobj(app.Fig, 'Type', 'uicontrol');
    c = c(strcmp({c.Style}, 'popupmenu'));
    k = cellfun(@(s) iscell(s) && ~isempty(s) && startsWith(s{1}, 'L'), {c.String});
    h = c(find(k, 1));
end

% Author: Diellor Basha, 2026
