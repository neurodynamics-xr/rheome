function buildUI(app)
% BUILDUI  The window: a control strip and a tab group.
%
% Keyboard is bound at the FIGURE, so it works whichever tab has focus:
%   left/right       step one frame        shift+left/right   step ten
%   up/down          previous/next page    space              play/pause
%   m                toggle signed/magnitude
%
% Author: Diellor Basha, 2026

    app.Fig = uifigure('Name', sprintf('cortical flow — %s', app.Bundle.name), ...
                       'Position', [60 60 1420 900], ...
                       'KeyPressFcn', @(~,e) app.onKey(e));
    g = uigridlayout(app.Fig, [2 1]);
    g.RowHeight = {46, '1x'};  g.ColumnWidth = {'1x'};

    % ---- control strip ----
    app.Ctl = uigridlayout(g, [1 14]);
    app.Ctl.Layout.Row = 1;
    app.Ctl.ColumnWidth = {42, 100, 40, 80, 34, 34, 34, 46, '1x', 100, 44, 230};
    app.Ctl.Padding = [8 4 8 4];  app.Ctl.ColumnSpacing = 5;

    uilabel(app.Ctl, 'Text', 'Band');
    app.BandDrop = uidropdown(app.Ctl, 'Items', fieldnames(app.Bands), ...
        'Value', app.BandName, 'ValueChangedFcn', @(s,~) app.onBand(s.Value));

    uilabel(app.Ctl, 'Text', 'Page');
    app.PageSpin = uispinner(app.Ctl, 'Limits', [1 app.Bundle.pager.NumPages], ...
        'Value', app.PageIndex, 'Step', 1, 'RoundFractionalValues', 'on', ...
        'ValueChangedFcn', @(s,~) app.onPage(s.Value));

    % frame stepping -- a slider cannot land on a single frame reliably
    uibutton(app.Ctl, 'Text', '◀', 'Tooltip', 'previous frame  (←)', ...
        'ButtonPushedFcn', @(~,~) app.stepFrame(-1));
    uibutton(app.Ctl, 'Text', '▶', 'Tooltip', 'next frame  (→)', ...
        'ButtonPushedFcn', @(~,~) app.stepFrame(+1));
    app.PlayBtn = uibutton(app.Ctl, 'state', 'Text', '▶▶', 'Tooltip', 'play  (space)', ...
        'ValueChangedFcn', @(s,~) app.onPlay(s.Value));

    app.FrameLbl = uilabel(app.Ctl, 'Text', '', 'HorizontalAlignment', 'right');
    app.FrameSlider = uislider(app.Ctl, 'Limits', [1 2], 'Value', 1, ...
        'MajorTicks', [], 'MinorTicks', [], ...
        'ValueChangingFcn', @(~,e) app.setFrame(round(e.Value)));

    app.ModeSwitch = uidropdown(app.Ctl, 'Items', {'signed','magnitude'}, ...
        'Value', app.MapMode, 'ValueChangedFcn', @(s,~) app.onMode(s.Value));
    uilabel(app.Ctl, 'Text', '(m)');
    app.StatusLbl = uilabel(app.Ctl, 'Text', '', 'FontColor', [0.25 0.25 0.25]);

    % ---- tabs ----
    app.Tabs = uitabgroup(g, 'SelectionChangedFcn', @(~,~) app.refresh());
    app.Tabs.Layout.Row = 2;

    % Cortex: scale surfaces on top, a few sensor traces underneath
    tCortex = uitab(app.Tabs, 'Title', 'Cortex');
    gc = uigridlayout(tCortex, [2 1]);
    gc.RowHeight = {'3.4x', '1x'};  gc.RowSpacing = 2;  gc.Padding = [2 2 2 2];
    app.CortexAx  = uigridlayout(gc, [1 1]);
    app.SensorAx  = uiaxes(gc);

    tSpec = uitab(app.Tabs, 'Title', 'Spectrum');
    gs = uigridlayout(tSpec, [1 2]);
    app.SpecAxMode  = uiaxes(gs);
    app.SpecAxScale = uiaxes(gs);

    tTime = uitab(app.Tabs, 'Title', 'Timeline');
    gt = uigridlayout(tTime, [3 1]);
    gt.RowHeight = {'1.6x', '1x', '0.8x'};
    app.ScalogramAx = uiaxes(gt);
    app.TimeAxPower = uiaxes(gt);
    app.TimeAxScale = uiaxes(gt);
end

% Author: Diellor Basha, 2026
