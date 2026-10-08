function layout(app)
% LAYOUT  Place the axes and the channel list for the current view.
%
% The two views want different furniture, so rather than building two figures the same
% objects are moved and hidden: the channels view gives almost everything to the stack and
% shows the channel list; the detail view hides the list and restores the three panels.
%
% Author: Diellor Basha, 2026

    if isempty(app.Fig) || ~isvalid(app.Fig), return; end
    ch = strcmp(app.View, 'channels');
    i_vis(app.ChanList, ch);
    i_vis(app.TileAx, ~ch);
    if ch
        set(app.StatAx, 'Position', [0.115 0.215 0.860 0.660]);
        set(app.RawAx,  'Position', [0.115 0.065 0.860 0.110]);
        set(app.ChanList, 'Position', [0.008 0.065 0.095 0.810]);
    else
        set(app.StatAx, 'Position', [0.068 0.630 0.870 0.245]);
        set(app.TileAx, 'Position', [0.068 0.330 0.870 0.265]);
        set(app.RawAx,  'Position', [0.068 0.065 0.870 0.205]);
    end
    i_vis(app.BandDrop, ~ch);
    i_vis(app.StatDrop, ~ch);
    i_vis(app.UnitDrop, ~ch);
    i_vis(app.EnvBox,   ~ch);
end

function i_vis(h, on)
    if isempty(h) || ~isvalid(h), return; end
    if on, set(h, 'Visible', 'on'); else, set(h, 'Visible', 'off'); end
    if isa(h, 'matlab.graphics.axis.Axes') && ~on
        cla(h);  set(h, 'XTick', [], 'YTick', []);
    end
end
% Author: Diellor Basha, 2026
