function onKey(app, e)
% ONKEY  Keyboard: pan, zoom, full extent, raw on/off.
%
%   left / right    pan a quarter window (shift: a whole one)
%   up / down       coarser / finer level, pinning it
%   + / -           zoom in / out by 2          0  the whole record
%   l               LOAD the raw samples of this window (the only read of the recording)
%   r               raw panel on / off          a  hand the level back to the rule
%   [ / ]           previous / next band (with FollowBand, the window follows the band)
%   s               strip: each band at its own tile length, or all at the current level
%   e               the min/max envelope, or the statistic
%   v               switch between the channel stack and the detail view
%   , / .           the stack's gain down / up
%
% Author: Diellor Basha, 2026

    shift = iscell(e.Modifier) && any(strcmp(e.Modifier, 'shift'));
    step = 0.25;  if shift, step = 1; end
    switch e.Key
        case 'leftarrow',  app.panBy(-step);
        case 'rightarrow', app.panBy(+step);
        case 'uparrow',    app.setLevel(app.Level + 1);
        case 'downarrow',  app.setLevel(app.Level - 1);
        case {'equal','add','plus'},      app.zoomBy(0.5);
        case {'hyphen','subtract','minus'}, app.zoomBy(2);
        case '0',          app.fullExtent();
        case 'l',          app.loadRaw();
        case 'r',          app.ShowRaw = ~app.ShowRaw;  app.refresh();
        case 'a',          app.setAuto(true);
        case {'leftbracket','bracketleft'},   app.stepBand(-1);
        case {'rightbracket','bracketright'}, app.stepBand(+1);
        case 's'
            if strcmp(app.StripMode, 'natural'), app.setStripMode('level');
            else, app.setStripMode('natural'); end
        case 'e',          app.setEnvelope(~app.Envelope);
        case 'comma',  app.Gain = max(app.Gain / 1.5, 1e-3);  app.refresh();
        case 'period', app.Gain = min(app.Gain * 1.5, 1e3);   app.refresh();
        case 'v'
            if strcmp(app.View, 'channels'), app.setView('detail');
            else, app.setView('channels'); end
    end
end
% Author: Diellor Basha, 2026
