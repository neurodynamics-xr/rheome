function onKey(app, e)
% ONKEY  Keyboard shortcuts, bound at the figure so they work on any tab.
%
%   left/right  +/- 1 frame     shift+left/right  +/- 10
%   up/down     +/- 1 page      space             play/pause     m  toggle map mode
%
% Author: Diellor Basha, 2026

    big = any(strcmp(e.Modifier, 'shift'));
    switch e.Key
        case 'rightarrow', app.stepFrame(1 + 9*big);
        case 'leftarrow',  app.stepFrame(-(1 + 9*big));
        case 'uparrow'
            if app.PageIndex > 1
                app.PageSpin.Value = app.PageIndex - 1;  app.onPage(app.PageIndex - 1);
            end
        case 'downarrow'
            if app.PageIndex < app.Bundle.pager.NumPages
                app.PageSpin.Value = app.PageIndex + 1;  app.onPage(app.PageIndex + 1);
            end
        case 'space'
            app.PlayBtn.Value = ~app.PlayBtn.Value;
            app.onPlay(app.PlayBtn.Value);
        case 'm'
            m = 'signed';
            if strcmp(app.MapMode, 'signed'), m = 'magnitude'; end
            app.ModeSwitch.Value = m;  app.onMode(m);
    end
end

% Author: Diellor Basha, 2026
