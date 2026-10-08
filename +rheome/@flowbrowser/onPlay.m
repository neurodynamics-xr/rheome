function onPlay(app, on)
% ONPLAY  Animate forward, rolling into the next page at the end of the core.
%
% ⚠ STOPS AT THE END OF THE CORE, not the end of the span: the margin belongs to the next
% page, so playing through it would show the same seconds twice.
%
% Author: Diellor Basha, 2026

    app.Playing = logical(on);
    while app.Playing && isvalid(app.Fig)
        last = find(app.Page.Core, 1, 'last');
        if app.Frame >= last
            if app.PageIndex >= app.Bundle.pager.NumPages
                app.Playing = false;  app.PlayBtn.Value = false;  break;
            end
            app.PageSpin.Value = app.PageIndex + 1;
            app.onPage(app.PageIndex + 1);
        else
            app.onFrame(app.Frame + 1);
            app.FrameSlider.Value = app.Frame;
        end
        drawnow limitrate;
    end
end

% Author: Diellor Basha, 2026
