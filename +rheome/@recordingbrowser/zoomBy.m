function zoomBy(app, f)
% ZOOMBY  Scale the window about its centre; f < 1 zooms in, f > 1 zooms out.
%
% The level follows from the new width through the level-of-detail rule, so zooming in twice
% by 2 is exactly one coarser level's worth of detail gained -- the view and the store stay
% in step without either one being told about the other.
%
% Author: Diellor Basha, 2026

    f = double(f);
    if ~(f > 0), error('recordingbrowser:zoom', 'Zoom factor must be positive.'); end
    c = mean(app.Window);  w = diff(app.Window) * f;
    app.setWindow(c - w/2, c + w/2);
end
% Author: Diellor Basha, 2026
