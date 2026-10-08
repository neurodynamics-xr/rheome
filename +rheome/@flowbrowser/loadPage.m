function loadPage(app)
% LOADPAGE  Compute the current (band, page). The only slow step -- ~2-3 s.
% Author: Diellor Basha, 2026

    band = app.Bands.(app.BandName);
    app.setStatus(sprintf('computing %s, page %d …', app.BandName, app.PageIndex));
    drawnow limitrate;

    t0 = tic;
    app.Page = rheome.flowpage(app.Bundle, 'Band', band, 'PageIndex', app.PageIndex, ...
        'FreqLimits', app.Bundle.freqLimits, 'VoicesPerOctave', app.Bundle.voicesOct);
    dt = toc(t0);

    % land on the first CORE frame -- the margin is not this page's to show by default
    iCore = find(app.Page.Core, 1);
    app.FrameSlider.Limits = [1 max(2, app.Page.NumFrames)];

    app.ScalogramData = scalogram(app.Page);      % once per page; every frame reads it
    app.buildCortexAxes();
    app.drawSensors();
    app.drawTimeline();
    app.Frame = 0;                       % force setFrame to do the work
    app.setFrame(max(1, iCore));
    app.setStatus(sprintf('%s  %d filters  %.0f Hz  %d frames  (%.1f s)', ...
        app.BandName, app.Page.NumFilters, app.Page.Rate, app.Page.NumFrames, dt));
end

% Author: Diellor Basha, 2026
