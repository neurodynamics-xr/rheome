function onBand(app, name)
% ONBAND  Switch band: recompute the page. The page GRID does not change -- the margin was
% sized for the whole master range in rheome.flowpage.prepare -- so only the band's filters and its
% derived rate change.
% Author: Diellor Basha, 2026
    app.BandName = name;
    app.loadPage();
end

% Author: Diellor Basha, 2026
