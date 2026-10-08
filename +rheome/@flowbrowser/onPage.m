function onPage(app, k)
% ONPAGE  Advance to another page of the recording (the slow step: a fresh CWT).
% Author: Diellor Basha, 2026
    app.PageIndex = round(k);
    app.loadPage();
end

% Author: Diellor Basha, 2026
