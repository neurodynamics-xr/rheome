function stepBand(app, d)
% STEPBAND  Move one band up or down the ladder, carrying the window with it.
%
% Stepping down an octave doubles the tile and doubles the window, so the view keeps the
% same number of tiles and the same number of cycles on screen while the frequency halves.
% That is the constant-Q tiling made navigable.
%
% Author: Diellor Basha, 2026

    js = app.Db.bands.j(:)';
    i = find(js == app.Band, 1);
    if isempty(i), i = 1; end
    app.setBand(js(min(max(i + round(d), 1), numel(js))));
end
% Author: Diellor Basha, 2026
