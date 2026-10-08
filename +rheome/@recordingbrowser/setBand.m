function setBand(app, b)
% SETBAND  Choose the band -- and, with FollowBand, let the band choose the time window.
%
%   app.setBand(6)        % 8-16 Hz: level 3, 2 s tiles, a window of TilesPerView of them
%
% ⭐ THE BAND IS THE NATURAL HANDLE, not the level. A band fixes its own support, so it
% fixes the shortest tile that can hold it and therefore the level; asking to look at
% 8-16 Hz is asking for 2 s tiles, and the window that shows a useful number of them. This
% is the level-of-detail rule read backwards, and it is why the ladder (rheome.select.ladder) is
% the same table either way.
%
% ⚠ THE CHANNEL FLOOR STILL WINS. On a store whose channel rows start at 16 s, picking an
% 8-16 Hz band in channel scope gives 16 s tiles, not 2 s: the band's own level exists only
% for sensor-tree nodes. The readout says which constraint bit.
%
% With FollowBand off the band only moves the strip's cursor and, for a per-band statistic,
% the floor.
%
% Author: Diellor Basha, 2026

    b = double(b);
    if isempty(app.Db.bands)
        error('recordingbrowser:band', 'This is a preview store: it has no bands (rheome.ingest.preview).');
    end
    if ~ismember(b, app.Db.bands.j)
        error('recordingbrowser:band', 'No band %d in this store (bands are %s).', b, mat2str(app.Db.bands.j(:)'));
    end
    app.Band = b;
    if ~app.FollowBand
        app.refresh();
        return
    end
    g = app.Db.grid;
    j = find(app.Db.bands.j == b, 1);
    base = 0;
    if ~strcmp(app.Scope, 'group'), base = app.Db.channelLevel; end
    L = min(max(app.Db.bands.naturalLevel(j), base), g.Lmax);
    w = min(app.TilesPerView * g.tExtent(L+1), app.Db.meta.duration);
    c = mean(app.Window);
    app.Auto  = false;                       % the band owns the level now
    app.Level = L;
    app.Window = i_clamp(c - w/2, w, app.Db.meta.duration);
    app.refresh();
    app.setStatus(sprintf('band %d (%.4g-%.4g Hz): support %.3g s -> %g s tiles at level %d, %.1f cycles per tile', ...
        b, app.Db.bands.fLo(j), app.Db.bands.fHi(j), app.Db.bands.tSupport(j), g.tExtent(L+1), L, ...
        g.tExtent(L+1) * app.Db.bands.fCenter(j)));
end

function w = i_clamp(t0, width, D)
    if width >= D, w = [0 D]; return; end
    t0 = min(max(t0, 0), D - width);
    w = [t0, t0 + width];
end
% Author: Diellor Basha, 2026
