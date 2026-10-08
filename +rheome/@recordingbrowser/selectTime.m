function selectTime(app, t)
% SELECTTIME  Select the tile containing a time, and report its place in the pyramid.
%
%   app.selectTime(137.5)        % or click the strip
%
% The key arithmetic IS the hierarchy: the parent of (L, k) is (L+1, floor((k-1)/2)+1) and
% the children are (L-1, 2k-1) and (L-1, 2k). No index is consulted -- which is why a
% descent costs nothing and why a database can do the same join with integer arithmetic.
%
% Author: Diellor Basha, 2026

    g = app.Db.grid;  L = app.Level;
    k = min(max(floor(double(t) / g.tExtent(L+1)) + 1, 1), g.K(L+1));
    s = struct();
    s.level = L;  s.k = k;
    s.t_lo = (k-1) * g.tExtent(L+1);
    s.t_hi = min(s.t_lo + g.tExtent(L+1), app.Db.meta.duration);
    s.t_center = g.tCenter{L+1}(k);
    if L < g.Lmax
        s.parent = [L+1, min(floor((k-1)/2) + 1, g.K(L+2))];
    else
        s.parent = [];
    end
    if L > 0
        s.children = [L-1, 2*k-1; L-1, min(2*k, g.K(L))];
    else
        s.children = [];
    end
    s.unit = app.Unit;  s.scope = app.Scope;          % a tile of WHICH channel: the write key
    app.Selected = s;
    % ⚠ the channels view has no Rows table, so the statistic lookup is skipped rather than
    % indexing an empty table -- which is what a click in the stack used to do.
    i = [];
    if ismember('k', app.Rows.Properties.VariableNames), i = find(app.Rows.k == k, 1); end
    txt = sprintf('%s tile (L%d, k%d)  %.2f-%.2f s', i_unitname(app), L, k, s.t_lo, s.t_hi);
    if ~isempty(i)
        [y, lbl] = i_value(app, i);
        txt = sprintf('%s   %s = %.4g', txt, lbl, y);
    end
    if ~isempty(s.parent), txt = sprintf('%s   parent (L%d, k%d)', txt, s.parent(1), s.parent(2)); end
    if ~isempty(s.children), txt = sprintf('%s   children (L%d, k%d..%d)', txt, s.children(1,1), s.children(1,2), s.children(2,2)); end
    app.setStatus(txt);
    if ~isempty(app.Fig) && isvalid(app.Fig)
        if strcmp(app.View, 'channels'), app.drawChannels(); else, app.drawStat(); end
        app.updateLabels();
    end
end

function nm = i_unitname(app)
    if strcmp(app.Scope, 'group')
        nm = sprintf('node %d', app.Unit);
    else
        nm = char(string(app.Db.meta.ChannelName{app.Unit}));
    end
end

function [y, lbl] = i_value(app, i)
    R = app.Rows;  s = app.Stat;  lbl = s;
    v = R.(s);
    if size(v, 2) > 1
        b = find(R.Properties.CustomProperties.bands == app.Band, 1);
        if isempty(b), y = NaN; return; end
        lbl = sprintf('%s(band %d)', s, app.Band);
        y = double(v(i, b));
    else
        y = double(v(i));
    end
end
% Author: Diellor Basha, 2026
