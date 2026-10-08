function C = context(db, M, opts)
% SELECT.CONTEXT  Join measurements to the mergeable matrix: the state each one occurred in.
%
%   C = rheome.select.context(db, M)
%   C = rheome.select.context(db, M, Stats=["rms","share"], Ancestors=[0 2 4])
%
% M is a measurement table from rheome.select.measures (or any table with level, k and, for a unit
% scope, unit_id). Each row comes back with the derived statistics of ITS OWN node, and, for
% each offset in Ancestors, of the node that many levels above it. Columns from an ancestor
% are suffixed `_upN`.
%
% ⭐ THIS IS WHERE THE TWO MATRICES MEET. The measurement says what happened and where; the
% moments say what the signal was doing there, at that scale and at any coarser one. So
% "group the vortices by the alpha amplitude they occurred in" is a join, not a second pass
% over the data, and "...and by the state of the minute around them" is the same join at
% Ancestors = 4.
%
% ⚠ ONLY THE CONTEXT CROSSES LEVELS, NEVER THE MEASUREMENT. The ancestor's statistics are
% merged quantities and are exact at that level; the measurement itself stays where it was
% made. Nothing here implies a measurement could be rolled up.
%
% A row whose band_id is set picks that band's column out of a per-band statistic (share,
% bandPower), so a vortex found in one band is scored against its own band by default.
%
% ⚠ NaN MEANS THE CONTEXT DOES NOT EXIST, NOT THAT IT IS ZERO. A band below its natural
% level at that node, or a channel below the store's channel floor, has no row to join to;
% the column is left NaN rather than filled from a neighbouring band or level. Raising the
% Ancestors offset is usually what makes it appear.
%
% See also: rheome.select.measures, rheome.select.derive, rheome.select.ladder
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        M table
        opts.Stats     string = ["rms","envMax","share","spectralCentroid"]
        opts.Ancestors double = 0
    end
    if isempty(M), C = M; return; end
    need = {'level','k'};
    if ~all(ismember(need, M.Properties.VariableNames))
        error('select:context:rows', 'M needs level and k.');
    end
    C = M;
    g = db.grid;
    if ~ismember('unit_id', C.Properties.VariableNames), C.unit_id = zeros(height(C), 1); end
    if ~ismember('scope', C.Properties.VariableNames), C.scope = repmat("none", height(C), 1); end
    if ~ismember('band_id', C.Properties.VariableNames), C.band_id = zeros(height(C), 1); end

    for d = unique(double(opts.Ancestors(:)))'
        if d < 0, error('select:context:ancestor', 'An ancestor offset must be >= 0.'); end
        suffix = '';  if d > 0, suffix = sprintf('_up%d', d); end
        cols = i_blank(opts.Stats, height(C));
        % one derive call per (scope, unit, ancestor level): the levels cache in the handle
        key = strcat(C.scope, "|", string(C.unit_id), "|", string(min(C.level + d, g.Lmax)));
        for kk = unique(key)'
            r = find(key == kk);
            L = min(C.level(r(1)) + d, g.Lmax);
            [ok, R] = i_derive(db, L, opts.Stats, C.scope(r(1)), C.unit_id(r(1)));
            if ~ok, continue; end
            ka = i_ancestor(C.k(r), d, g, L);
            [tf, pos] = ismember(ka, R.k);
            for s = opts.Stats
                v = R.(char(s));
                if size(v, 2) > 1                                  % a per-band column
                    bAt = R.Properties.CustomProperties.bands;
                    for i = find(tf(:))'
                        b = C.band_id(r(i));
                        j = find(bAt == b, 1);
                        if ~isempty(j), cols.(char(s))(r(i)) = v(pos(i), j); end
                    end
                else
                    cols.(char(s))(r(tf)) = double(v(pos(tf)));
                end
            end
        end
        for s = opts.Stats
            C.(sprintf('%s%s', char(s), suffix)) = cols.(char(s));
        end
        if d > 0
            C.(sprintf('k%s', suffix)) = i_ancestor(C.k, d, g, []);
            C.(sprintf('level%s', suffix)) = min(C.level + d, g.Lmax);
        end
    end
end

% The ancestor key is arithmetic, not a lookup: d levels up, k -> floor((k-1)/2^d) + 1.
function ka = i_ancestor(k, d, g, L)
    ka = floor((k - 1) / 2^d) + 1;
    if ~isempty(L), ka = min(ka, g.K(L + 1)); end
end

% A unit that carries no rows at this level (the diagonal) leaves NaNs rather than throwing:
% a context join must not fail because one row sits below the channel floor.
function [ok, R] = i_derive(db, L, stats, scope, unit)
    ok = true;  R = table();
    try
        if scope == "group"
            R = rheome.select.derive(db, Level=L, Stats=stats, Scope="group", Groups=unit);
        elseif scope == "channel"
            R = rheome.select.derive(db, Level=L, Stats=stats, Channels=unit);
        else
            R = rheome.select.derive(db, Level=L, Stats=stats, Channels=1);   % no unit: the first channel
        end
    catch
        ok = false;
    end
end

function s = i_blank(stats, n)
    s = struct();
    for k = stats, s.(char(k)) = nan(n, 1); end
end
% Author: Diellor Basha, 2026
