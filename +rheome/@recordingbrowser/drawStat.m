function drawStat(app)
% DRAWSTAT  The statistic on the top axis, one STEP per tile.
%
% ⚠ STAIRS, NOT A LINE. A tile's value holds over its whole extent; joining tile centres
% with a line draws a slope the data does not have and puts the transition in the wrong
% place -- at 64 s tiles that is a 32 s lie. The steps are the tile edges.
%
% Author: Diellor Basha, 2026

    ax = app.StatAx;  cla(ax, 'reset');
    R = app.Rows;
    if app.Envelope, i_envelope(app, ax, R); return; end
    [y, lbl] = i_column(app, R);
    e = [R.t_lo; R.t_hi(end)];                                  % tile edges: K+1 of them
    stairs(ax, e, [y; y(end)], 'Color', [0.10 0.25 0.55], 'LineWidth', 1.4);
    hold(ax, 'on');
    if any(isfinite(y))
        i = find(y == max(y(isfinite(y))), 1);
        plot(ax, R.t_center(i), y(i), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0.8 0.1 0.1]);
    end
    if isfield(app.Selected, 'k') && ~isempty(app.Selected.k) && app.Selected.level == app.Level
        xline(ax, app.Selected.t_center, 'Color', [0.8 0.1 0.1], 'LineWidth', 0.8);
    end
    grid(ax, 'on');
    xlim(ax, app.Window);
    yl = i_ylim(y);
    if ~isempty(yl), ylim(ax, yl); end
    ylabel(ax, lbl, 'Interpreter', 'none');
    title(ax, sprintf('%s  --  %s at level %d (%g s tiles), %d tiles in view', ...
        app.Name, lbl, app.Level, app.Db.grid.tExtent(app.Level+1), height(R)), 'Interpreter', 'none');
    set(ax, 'XTickLabel', []);
end

% ⭐ THE ENVELOPE. A filled band between the tile's min and max, one column per tile, at
% whatever level the zoom picked. Because both merge by extremum, the band CONTAINS every
% sample of the span it covers: it can be too wide, never too narrow. When raw samples are
% already in hand (a short window) they are drawn over it, which is the proof -- the trace
% must stay inside the band everywhere.
function i_envelope(app, ax, R)
    lo = double(R.min(:));  hi = double(R.max(:));
    x = [R.t_lo(:)'; R.t_hi(:)'];                            % two x per tile: a flat step
    % ⚠ repelem OF A SCALAR RETURNS A ROW, so a one-tile window (a long tile zoomed into,
    % which is the normal end of a zoom) gave the fill an x and a y of different shapes and
    % threw inside a callback. reshape(repmat(...)) is shape-stable for K = 1 and K > 1 alike.
    xv = x(:);
    dbl = @(v) reshape(repmat(v(:)', 2, 1), [], 1);          % each value twice: a flat step
    fill(ax, [xv; flipud(xv)], [dbl(lo); flipud(dbl(hi))], [0.42 0.58 0.80], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.85);
    hold(ax, 'on');
    if ~isempty(app.Raw.x)
        plot(ax, app.Raw.t, app.Raw.x, 'Color', [0.15 0.15 0.15], 'LineWidth', 0.5);
    end
    grid(ax, 'on');  xlim(ax, app.Window);
    yl = i_ylim([lo; hi]);
    if ~isempty(yl), ylim(ax, yl); end
    ylabel(ax, sprintf('min/max (%s)', app.Db.meta.units));
    g = app.Db.grid;
    extra = '';
    if ~isempty(app.Raw.x), extra = sprintf(', %d raw samples drawn over it', numel(app.Raw.x)); end
    title(ax, sprintf('%s  --  min/max envelope at level %d (%g s tiles), %d tiles in view%s', ...
        app.Name, app.Level, g.tExtent(app.Level+1), height(R), extra), 'Interpreter', 'none');
    set(ax, 'XTickLabel', []);
end

% A single tile in view gives a degenerate range, and MATLAB then pads it to something like
% [-1 1.5] around a value of 0.47 -- which reads as a plot of nothing. Pad by a tenth of the
% range instead, or by a twentieth of the value when the range is zero.
function yl = i_ylim(y)
    y = y(isfinite(y));
    if isempty(y), yl = []; return; end
    lo = min(y);  hi = max(y);
    if hi > lo
        m = 0.1 * (hi - lo);
    else
        m = max(abs(lo) / 20, realmin);
    end
    yl = [lo - m, hi + m];
end

% A per-band statistic is a matrix column per band; the selected band picks one, and a band
% the level does not carry is reported rather than silently drawn as something else.
function [y, lbl] = i_column(app, R)
    s = app.Stat;
    if ismember(s, R.Properties.VariableNames) && size(R.(s), 2) == 1
        y = double(R.(s));  lbl = s;  return
    end
    bAt = R.Properties.CustomProperties.bands;
    i = find(bAt == app.Band, 1);
    j = find(app.Db.bands.j == app.Band, 1);
    nm = sprintf('%s [%.3g-%.3g Hz]', s, app.Db.bands.fLo(j), app.Db.bands.fHi(j));
    if isempty(i)
        y = nan(height(R), 1);  lbl = [nm '  (not carried here)'];  return
    end
    M = double(R.(s));
    y = M(:, i);  lbl = nm;
end
% Author: Diellor Basha, 2026
