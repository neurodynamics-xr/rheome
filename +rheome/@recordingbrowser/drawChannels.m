function drawChannels(app)
% DRAWCHANNELS  The stack: one min/max band per selected channel, on a shared scale.
%
% ⭐ THIS IS THE PICTURE OF THE STORED MATRIX. Every row is exact at the level the zoom
% picked: the band contains every sample of the span it covers, because min and max merge by
% extremum and are stored as outward-rounded bounds. Nothing here is decimated or smoothed,
% and nothing here reads the recording.
%
% ⚠ ONE SCALE FOR ALL ROWS, so the rows stay comparable: a quiet channel looks quiet, and a
% loud one looks loud. Per-row normalisation would hide exactly what the view is for.
%
% ⚠ BUT NOT THE WIDEST ROW'S SCALE. One bad sensor is often ten times the rest, and scaling
% to it flattens every other row to a line -- measured on the reference data, where MLC11 spans 1.05e-11 T
% against a typical 1e-12. The scale is the 75th percentile of the drawn ranges times Gain,
% and rows that exceed their lane are CLIPPED and counted in the title, which is what a trace
% viewer has always done.
%
% Author: Diellor Basha, 2026

    ax = app.StatAx;  cla(ax, 'reset');
    S = app.Stack;
    n = numel(S.units);
    if n == 0
        text(ax, 0.5, 0.5, 'no channels selected', 'Units', 'normalized', ...
             'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45]);
        return
    end
    lo = double(S.min);  hi = double(S.max);                 % [K x n]
    mid = (max(hi, [], 1) + min(lo, [], 1)) / 2;             % each row about its own centre
    rng = max(hi, [], 1) - min(lo, [], 1);
    span = max(i_p75(rng), realmin);
    step = 1;                                                % rows one unit apart
    lane = 0.46 * step;                                      % half a lane, so rows never touch
    sc = app.Gain * lane / (span / 2);
    clipped = sum(rng > span * 1.001);

    x  = [S.tLo(:)'; S.tHi(:)'];  xv = x(:);
    dbl = @(v) reshape(repmat(v(:)', 2, 1), [], 1);
    clip = @(v) max(min(v, lane), -lane);
    hold(ax, 'on');
    col = [0.42 0.58 0.80];  focus = [0.85 0.33 0.20];
    for i = 1:n
        y0 = (n - i) * step;
        c = col;  if S.units(i) == app.Unit, c = focus; end
        a = clip(sc * (lo(:, i) - mid(i)));  b = clip(sc * (hi(:, i) - mid(i)));
        fill(ax, [xv; flipud(xv)], y0 + [dbl(a); flipud(dbl(b))], c, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.9);
    end
    % ⚠ THE LOADED TRACE BELONGS TO ITS OWN CHANNEL, not to whichever row has the focus now.
    % Drawing it on the focused row would relabel someone else's samples, which is the one
    % thing a viewer of a stored matrix must never do.
    if ~isempty(app.Raw.x) && isfield(app.Raw, 'unit')
        i = find(S.units == app.Raw.unit, 1);
        if ~isempty(i)
            y0 = (n - i) * step;
            plot(ax, app.Raw.t, y0 + clip(sc * (app.Raw.x - mid(i))), 'Color', [0.1 0.1 0.1], 'LineWidth', 0.5);
        end
    end
    % measurements stored on these tiles, as ticks on their own rows
    nm = i_marks(app, ax, S, n, step);
    set(ax, 'YTick', (0:n-1) * step, 'YTickLabel', flip(i_names(app, S.units)), 'FontSize', 8);
    ylim(ax, [-0.6*step, (n-1)*step + 0.6*step]);
    xlim(ax, app.Window);
    grid(ax, 'on');  set(ax, 'XTickLabel', [], 'YGrid', 'off');
    g = app.Db.grid;
    extra = '';
    if numel(app.Units) > n, extra = sprintf(' of %d selected', numel(app.Units)); end
    cl = '';  if clipped > 0, cl = sprintf(', %d clipped', clipped); end
    if nm > 0, cl = sprintf('%s, %d measured', cl, nm); end
    title(ax, sprintf('%s  --  min/max envelope of %d channel%s%s at level %d (%g s tiles), %d tiles, lane %.3g %s%s', ...
        app.Name, n, i_plural(n), extra, app.Level, g.tExtent(app.Level+1), numel(S.k), ...
        span / app.Gain, app.Db.meta.units, cl), 'Interpreter', 'none', 'FontSize', 9);
    set(ax, 'ButtonDownFcn', @(~, e) i_click(app, e));
end

% ⭐ WHAT IS STORED ON THESE TILES. A measurement is keyed by (scope, unit, level, k), so a
% mark belongs to one row and one tile; it is drawn at the tile's centre on that row's top
% edge. Measurements made at another level still show, placed by the span they cover, which
% is how a coarse note stays visible while you are zoomed in.
function n = i_marks(app, ax, S, nRows, step)
    n = 0;
    M = app.Marks;
    if isempty(M) || isempty(S.units), return; end
    w = app.Window;  g = app.Db.grid;
    ext = g.tExtent(M.level + 1)';
    lo = (M.k - 1) .* ext;  hi = min(lo + ext, app.Db.meta.duration);
    keep = hi > w(1) & lo < w(2) & ismember(M.unit_id, S.units) & M.scope == string(app.Scope);
    M = M(keep, :);  lo = lo(keep);  hi = hi(keep);
    n = height(M);
    if n == 0, return; end
    for i = 1:n
        r = find(S.units == M.unit_id(i), 1);
        y = (nRows - r) * step + 0.47 * step;
        plot(ax, [lo(i) hi(i)], [y y], 'Color', [0.10 0.45 0.15], 'LineWidth', 2);
        plot(ax, mean([lo(i) hi(i)]), y, 'v', 'MarkerSize', 4, ...
             'MarkerFaceColor', [0.10 0.45 0.15], 'MarkerEdgeColor', 'none');
    end
end

function s = i_plural(n), if n == 1, s = ''; else, s = 's'; end, end

% The 75th percentile of the drawn ranges: robust to a handful of bad sensors, and still a
% real number from the data rather than a constant. MATLAB's prctile is Statistics Toolbox,
% which this tree does not use.
function y = i_p75(x)
    x = sort(double(x(:)));  n = numel(x);
    if n == 0, y = 0; return; end
    pos = 100 * ((1:n)' - 0.5) / n;
    if 75 <= pos(1), y = x(1); elseif 75 >= pos(end), y = x(end); else, y = interp1(pos, x, 75); end
end

function nm = i_names(app, units)
    if strcmp(app.Scope, 'group')
        nm = arrayfun(@(u) sprintf('node %d', u), units, 'UniformOutput', false);
    else
        nm = cellstr(string(app.Db.meta.ChannelName(units)));
    end
end

% Clicking the stack focuses the channel under the pointer and selects the tile under it,
% so the raw load and the detail view follow the eye.
function i_click(app, e)
    if ~isprop(e, 'IntersectionPoint'), return; end
    p = e.IntersectionPoint;
    n = numel(app.Stack.units);
    i = min(max(round(n - p(2)), 1), n);
    app.Unit = app.Stack.units(i);
    app.selectTime(p(1));
end
% Author: Diellor Basha, 2026
