function drawTiles(app)
% DRAWTILES  The constant-Q strip: one row per band, cut into that band's own tiles.
%
% ⭐ THE STAIRCASE IS THE STRUCTURE. In 'natural' mode every row is drawn at its own level,
% so 64-128 Hz is cut into 0.25 s tiles and 1-2 Hz into 16 s tiles in the same picture, and
% the row label carries the tile length. That is the constant-Q tiling itself: halve the
% frequency, double the support, double the tile, and the number of CYCLES in a tile stays
% the same (22.6 on this bank). In 'level' mode every row is cut at the level the statistic
% above uses, so the columns line up with the axis above and nothing mixes.
%
% ⭐ ONE IMAGE PER BAND, NOT ONE PATCH PER TILE. Each row is a 1 x K image placed at its own
% octave -- fifteen graphics objects for a few thousand tiles, which is what keeps a zoom
% instant. Patching every tile costs K x nBands objects and visibly stalls at 400 tiles.
%
% ⚠ IN 'level' MODE THE EMPTY REGION IS THE DIAGONAL, NOT MISSING DATA. Bands whose support
% is longer than that level's tile are not stored there; they are shaded, and the label says
% which level would carry them. In 'natural' mode there is nothing to shade, because every
% band is drawn at a level that does carry it.
%
% ⚠ THE COMPLETING MEMBERS ARE CLAMPED TO ONE OCTAVE. The bank's 'below' member reaches down
% to a few millihertz and its 'above' member up to Nyquist; drawn literally on a log axis the
% low one swallows the plot and crushes every octave that matters into a sliver. Each is
% given one octave of room at the edge, which is a drawing decision and changes no number.
%
% Author: Diellor Basha, 2026

    ax = app.TileAx;  cla(ax, 'reset');
    S = app.Strip;  b = app.Db.bands;
    if isempty(b)
        % a preview store: there is no spectrum to draw, and saying so is better than an
        % empty pane. The window's band SUPPORT is still meaningful, so it goes here.
        i_nobands(app, ax);  return
    end
    oct = strcmp(string(b.kind), "octave");
    ylo = min(b.fLo(oct));  yhi = max(b.fHi);
    ybot = ylo / 2;  ytop = max(yhi, 2 * max(b.fHi(oct)));       % one octave for each edge member

    % the colour axis spans every row, so rows are comparable with each other
    v = cellfun(@(p) 10*log10(max(p, realmin) / 1e-26), {S.rows.power}, 'UniformOutput', false);
    fin = cellfun(@(x) x(isfinite(x)), v, 'UniformOutput', false);
    allv = [fin{:}];
    if isempty(allv), cl = [0 1]; else, cl = [min(allv) max(allv)]; end
    if diff(cl) <= 0, cl = cl(1) + [0 1]; end

    hold(ax, 'on');
    yt = zeros(1, numel(S.rows));  lbl = cell(1, numel(S.rows));
    for i = 1:numel(S.rows)
        r = S.rows(i);
        k = find(b.j == r.band, 1);
        y = log2([max(b.fLo(k), ybot) min(b.fHi(k), ytop)]);
        image(ax, 'XData', [r.tLo(1) r.tHi(end)], 'YData', y, ...
                  'CData', v{i}, 'CDataMapping', 'scaled');
        yt(i)  = mean(y);
        lbl{i} = sprintf('%.3g Hz · %s', b.fLo(k), i_dur(app.Db.grid.tExtent(r.level+1)));
    end

    % the bands this level does not carry: shaded, so the gap reads as the diagonal
    miss = setdiff(b.j(:)', S.bands);
    for j = miss
        k = find(b.j == j, 1);
        y = log2([max(b.fLo(k), ybot) min(b.fHi(k), ytop)]);
        patch(ax, 'XData', app.Window([1 2 2 1]), 'YData', y([1 1 2 2]), ...
              'FaceColor', [0.90 0.90 0.90], 'EdgeColor', [0.97 0.97 0.97]);
    end
    colormap(ax, parula);  clim(ax, cl);
    if ~isempty(app.CB_) && isvalid(app.CB_), delete(app.CB_); end
    app.CB_ = colorbar(ax, 'Position', [0.948 0.330 0.011 0.265]);   % fixed, so the axes never shift

    ylim(ax, log2([ybot ytop]));
    [yt, i] = sort(yt);  lbl = lbl(i);
    set(ax, 'YTick', yt, 'YTickLabel', lbl, 'FontSize', 8);
    if S.hasBand
        k = find(b.j == app.Band, 1);
        yline(ax, log2(max(b.fLo(k), ybot)), 'Color', [0.9 0.2 0.2], 'LineWidth', 1.2);
        yline(ax, log2(min(b.fHi(k), ytop)), 'Color', [0.9 0.2 0.2], 'LineWidth', 1.2);
    end
    xlim(ax, app.Window);
    ylabel(ax, 'band · tile');  xlabel(ax, 'time (s)');
    title(ax, sprintf('band power per tile, dB re 1e-26 T^2  --  %s', i_caption(app, S, miss, b)), ...
          'Interpreter', 'none', 'FontSize', 9);
    set(ax, 'ButtonDownFcn', @(~, e) i_click(app, e));
end

% What a window of this length could resolve, which is the question a preview is used to
% answer: pick a stretch, read off the bands that fit in it.
function i_nobands(app, ax)
    W = diff(app.Window);
    T = rheome.ingest.support(W);
    text(ax, 0.5, 0.62, 'preview store: min and max only, no bands', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'Color', [0.45 0.45 0.45], 'Interpreter', 'none');
    text(ax, 0.5, 0.38, sprintf('this %.3g s window supports %.3g Hz and above (lowest full octave %.3g-%.3g Hz, %.0f cycles)', ...
         W, T.f_lowest, T.band_lo, T.band_hi, T.cycles_in_window), 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'Color', [0.25 0.25 0.25], 'Interpreter', 'none');
    set(ax, 'XTick', [], 'YTick', []);
    xlim(ax, app.Window);
    title(ax, 'band support of the selected window  --  rheome.ingest.support', 'Interpreter', 'none', 'FontSize', 9);
end

function s = i_caption(app, S, miss, b)
    if strcmp(S.mode, 'natural')
        lv = unique(S.levels);
        s = sprintf('each band at its own tile: levels %d-%d, %g s to %g s, %.1f cycles per tile', ...
            min(lv), max(lv), app.Db.grid.tExtent(min(lv)+1), app.Db.grid.tExtent(max(lv)+1), ...
            median(app.Ladder.cycles_per_tile(app.Ladder.kind == "octave")));
    elseif isempty(miss)
        s = sprintf('every band at level %d, %g s tiles', S.level, app.Db.grid.tExtent(S.level+1));
    else
        s = sprintf('level %d only: %d band(s) need a coarser level (from %d up)', S.level, ...
            numel(miss), min(b.naturalLevel(ismember(b.j, miss))));
    end
end

% Seconds throughout, never minutes: the tile lengths are dyadic in seconds (0.25, 0.5, 1,
% … 1024) and converting the long ones to minutes turns a clean ladder into 1.06667m.
function s = i_dur(t)
    s = sprintf('%gs', t);
end

function i_click(app, e)
    if isprop(e, 'IntersectionPoint'), app.selectTime(e.IntersectionPoint(1)); end
end
% Author: Diellor Basha, 2026
