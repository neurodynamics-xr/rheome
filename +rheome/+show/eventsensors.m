function [hFig, info] = eventsensors(X, t, E, pos2, varargin)
% SHOW.EVENTSENSORS  A tracked event on the sensor traces and on a sensor topography.
%
%   [hFig, info] = rheome.show.eventsensors(X, t, E, pos2)
%   rheome.show.eventsensors(X, t, E, pos2, 'NumChannels', 8, 'Labels', names, 'Visible', 'off')
%       X     [C x nS] sensor data (the traces shown)      t [1 x nS] times (s)
%       E     rheome.detect.eventsensors output             pos2 [C x 2] flattened sensor positions
%
% Left: the NumChannels channels with the largest event footprint, stacked; trace-neutral (#8a94a3) =
% the recording, trace-event (#d1495b, 1.5x width, on the same line) = the samples where that channel
% CARRIES the event (E.mask). No shaded spans and no second data line on the trace (house rule,
% library/styles/cortical-flow-visual-language.md §5c); E.Xevent stays in E. Right: E.footprint on the array,
% carrying channels ringed and the traced ones numbered. So a reader can find in the raw traces what
% the tracker followed on the cortex.
%
% OPTIONS: 'NumChannels' (8)  'Labels' {1 x C} channel names  'Color' (trace-event [0.820 0.286 0.357])
%          'Visible' ('on')  'Title' ('')
% OUTPUT:  info .channels [NumChannels x 1] the traced channels   .carrying [n x 1] union over steps
%
% See also: rheome.detect.eventsensors, rheome.detect.tiletrack
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('NumChannels', 8);
    p.addParameter('Labels', {});
    p.addParameter('Color', [209 73 91] / 255);     % trace-event #d1495b (§5c)
    p.addParameter('Visible', 'on');
    p.addParameter('Title', '');
    p.parse(varargin{:});  o = p.Results;

    X = real(double(X));  t = double(t(:)).';  C = size(X, 1);
    [~, ord] = sort(E.footprint, 'descend');  ch = ord(1:min(o.NumChannels, C));
    carrying = unique(vertcat(E.channels{:}));
    lab = o.Labels;  if isempty(lab), lab = arrayfun(@(c) sprintf('ch %d', c), 1:C, 'uni', 0); end

    hFig = figure('Visible', o.Visible, 'Color', 'w', 'Position', [100 100 1100 480]);
    tl = tiledlayout(hFig, 1, 3, 'TileSpacing', 'compact');
    if ~isempty(o.Title), title(tl, o.Title); end
    ax = nexttile(tl, [1 2]);  hold(ax, 'on');
    sp = 3 * max(std(X(ch, :), 0, 2));  yl = [-sp, numel(ch) * sp];
    for i = 1:numel(ch)
        c = ch(i);  off = (numel(ch) - i) * sp;
        plot(ax, t, X(c, :) + off, 'Color', [138 148 163] / 255, 'LineWidth', 1);   % trace-neutral #8a94a3
        y = X(c, :) + off;  y(~E.mask(c, :)) = NaN;
        plot(ax, t, y, 'Color', o.Color, 'LineWidth', 1.5);
        text(ax, t(1), off, sprintf('%d  %s', i, lab{c}), 'HorizontalAlignment', 'right', 'FontSize', 8);
    end
    set(ax, 'YTick', [], 'YLim', yl);  xlabel(ax, 'time (s)');  box(ax, 'off');
    title(ax, 'sensor traces: red = the samples the tracked pattern occupies on that sensor');

    ax2 = nexttile(tl);  hold(ax2, 'on');  axis(ax2, 'equal', 'off');
    scatter(ax2, pos2(:, 1), pos2(:, 2), 36, E.footprint, 'filled');
    scatter(ax2, pos2(carrying, 1), pos2(carrying, 2), 70, o.Color, 'LineWidth', 1.2);
    text(ax2, pos2(ch, 1), pos2(ch, 2), compose(' %d', (1:numel(ch))'), 'FontSize', 8);
    colormap(ax2, parula);  cb = colorbar(ax2);  cb.Label.String = 'event footprint (mean |x_{ev}|^2)';
    title(ax2, 'topography: ringed = carrying channels');
    info = struct('channels', ch, 'carrying', carrying);
end

% Author: Diellor Basha, 2026
