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
% OPTIONS: 'NumChannels' (8)  'Labels' {1 x C} channel names  'Color' (trace-event, figure tokens)
%          'Unit' ('a.u.') unit of X, for the footprint colourbar  'Visible' ('on')
%          'Title' ('') the figure's name (the legend carries it)  'Export' file stem: <stem>.pdf/.svg/.png
% Sizes and colours are the figure tokens (library/styles/tokens/figure-tokens.json), 190 x 70 mm.
% OUTPUT:  info .channels [NumChannels x 1] the traced channels   .carrying [n x 1] union over steps
%
% See also: rheome.detect.eventsensors, rheome.detect.tiletrack
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('NumChannels', 8);
    p.addParameter('Labels', {});
    p.addParameter('Color', T.rgb.trace_event);     % trace-event (§5c)
    p.addParameter('Unit', 'a.u.');
    p.addParameter('Visible', 'on');
    p.addParameter('Title', '');
    p.addParameter('Export', '');
    p.parse(varargin{:});  o = p.Results;

    X = real(double(X));  t = double(t(:)).';  C = size(X, 1);
    [~, ord] = sort(E.footprint, 'descend');  ch = ord(1:min(o.NumChannels, C));
    carrying = unique(vertcat(E.channels{:}));
    lab = o.Labels;  if isempty(lab), lab = arrayfun(@(c) sprintf('ch %d', c), 1:C, 'uni', 0); end

    F = nxr_figure('double', 70, {'AAB'}, 'Kinds', struct('B', 'image'), ...
        'Titles', struct('A', 'Sensor traces, event samples in colour', 'B', 'Event footprint, carrying sensors ringed'));
    hFig = F.fig;  hFig.Visible = o.Visible;
    ax = F.ax.A;  hold(ax, 'on');
    sp = 3 * max(std(X(ch, :), 0, 2));  if sp <= 0, sp = 1; end
    for i = 1:numel(ch)
        c = ch(i);  off = (numel(ch) - i) * sp;
        plot(ax, t, X(c, :) + off, 'Color', T.rgb.trace_neutral, 'LineWidth', T.line.data_pt);
        y = X(c, :) + off;  y(~E.mask(c, :)) = NaN;
        plot(ax, t, y, 'Color', o.Color, 'LineWidth', T.line.event_pt);
        text(ax, t(1), off, sprintf('%d  %s ', i, lab{c}), 'HorizontalAlignment', 'right', 'Interpreter', 'none');
    end
    nxr_axes_style(ax, 'time (s)');
    set(ax, 'YLim', [-sp, numel(ch) * sp], 'XLim', [t(1) t(end)]);  ax.YAxis.Visible = 'off';

    ax2 = F.ax.B;  hold(ax2, 'on');  axis(ax2, 'equal', 'off');
    ms = T.marker.size_pt;
    scatter(ax2, pos2(:, 1), pos2(:, 2), ms^2, E.footprint, 'filled');
    scatter(ax2, pos2(carrying, 1), pos2(carrying, 2), (2*ms)^2, o.Color, 'LineWidth', T.marker.edge_pt);
    text(ax2, pos2(ch, 1), pos2(ch, 2), compose(' %d', (1:numel(ch))'));
    hi = max(E.footprint);  if ~(hi > 0), hi = 1; end
    nicebar(nxr_colorbar(ax2, 'event footprint, mean |x|²', [o.Unit '²'], 'magnitude', [0 hi]), ax2);
    info = struct('channels', ch, 'carrying', carrying);
    finish(F, o.Title, o.Export);
end

% Author: Diellor Basha, 2026
