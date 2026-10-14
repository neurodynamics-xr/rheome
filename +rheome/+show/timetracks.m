function [hFig, V] = timetracks(t, tracks, varargin)
% SHOW.TIMETRACKS  A pattern's parameters over time, one track per row on a shared time axis.
%
%   hFig = rheome.show.timetracks(t, tracks)
%   [hFig, V] = rheome.show.timetracks(t, tracks, 'Width', 'single', 'Export', stem, 'Visible', 'off')
%
% What a travelling, rotating or standing pattern DOES over time: where its extremum is, its scale, its heading, its
% rotation and chirality, its divergence and curl, its amplitude, each a row. A row holds the PLANT (the analytic
% course an analytic field was made with: rheome.show.track, the Editor's track) as a line and/or the DETECTED values
% (what an analysis measured from the data) as markers, so the same figure shows a planted pattern, its detection, and
% how far one is from the other. The key (plant / detected) is in the panel's key band; every row's axis says its
% quantity and unit.
%
% INPUTS:
%   t       [1 x nT] the time axis (s); a plant's support is [t(1) t(end)]
%   tracks  struct array (or cell of structs), one row each, with
%     Label   the parameter's name on its axis ("heading", "scale", "divergence") REQUIRED
%     Unit    its unit ("rad", "mm", "1/s"; "sign" for a chirality)          REQUIRED
%     Kind    'value' | 'angle' (wrapped to (-pi, pi], drawn without the wrap's jump, ticks -pi 0 pi) |
%             'sign' (+1 / -1 steps: chirality, counter-clockwise +1) (default 'value'; 'angle' when Unit is "rad"
%             and Label names an angle: heading, bearing, tilt, direction, phase)
%     Spec    an analytic track (rheome.show.track): the plant, evaluated over t
%     Values  [1 x n] measured values (NaN: none at that time), at Time
%     Time    [1 x n] the times of Values (default t)
%
% OPTIONS (name/value):
%   'Width'   'single' (default) | 'onehalf' | 'double'
%   'RowHeight' mm per row (default 14)
%   'Title'   the figure's name (not drawn: a journal figure's title is its legend's first sentence)
%   'Export'  file stem: <stem>.pdf/.svg/.png at print size (nxr_save)
%   'Visible' 'on' (default) | 'off'
%
% OUTPUT:
%   hFig  the figure;  V  struct array per row: .label .unit .kind .plant (on t, [] without Spec) .t .values
%
% Example: a vortex whose heading turns at 2 Hz while it grows and then reverses its chirality
%   t = linspace(0, 1, 500);
%   tr = struct('Label', {'heading', 'scale', 'chirality'}, 'Unit', {'rad', 'mm', 'sign'}, ...
%       'Spec', {struct('Start', 0, 'To', 4*pi, 'Profile', 'ramp'), struct('Start', 10, 'To', 30, 'Profile', 'ease'), ...
%                struct('Keys', struct('x', {0, 0.5}, 'v', {1, -1}, 'seg', {'hold', 'hold'}))});
%   rheome.show.timetracks(t, tr);
%
% See also: rheome.show.track, rheome.show.phasor, rheome.show.levelstack
%
% Author: Diellor Basha, 2026
    T = tokens();
    p = inputParser;
    p.addParameter('Width', 'single');
    p.addParameter('RowHeight', 14, @(x) isscalar(x) && x > 0);
    p.addParameter('Title', 'Time-tracks');
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;

    t = double(t(:)');
    if numel(t) < 2 || any(diff(t) <= 0), error('rheome:show:timetracks', 't: increasing times (s), at least two.'); end
    if isstruct(tracks), tracks = num2cell(tracks); end
    n = numel(tracks);
    if n == 0, error('rheome:show:timetracks', 'tracks: at least one.'); end
    angleWords = {'heading', 'bearing', 'tilt', 'direction', 'phase', 'angle'};
    V = struct('label', {}, 'unit', {}, 'kind', {}, 'plant', {}, 't', {}, 'values', {});
    for i = 1:n
        a = tracks{i};
        if ~isfield(a, 'Label') || isempty(a.Label) || ~isfield(a, 'Unit') || isempty(a.Unit)
            error('rheome:show:timetracks', 'tracks{%d}: Label and Unit are required (an axis says its quantity and unit).', i);
        end
        kind = 'value';
        if strcmp(a.Unit, 'rad') && any(contains(lower(a.Label), angleWords)), kind = 'angle'; end
        if strcmp(a.Unit, 'sign'), kind = 'sign'; end
        if isfield(a, 'Kind') && ~isempty(a.Kind), kind = char(a.Kind); end
        if ~any(strcmp(kind, {'value', 'angle', 'sign'})), error('rheome:show:timetracks', 'tracks{%d}.Kind: value, angle or sign.', i); end
        plant = [];
        if isfield(a, 'Spec') && ~isempty(a.Spec), plant = rheome.show.track(a.Spec, (t - t(1)) / (t(end) - t(1))); end
        val = [];  tv = t;
        if isfield(a, 'Values') && ~isempty(a.Values)
            val = double(a.Values(:)');
            if isfield(a, 'Time') && ~isempty(a.Time), tv = double(a.Time(:)'); end
            if numel(tv) ~= numel(val), error('rheome:show:timetracks', 'tracks{%d}: Values and Time differ in length.', i); end
        end
        if isempty(plant) && isempty(val), error('rheome:show:timetracks', 'tracks{%d}: a Spec (plant) or Values (detected).', i); end
        V(i) = struct('label', char(a.Label), 'unit', char(a.Unit), 'kind', kind, 'plant', plant, 't', tv, 'values', val);
    end

    hasPlant = any(~cellfun(@isempty, {V.plant}));  hasDet = any(~cellfun(@isempty, {V.values}));
    H = n * o.RowHeight + 22;
    F = nxr_figure(o.Width, H, {'A'}, 'Kinds', struct('A', 'plot'), 'Keys', struct('A', 'g'));
    hFig = F.fig;  hFig.Visible = o.Visible;
    P = F.ax.A.Position;  axis(F.ax.A, 'off');
    gap = 2.5 / F.size_mm(2);  bottom = 8 / F.size_mm(2);  left = 26 / F.size_mm(1);   % horizontal row labels: they never reach a neighbouring row
    h = (P(4) - bottom - (n - 1) * gap) / n;
    for i = 1:n
        ax = axes(hFig, 'Position', [P(1) + left, P(2) + bottom + (n - i) * (h + gap), P(3) - left, h]);
        hold(ax, 'on');
        v = V(i);
        if ~isempty(v.plant), plot(ax, t, shown(v.plant, v.kind), '-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt); end
        if ~isempty(v.values)
            plot(ax, v.t, wrapped(v.values, v.kind), 'o', 'MarkerSize', T.marker.size_pt, 'MarkerFaceColor', T.rgb.trace_event, ...
                'MarkerEdgeColor', T.rgb.neutral_ink, 'LineWidth', T.marker.edge_pt);
        end
        ax.XLim = t([1 end]);
        switch v.kind
            case 'angle'
                ax.YLim = [-pi pi];  ax.YTick = [-pi 0 pi];  ax.YTickLabel = {'$-\pi$', '0', '$\pi$'};  ax.TickLabelInterpreter = 'latex';
            case 'sign'
                ax.YLim = [-1.4 1.4];  ax.YTick = [-1 1];  ax.YTickLabel = {[char(8722) '1'], '+1'};
        end
        xl = '';  if i == n, xl = 'time (s)'; end
        nxr_axes_style(ax, xl, sprintf('%s (%s)', v.label, v.unit));
        ax.YLabel.Rotation = 0;  ax.YLabel.HorizontalAlignment = 'right';  ax.YLabel.VerticalAlignment = 'middle';
        if i < n, ax.XTickLabel = {}; end
        hold(ax, 'off');
    end
    keys = {};
    if hasPlant, keys{end+1} = struct('kind', 'line', 'text', 'plant (analytic course)', 'color', T.rgb.neutral_ink); end
    if hasDet, keys{end+1} = struct('kind', 'marker', 'text', 'detected', 'color', T.rgb.trace_event); end
    nxr_key(F, 'A', keys);
    finish(F, o.Title, o.Export);
end

function y = wrapped(x, kind)
    y = x;
    if strcmp(kind, 'angle'), y = mod(x + pi, 2 * pi) - pi; end
end

function y = shown(x, kind)
% a line: an angle wrapped, broken where it wraps (no vertical jump across the axis)
    y = wrapped(x, kind);
    if strcmp(kind, 'angle'), y([false abs(diff(y)) > pi]) = NaN; end
end
