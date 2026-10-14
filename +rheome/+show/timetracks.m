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
% how far one is from the other. A keyed plant's keys are ink diamonds. The key (plant / keys / detected) is in the
% panel's key band; every lane says its quantity and unit at its left (tokens dynamics.timeTrack, as
% style.time_tracks in Python).
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
%     Driver  the analytic driver, named once at the lane's end (default from Spec: 'ramp 0→6.28', '+ oscillation
%             2 cycles'; '' for none)
%
% OPTIONS (name/value):
%   'Width'   'single' (default) | 'onehalf' | 'double'
%   'RowHeight' mm per lane (default: the tokens' largest lane, dynamics.timeTrack.layout.lane_mm)
%   'Playhead'  times (s) of frames shown elsewhere in the figure: a thin slate line on every lane
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
    p.addParameter('RowHeight', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('Playhead', []);
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
    V = struct('label', {}, 'unit', {}, 'kind', {}, 'plant', {}, 't', {}, 'values', {}, 'driver', {}, 'keyT', {}, 'keyV', {});
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
        drv = '';  keysT = [];  keysV = [];
        if isfield(a, 'Spec') && ~isempty(a.Spec), [drv, keysT, keysV] = driverOf(a.Spec, t); end
        if isfield(a, 'Driver'), drv = char(a.Driver); end
        V(i) = struct('label', char(a.Label), 'unit', char(a.Unit), 'kind', kind, 'plant', plant, 't', tv, 'values', val, ...
            'driver', drv, 'keyT', keysT, 'keyV', keysV);
    end

    hasPlant = any(~cellfun(@isempty, {V.plant}));  hasDet = any(~cellfun(@isempty, {V.values}));
    hasKeys = any(~cellfun(@isempty, {V.keyT}));  hasDrv = any(~cellfun(@isempty, {V.driver}));
    lay = T.dynamics.timeTrack.layout;
    lane = o.RowHeight;  if isempty(lane), lane = max(lay.lane_mm); end
    nk = hasPlant + hasKeys + hasDet;  krows = max(1, ceil(nk / 2));
    H = n * lane + (n - 1) * lay.gap_mm + 24 + 5 * krows;
    F = nxr_figure(o.Width, H, {'A'}, 'Kinds', struct('A', 'plot'), 'Keys', struct('A', repmat('g', 1, krows)));
    hFig = F.fig;  hFig.Visible = o.Visible;
    P = F.ax.A.Position;  axis(F.ax.A, 'off');
    gap = lay.gap_mm / F.size_mm(2);  bottom = 8 / F.size_mm(2);  left = 26 / F.size_mm(1);   % horizontal lane labels: they never reach a neighbouring lane
    right = (8 + 22 * hasDrv) / F.size_mm(1);                                                  % the end values; the drivers' names
    h = lane / F.size_mm(2);                                    % the token's lane, laid from the top of the panel
    if n * h + (n - 1) * gap + bottom > P(4) + 1e-9, error('rheome:show:timetracks', 'the lanes do not fit the panel (RowHeight).'); end
    top = P(2) + P(4);
    for i = 1:n
        ax = axes(hFig, 'Position', [P(1) + left, top - i * h - (i - 1) * gap, P(3) - left - right, h]);
        hold(ax, 'on');
        v = V(i);
        if ~isempty(v.plant), plot(ax, t, shown(v.plant, v.kind), '-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt); end
        if ~isempty(v.keyT)
            plot(ax, v.keyT, wrapped(v.keyV, v.kind), 'd', 'MarkerSize', T.marker.size_pt, 'MarkerFaceColor', T.rgb.neutral_ink, ...
                'MarkerEdgeColor', T.rgb.neutral_paper, 'LineWidth', T.marker.edge_pt);
        end
        if ~isempty(v.values)
            plot(ax, v.t, wrapped(v.values, v.kind), 'o', 'MarkerSize', 0.7 * T.marker.size_pt, 'MarkerFaceColor', T.rgb.atlas_event, ...
                'MarkerEdgeColor', 'none');
        end
        for ph = double(o.Playhead(:)')
            xline(ax, ph, '-', 'Color', T.rgb.neutral_slate, 'LineWidth', T.line.outline_pt);
        end
        ax.XLim = t([1 end]);
        if strcmp(v.kind, 'value')   % the lane's two ends ticked: lanes are short
            y = [v.plant(:); v.values(:)];  y = y(isfinite(y));  lo = min(y);  hi = max(y);
            if hi <= lo, hi = lo + 1; lo = lo - 1; end
            ax.YLim = [lo hi] + [-0.08 0.08] * (hi - lo);  ax.YTick = [lo hi];  ax.YTickLabel = compose('%.3g', [lo hi]);
        end
        if ~isempty(v.driver)   % right of the end values
            text(ax, 1 + 9 / (F.size_mm(1) * ax.Position(3)), 0.5, v.driver, 'Units', 'normalized', 'FontSize', T.type.annotation.size_pt, 'Color', T.rgb.neutral_slate, ...
                'VerticalAlignment', 'middle', 'Interpreter', 'none');
        end
        switch v.kind
            case 'angle'
                ax.YLim = [-pi pi] * 1.08;  ax.YTick = [-pi pi];  ax.YTickLabel = {'$-\pi$', '$\pi$'};  ax.TickLabelInterpreter = 'latex';
            case 'sign'
                ax.YLim = [-1.4 1.4];  ax.YTick = [-1 1];  ax.YTickLabel = {[char(8722) '1'], '+1'};
        end
        xl = '';  if i == n, xl = 'time (s)'; end
        nxr_axes_style(ax, xl);
        % as style.time_tracks: the lane's end values right of it, each aligned inward (short lanes never collide), its
        % name and unit left of it
        tk = ax.YTick;  tl = ax.YTickLabel;  ip = ax.TickLabelInterpreter;
        ax.YAxisLocation = 'right';  ax.YTickLabel = {};
        al = {'bottom', 'top'};
        for q = 1:numel(tk)
            text(ax, 1 + 1.6 / (F.size_mm(1) * ax.Position(3)), (tk(q) - ax.YLim(1)) / diff(ax.YLim), tl{q}, 'Units', 'normalized', 'Interpreter', ip, ...
                'FontSize', T.type.tickLabel.size_pt, 'Color', T.rgb.neutral_ink, 'VerticalAlignment', al{min(q, 2)});
        end
        text(ax, -0.012, 0.5, sprintf('%s (%s)', v.label, v.unit), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'FontSize', T.type.axisTitle.size_pt, 'Color', T.rgb.neutral_ink, 'Interpreter', 'none');
        if i < n, ax.XTickLabel = {}; end
        hold(ax, 'off');
    end
    keys = {};
    if hasPlant, keys{end+1} = struct('kind', 'line', 'text', 'plant (analytic course)', 'color', T.rgb.neutral_ink); end
    if hasKeys, keys{end+1} = struct('kind', 'marker', 'text', sprintf('keyframes (%s, as the Editor''s)', char(9670)), 'color', T.rgb.neutral_ink); end
    if hasDet, keys{end+1} = struct('kind', 'marker', 'text', 'detected', 'color', T.rgb.atlas_event); end
    nxr_key(F, 'A', keys);
    finish(F, o.Title, o.Export);
end

function [d, kt, kv] = driverOf(spec, t)
% the driver's name at the lane's end, and a keyed track's keys in time
    d = '';  kt = [];  kv = [];
    if isfield(spec, 'Keys') && ~isempty(spec.Keys)
        K = spec.Keys;  if isnumeric(K), K = struct('x', num2cell(K(:, 1)), 'v', num2cell(K(:, 2))); end
        kt = t(1) + [K.x] * (t(end) - t(1));  kv = [K.v];
    elseif isfield(spec, 'To') && ~isempty(spec.To)
        s0 = 0;  if isfield(spec, 'Start') && ~isempty(spec.Start), s0 = spec.Start; end
        pr = 'ramp';  if isfield(spec, 'Profile') && ~isempty(spec.Profile), pr = char(spec.Profile); end
        d = sprintf('%s %.3g%s%.3g', pr, s0, char(8594), spec.To);
    end
    if isfield(spec, 'Osc') && ~isempty(spec.Osc)
        d = strtrim(sprintf('%s + oscillation %g cycles', d, spec.Osc.rate));
        if startsWith(d, '+ '), d = d(3:end); end
    end
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
