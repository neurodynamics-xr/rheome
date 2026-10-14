function g = phasor(ax, phase, varargin)
% SHOW.PHASOR  Phase as a phasor glyph: a unit circle with the phase's vector(s), coloured on the cyclic phase map.
%
%   g = rheome.show.phasor(ax, phase)                         one glyph: the phases as points on the circle and their
%                                                              mean resultant vector (a polar plot of phase)
%   g = rheome.show.phasor(ax, phase, 'Amplitude', a)         the phases as vectors of length a (a phasor diagram)
%   g = rheome.show.phasor(ax, phase, 'At', xy, 'Radius', r)  one glyph per row of xy, each with its phase(k) (and
%                                                              a(k)): a phasor MARKER on a map, a scalogram or a trace
%
% Phase is an angle on a circle, not a number on a line: a glyph shows it as a direction (0 to the right,
% counter-clockwise positive, mathematical convention) and its concentration as the length of the mean resultant
% vector, R = |mean(a e^{i phase})| / mean(a) (1: all phases equal, 0: uniform or cancelling), as the sleep study's polar
% plots do. The colour of each phase is the figure tokens' cyclic phase map (quantity "phase", flow-phase), the map
% a phase colourbar uses, so glyphs and maps of phase read alike. Everything drawn needs a key (the legend says: the
% circle is |phasor| = Radius, the dark arrow the mean resultant).
%
% INPUTS:
%   ax     axes to draw in ([]: a new single-column figure, equal aspect, axes off)
%   phase  [n x 1] rad
%
% OPTIONS (name/value):
%   'Amplitude'  [n x 1] >= 0: vectors of length a / max(a) * Radius (default: none, points on the circle)
%   'At'         [1 x 2] the glyph's centre, or [n x 2] one glyph per phase (default [0 0])
%   'Radius'     the circle's radius in ax's data units (default 1)
%   'Mean'       draw the mean resultant vector (default true for one glyph, false for markers)
%   'Labels'     write 0, pi/2, pi, -pi/2 at the circle (default: true for one glyph without ax, else false)
%   'Theme'      '' (the tokens' default, Pergammon) | 'neuroimage' (MS1) | 'nature': colours, and the new figure's size
%
% OUTPUT:
%   g  struct: .R (mean resultant length) .mu (its angle, rad) per glyph, .handles
%
% Example:
%   rheome.show.phasor([], pi/3 + 0.4 * randn(60, 1));
%
% See also: rheome.show.timetracks
%
% Author: Diellor Basha, 2026
    p = inputParser;
    p.addParameter('Amplitude', []);
    p.addParameter('At', [0 0]);
    p.addParameter('Radius', 1, @(x) isscalar(x) && x > 0);
    p.addParameter('Mean', []);
    p.addParameter('Labels', []);
    p.addParameter('Theme', '');
    p.parse(varargin{:});
    o = p.Results;
    T = tokens(o.Theme);
    phase = double(phase(:));
    if isempty(phase) || any(~isfinite(phase)), error('rheome:show:phasor', 'phase: finite angles (rad), at least one.'); end
    a = double(o.Amplitude(:));
    if ~isempty(a) && (numel(a) ~= numel(phase) || any(a < 0 | ~isfinite(a))), error('rheome:show:phasor', 'Amplitude: one finite value >= 0 per phase.'); end
    At = double(o.At);
    markers = size(At, 1) > 1;
    if markers && size(At, 1) ~= numel(phase), error('rheome:show:phasor', 'At: one row per phase, or one centre.'); end
    newAxes = isempty(ax);
    if newAxes
        F = nxr_figure('single', 70, {'A'}, 'Kinds', struct('A', 'plot'), 'Letters', false, 'Theme', o.Theme);
        ax = F.ax.A;  axis(ax, 'equal', 'off');
    end
    useMean = o.Mean;  if isempty(useMean), useMean = ~markers; end
    labels = o.Labels;  if isempty(labels), labels = newAxes && ~markers; end
    cmap = nxr_cmap(T.quantity.phase.colormap);
    col = @(ph) cmap(1 + round((mod(ph + pi, 2*pi)) / (2*pi) * (size(cmap, 1) - 1)), :);
    held = ishold(ax);  hold(ax, 'on');
    r = o.Radius;  th = linspace(0, 2*pi, 97);
    if markers, groups = num2cell((1:numel(phase))'); else, groups = {(1:numel(phase))'}; At = At(1, :); end
    g = struct('R', zeros(numel(groups), 1), 'mu', zeros(numel(groups), 1), 'handles', gobjects(0));
    for k = 1:numel(groups)
        idx = groups{k};  c = At(k, :);  ph = phase(idx);
        w = ones(size(ph));  if ~isempty(a), w = a(idx); end
        z = sum(w .* exp(1i * ph)) / max(sum(w), eps);
        g.R(k) = abs(z);  g.mu(k) = angle(z);
        hs = plot(ax, c(1) + r * cos(th), c(2) + r * sin(th), '-', 'Color', T.rgb.neutral_slate, 'LineWidth', T.line.outline_pt);
        if ~markers   % quarter ticks
            for q = 0:3
                e = [cos(q*pi/2) sin(q*pi/2)];
                hs(end+1) = plot(ax, c(1) + r * [0.92 1.08] * e(1), c(2) + r * [0.92 1.08] * e(2), '-', 'Color', T.rgb.neutral_slate, 'LineWidth', T.line.tick_pt); %#ok<AGROW>
            end
        end
        if isempty(a)
            for j = 1:numel(ph)
                hs(end+1) = plot(ax, c(1) + r * cos(ph(j)), c(2) + r * sin(ph(j)), 'o', 'MarkerSize', T.marker.size_pt, ...
                    'MarkerFaceColor', col(ph(j)), 'MarkerEdgeColor', T.rgb.neutral_ink, 'LineWidth', T.marker.edge_pt); %#ok<AGROW>
            end
        else
            L = w / max(a);   % lengths relative to the largest amplitude shown
            for j = 1:numel(ph)
                hs = [hs arrow(ax, c, r * L(j), ph(j), col(ph(j)), T.line.arrow_pt, r)]; %#ok<AGROW>
            end
        end
        if useMean, hs = [hs arrow(ax, c, r * g.R(k), g.mu(k), T.rgb.neutral_ink, 1.6 * T.line.arrow_pt, r)]; end %#ok<AGROW>
        if labels
            lab = {'0', '$\pi/2$', '$\pi$', '$-\pi/2$'};  al = {'left', 'center', 'right', 'center'};  vl = {'middle', 'bottom', 'middle', 'top'};
            for q = 0:3
                hs(end+1) = text(ax, c(1) + 1.14 * r * cos(q*pi/2), c(2) + 1.14 * r * sin(q*pi/2), lab{q+1}, 'Interpreter', 'latex', ...
                    'FontSize', T.type.tickLabel.size_pt, 'Color', T.rgb.neutral_ink, 'HorizontalAlignment', al{q+1}, 'VerticalAlignment', vl{q+1}); %#ok<AGROW>
            end
        end
        g.handles = [g.handles hs];
    end
    if ~held, hold(ax, 'off'); end
end

function h = arrow(ax, c, len, ang, colour, lw, r)
% a shaft and a filled head; the head's size is a fraction of the glyph's radius, not of the vector
    e = [cos(ang) sin(ang)];  n = [-e(2) e(1)];
    hl = min(0.22 * r, 0.6 * len);
    tip = c + len * e;  base = tip - hl * e;
    h = plot(ax, [c(1) base(1)], [c(2) base(2)], '-', 'Color', colour, 'LineWidth', lw);
    if len > 0
        v = [tip; base + 0.45 * hl * n; base - 0.45 * hl * n];
        h(2) = patch(ax, v(:, 1), v(:, 2), colour, 'EdgeColor', 'none');
    end
end
