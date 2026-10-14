function hFig = critical(surf, out, varargin)
% SHOW.CRITICAL  Classified critical points on the cortex, with the Poincare-Hopf check.
%
%   rheome.show.critical(surf, out)          % out from rheome.flow.critical
%
% Left: the points of one frame on the surface, coloured by type and sized by strength.
% Middle: counts per type across frames -- how the population changes in time.
% Right: the index sum per hemisphere against the Euler characteristic. That panel is a TEST of
% the mesh and the operator, not a finding: sum(index) must equal chi on a closed manifold
% surface, so a departure means mesh trouble rather than anything about the data.
%
% OPTIONS: 'Frame' which frame to draw (default 1), 'View', 'Visible', 'Title' (the figure's name),
%          'Export' file stem (<stem>.pdf/.svg/.png). Sizes and colours are the figure tokens; 190 x 65 mm.
% Vortex = flow-vortex circles, source = flow-source triangles, sink = open flow-source triangles
% (pointing down), saddle = slate squares.
%
% See also: rheome.flow.critical, rheome.detect.criticalPoints
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('Frame', 1);
    p.addParameter('View', [-90 10]);
    p.addParameter('Title','');
    p.addParameter('Export','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    % Helmholtz colours of the tokens; shape doubles colour (source/sink share a colour, not a shape)
    sty = struct('vortex', {{'o', T.rgb.flow_vortex, T.rgb.flow_vortex}}, ...
                 'source', {{'^', T.rgb.flow_source, T.rgb.flow_source}}, ...
                 'sink',   {{'v', T.rgb.neutral_paper, T.rgb.flow_source}}, ...
                 'saddle', {{'s', T.rgb.neutral_slate, T.rgb.neutral_ink}});
    c = out.cp(o.Frame);
    ok = ~isempty(out.chi) && all(abs(out.chi(:) - out.chiExpected) < 0.5);
    F = nxr_figure('double', 65, {'ABC'}, 'Kinds', struct('A', 'cortex'), 'Titles', struct( ...
        'A', sprintf('Frame %d: %d critical points', out.frames(o.Frame), numel(c.type)), ...
        'B', 'Population by type', ...
        'C', sprintf('Poincaré–Hopf: %s', ternary(ok, 'pass', 'mismatch, check the mesh'))));
    hFig = F.fig;  hFig.Visible = o.Visible;

    ax1 = F.ax.A;
    rheome.show.surface(surf, [], 'Parent', ax1, 'View', o.View, 'Title', '');
    hold(ax1,'on');
    if ~isempty(c.pos)
        ms = T.marker.size_pt;
        sz = (ms * (1 + 2 * (c.strength - min(c.strength)) / max(max(c.strength) - min(c.strength), eps))).^2;
        for i = 1:numel(c.type)
            s = sty.(lower(c.type{i}));
            scatter3(ax1, c.pos(i,1), c.pos(i,2), c.pos(i,3), sz(i), s{1}, 'MarkerFaceColor', s{2}, ...
                'MarkerEdgeColor', s{3}, 'LineWidth', T.marker.edge_pt);
        end
    end

    ax2 = F.ax.B;
    hb = bar(ax2, out.counts, 'stacked', 'LineWidth', T.line.outline_pt);
    for i = 1:numel(hb)
        s = sty.(lower(out.kinds{i}));  hb(i).FaceColor = s{2};  hb(i).EdgeColor = s{3};
    end
    set(ax2,'XTick',1:numel(out.frames),'XTickLabel',string(out.frames));
    ylim(ax2, [0 1.7*max([sum(out.counts, 2); 1])]);   % empty space top-right for the key (tokens: legend)
    nxr_axes_style(ax2, 'frame (index)', 'critical points (count)');
    legend(ax2, out.kinds, 'Location','northeast', 'Box', 'off', 'FontSize', T.type.legendText.size_pt);

    ax3 = F.ax.C;
    if isempty(out.chi)
        text(ax3,0.5,0.5,'no χ reported','HorizontalAlignment','center'); axis(ax3,'off');
    else
        plot(ax3, out.chi, 'o-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt, ...
            'MarkerSize', T.marker.size_pt, 'MarkerFaceColor', T.rgb.neutral_paper); hold(ax3,'on');
        yline(ax3, out.chiExpected, '--', sprintf('χ = %d expected', out.chiExpected), 'Color', T.rgb.neutral_slate, ...
            'LineWidth', T.line.axis_pt, 'FontSize', T.type.annotation.size_pt);
        nxr_axes_style(ax3, 'frame (index)', 'Σ index (count)');
    end
    finish(F, o.Title, o.Export);
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
