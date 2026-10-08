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
% OPTIONS: 'Frame' which frame to draw (default 1), 'View', 'Title', 'Visible'
%
% See also: rheome.flow.critical, rheome.detect.criticalPoints
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Frame', 1);
    p.addParameter('View', [-90 10]);
    p.addParameter('Title','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    col = struct('vortex',[0.85 0.15 0.15], 'source',[0.95 0.65 0.10], ...
                 'sink',[0.15 0.35 0.85], 'saddle',[0.35 0.35 0.35]);
    hFig = figure('Color','w','Position',[50 50 1330 420],'Visible',o.Visible);

    ax1 = subplot(1,3,1);
    rheome.show.surface(surf, zeros(size(surf.Vertices,1),1), 'Parent', ax1, 'View', o.View);
    hold(ax1,'on');
    c = out.cp(o.Frame);
    if ~isempty(c.pos)
        sz = 20 + 130 * (c.strength - min(c.strength)) / max(max(c.strength) - min(c.strength), eps);
        for i = 1:numel(c.type)
            scatter3(ax1, c.pos(i,1), c.pos(i,2), c.pos(i,3), sz(i), ...
                col.(lower(c.type{i})), 'filled', 'MarkerEdgeColor','k','LineWidth',0.4);
        end
    end
    title(ax1, sprintf('frame %d — %d critical points', out.frames(o.Frame), numel(c.type)), ...
        'FontWeight','normal');

    ax2 = subplot(1,3,2);
    bar(ax2, out.counts, 'stacked');
    set(ax2,'XTick',1:numel(out.frames),'XTickLabel',string(out.frames));
    legend(ax2, out.kinds, 'Location','eastoutside');
    xlabel(ax2,'frame'); ylabel(ax2,'count'); grid(ax2,'on');
    title(ax2,'population by type','FontWeight','normal');

    ax3 = subplot(1,3,3);
    if isempty(out.chi)
        text(ax3,0.5,0.5,'no \chi reported','HorizontalAlignment','center'); axis(ax3,'off');
    else
        plot(ax3, out.chi, 'o-','LineWidth',1.4,'MarkerFaceColor','w'); hold(ax3,'on');
        yline(ax3, out.chiExpected, 'r--', sprintf('\\chi = %d expected', out.chiExpected), 'LineWidth',1.4);
        xlabel(ax3,'frame'); ylabel(ax3,'\Sigma index'); grid(ax3,'on');
        ok = all(abs(out.chi(:) - out.chiExpected) < 0.5);
        title(ax3, sprintf('Poincaré–Hopf: %s', ternary(ok,'PASS','MISMATCH — check the mesh')), ...
            'FontWeight','normal');
    end
    if ~isempty(o.Title), sgtitle(hFig, o.Title, 'FontWeight','bold','FontSize',13); end
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
