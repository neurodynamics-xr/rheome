function drawTimeline(app)
% DRAWTIMELINE  Band power, lambda-centroid and the global index across the page.
%
% ⭐ THE HYPOTHESIS IS READ HERE. If recurrent patterns are coherent at high band power and
% fragment at desynchronisation, the centroid (upper trace, right axis) should rise -- move
% to finer scale -- as power falls. The two traces are drawn on one time axis so the
% relationship is visible directly rather than inferred from a fit.
%
% ⚠ THE MARGIN IS SHADED. Frames outside the page's core carry the wavelet cone and belong
% to the neighbouring pages too; summing statistics over them would double-count. They are
% shown, not hidden, so the boundary is visible while scrubbing.
%
% Author: Diellor Basha, 2026

    app.drawScalogram();

    fp = app.Page;
    t  = fp.Time;
    P  = bandpower_(fp);
    gi = globalIndex(fp);
    kb = arrayfun(@(i) centroid(fp, i), 1:fp.NumFrames);

    ax = app.TimeAxPower;  cla(ax);
    yyaxis(ax, 'left');
    plot(ax, t, 10*log10(max(P, realmin)), '-', 'LineWidth', 1.0);
    ylabel(ax, 'band power (dB)');
    yyaxis(ax, 'right');
    plot(ax, t, 1000*2*pi./kb, '-', 'LineWidth', 1.0);
    ylabel(ax, 'scale 2\pi/\bar{k}  (mm)');
    app.shadeMargin(ax);
    grid(ax, 'on');
    title(ax, 'power and spatial scale — does the scale fall as power falls?');

    ax2 = app.TimeAxScale;  cla(ax2);
    plot(ax2, t, gi, '-', 'LineWidth', 1.0, 'Color', [0.4 0.3 0.6]);
    ylim(ax2, [0 1]);
    app.shadeMargin(ax2);
    grid(ax2, 'on');
    xlabel(ax2, 'time (s)');
    ylabel(ax2, 'global index');
    title(ax2, 'fraction of energy in low spatial modes  (1 = global, 0 = local)');

    app.TimeCursors = [xline(app.ScalogramAx, t(app.Frame), 'r-', 'LineWidth', 1.2), ...
                       xline(ax,  t(app.Frame), 'r-', 'LineWidth', 1.2), ...
                       xline(ax2, t(app.Frame), 'r-', 'LineWidth', 1.2)];
end

% Author: Diellor Basha, 2026
