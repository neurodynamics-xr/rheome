function drawSpectrum(app)
% DRAWSPECTRUM  The lambda spectrum at the current frame -- the GLOBAL/LOCAL readout.
%
% ⭐ Low lambda = globally distributed structure; high lambda = spatially confined. Plotted
% against sqrt(lambda), which is a WAVENUMBER, so the x axis is a physical inverse length
% and the top axis can carry it in mm.
%
% ⚠ imagesc/plot against lambda directly would compress the low end into invisibility --
% lambda spans four decades here. sqrt(lambda) is the natural axis and the one the wavelength
% is defined from.
%
% ⚠ SORT FIRST. rheome.flow.context assembles the basis block-diagonally per hemisphere, so mode
% index is NOT lambda order; plotting a line against the unsorted axis draws a zigzag that
% looks like structure and is pure ordering.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    S  = modeSpectrum(fp, app.Frame);
    k  = sqrt(fp.Lambda);
    [k, ord] = sort(k);  S = S(ord);

    ax = app.SpecAxMode;  cla(ax);
    semilogy(ax, k, max(S, eps), '-', 'Color', [0.25 0.35 0.7], 'LineWidth', 0.9);
    hold(ax, 'on');
    kb = centroid(fp, app.Frame);
    xline(ax, kb, '-', sprintf('%.0f mm', 2000*pi/kb), 'Color', [0.85 0.2 0.2], ...
        'LineWidth', 1.4, 'LabelVerticalAlignment', 'top');
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'wavenumber \surd\lambda  (rad/m)   \leftarrow global    local \rightarrow');
    ylabel(ax, 'mode energy');
    title(ax, sprintf('\\lambda spectrum   t = %.3f s', fp.Time(app.Frame)));

    ax2 = app.SpecAxScale;  cla(ax2);
    E  = scaleEnergy(fp, app.Frame);
    kc = centerWavenumbers(fp.GraphBank);
    bar(ax2, 1:numel(E), E, 0.7, 'FaceColor', [0.35 0.55 0.45]);
    ax2.XTick = 1:numel(E);
    ax2.XTickLabel = compose('%.0f', 1000*2*pi./kc(:));
    grid(ax2, 'on');
    xlabel(ax2, 'spatial scale (mm)');
    ylabel(ax2, 'energy');
    gi = globalIndex(fp);
    title(ax2, sprintf('per scale   global index %.2f', gi(app.Frame)));
end

% Author: Diellor Basha, 2026
