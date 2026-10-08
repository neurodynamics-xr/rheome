function drawCortex(app)
% DRAWCORTEX  Update the per-scale surfaces at the current frame.
%
% ⭐ THE FAST PATH. The patch objects are created once by buildCortexAxes; this only rewrites
% FaceVertexCData, so nothing re-renders and scrubbing stays interactive.
%
% ⚠ ONE COLOUR SCALE ACROSS SCALES, per frame. Giving each scale its own limits would make
% every scale look equally active and destroy exactly the comparison the tab exists for --
% which spatial scale currently carries the flow.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    nS = fp.NumScales;
    % ⭐ ONE synthesis per frame. Signed, envelope and dominant-scale all derive from the
    % same complex Z; computing them separately doubled the redraw cost.
    Z  = complexmaps(fp, app.Frame);                     % [V x (nS+1)] complex, one GEMM
    switch lower(app.MapMode)
        case 'signed',    M = double(real(Z));
        otherwise,        M = double(abs(Z));
    end

    switch lower(app.MapMode)
        case 'signed'
            lim = max(abs(M(:)));  if lim <= 0, lim = 1; end
            clim = [-lim lim];  cmap = rheome.flowbrowser.coolwarm(256);
        otherwise
            lim = max(M(:));       if lim <= 0, lim = 1; end
            clim = [0 lim];     cmap = parula(256);
    end

    % panel 1 = all scales; the dominant map takes panel 2 when enabled; then the scales
    lead = 1 + app.ShowDominant;
    for j = 0:nS
        if j == 0, g = 0; idx = 1; else, g = app.ScaleOrder(j); idx = j + lead; end
        h = app.CortexPatch(idx);
        if ~isvalid(h), continue; end
        set(h, 'FaceVertexCData', M(:, g+1));
        ax = ancestor(h, 'axes');
        set(ax, 'CLim', clim);  colormap(ax, cmap);
    end
    if app.ShowDominant, app.drawDominant(Z); end

    % ---- label each scale with its wavelength AND its share of the energy ----
    % ⚠ THE SHARE IS WHAT IDENTIFIES THE PEAK, not the colour scale. Every surface is drawn
    % on the SAME limits (so scales are comparable), which means a weak scale still fills its
    % axes with colour -- it just spans less of the range. Without a number the eye cannot
    % rank them, and "which scale carries the flow" is the question this tab exists for.
    kc = centerWavenumbers(fp.GraphBank);
    E  = scaleEnergy(fp, app.Frame);
    frac = E / max(sum(E), realmin);
    [~, gTop] = max(E);
    for j = 1:nS
        g  = app.ScaleOrder(j);
        mm = 1000 * 2*pi / kc(g);
        L  = app.CortexTitle(j + lead);
        if g == gTop
            L.Text = sprintf('● %.0f mm — %.0f%%', mm, 100*frac(g));
            L.FontColor = app.ScaleColors(g, :) * 0.75;
            L.FontWeight = 'bold';
        else
            L.Text = sprintf('%.0f mm — %.0f%%', mm, 100*frac(g));
            L.FontColor = [0.35 0.35 0.35];
            L.FontWeight = 'normal';
        end
    end
    if app.ShowDominant
        app.CortexTitle(2).Text = 'dominant scale  (warm coarse → cool fine)';
        app.CortexTitle(2).FontColor = [0 0 0];
        app.CortexTitle(2).FontWeight = 'bold';
    end
    kb = centroid(fp, app.Frame);
    app.CortexTitle(1).Text = sprintf('all scales   t = %.3f s   centroid %.0f mm', ...
        fp.Time(app.Frame), 1000*2*pi/kb);
    app.CortexTitle(1).FontColor = [0 0 0];
    app.CortexTitle(1).FontWeight = 'bold';
end

% Author: Diellor Basha, 2026
