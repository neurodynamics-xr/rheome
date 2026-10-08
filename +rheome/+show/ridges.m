function hFig = ridges(surf, out, varargin)
% SHOW.RIDGES  Summary figure for scale-space ridge tracking (rheome.detect.ridges).
%
%   rheome.show.ridges(surf, out)
%   hFig = rheome.show.ridges(surf, out, name, value, ...)
%
% Four panels that together say what a ridge run found: WHERE the motifs are and how they moved,
% at WHAT SIZE, for HOW LONG, and how fast. The scale axis is a first-class coordinate here --
% trajectories are coloured by their modal sigma, so size is readable off the surface directly.
%
% INPUTS:
%   surf  surface struct (.Vertices, .Faces, .nV) or a path
%   out   output of rheome.detect.ridges
%
% OPTIONS:
%   'MaxRidges'  how many trajectories to draw, longest first (default 60)
%   'View'       camera [az el] (default [-90 10])
%   'Background' [nV x 1] scalar to shade the surface (default flat grey)
%   'Truth'      struct for validation overlays, any of:
%                  .sigma  true scale (m)     -> marked on the scale histogram
%                  .speed  true speed (m/s)   -> marked on the speed panel
%                  .path   [n x 3] true core path -> drawn dashed black on the surface
%   'Lift'       trajectory offset along the normal, as a fraction of mesh extent (default 0.012)
%   'Title'      overall title
%   'Visible'    'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.detect.ridges, rheome.show.surface
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('MaxRidges', 60);
    p.addParameter('View', [-90 10]);
    p.addParameter('Background', []);
    p.addParameter('Truth', struct());
    p.addParameter('Title', 'scale-space ridges');
    p.addParameter('Lift', 0.012);   % radial offset for trajectories, fraction of the mesh extent
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;

    R = out.ridges;
    if isempty(R), error('show:ridges:empty', 'No ridges to show.'); end
    sig  = out.summary.sigmas;
    life = [R.lifetime];
    spd  = [R.speedNet];
    smod = [R.sigmaMode];

    hFig = figure('Color','w','Position',[60 60 1600 400],'Visible',o.Visible);

    % ---- panel 1: trajectories on the surface, coloured by modal scale (wider than the rest) ----
    ax1 = axes('Parent', hFig, 'Position', [0.02 0.08 0.30 0.80]);
    bg = o.Background;
    flatBg = isempty(bg);
    if flatBg, bg = zeros(surf.nV,1); end
    rheome.show.surface(surf, bg, 'Parent', ax1, 'View', o.View, 'Title', '');
    if flatBg
        cbs = findobj(hFig, 'Type', 'colorbar');   % a flat background has nothing to scale
        delete(cbs);
        colormap(ax1, repmat([0.82 0.82 0.84], 64, 1));
    end
    hold(ax1,'on');
    % Trajectories lie ON the mesh, so they z-fight with it and disappear behind folds.
    % Lift them a little along the outward normal (radially, if no normals are stored).
    if isfield(surf,'VertNormals') && ~isempty(surf.VertNormals)
        Nrm = surf.VertNormals ./ max(vecnorm(surf.VertNormals,2,2), eps);
    else
        Nrm = surf.Vertices ./ max(vecnorm(surf.Vertices,2,2), eps);
    end
    span = max(max(surf.Vertices, [], 1) - min(surf.Vertices, [], 1));
    lift = @(P, idx) P + o.Lift * span * Nrm(idx, :);

    [~, ord] = sort([R.nFrames], 'descend');
    ord = ord(1:min(o.MaxRidges, numel(ord)));
    cmap = parula(256);
    lo = min(sig); hi = max(sig);
    for k = flip(ord)                                   % longest drawn last, so it sits on top
        f = (R(k).sigmaMode - lo) / max(hi - lo, eps);
        col = cmap(max(1, min(256, round(1 + f*255))), :);
        P = lift(R(k).pos, R(k).v);
        isLongest = (k == ord(1));
        lw = 1.4;  if isLongest, lw = 3.5; end
        plot3(ax1, P(:,1), P(:,2), P(:,3), '-', 'LineWidth', lw, 'Color', col);
        plot3(ax1, P(1,1), P(1,2), P(1,3), 'o', 'Color', col, 'MarkerSize', 5, 'MarkerFaceColor', col);
    end
    if isfield(o.Truth,'path') && ~isempty(o.Truth.path)
        T = o.Truth.path;
        T = T .* (1 + 1.4*o.Lift*span./max(vecnorm(T,2,2),eps));
        plot3(ax1, T(:,1), T(:,2), T(:,3), 'k--', 'LineWidth', 2.5);
    end
    title(ax1, sprintf('%d ridges | longest %d frames (colour = \\sigma %.0f-%.0f mm)', ...
        numel(R), R(ord(1)).nFrames, 1000*lo, 1000*hi), 'FontWeight','normal');

    % ---- panel 2: selected-scale histogram ----
    ax2 = axes('Parent', hFig, 'Position', [0.385 0.15 0.17 0.70]);
    bar(ax2, 1000*sig, out.summary.scaleHist, 'FaceColor', [0.30 0.45 0.70], 'EdgeColor','none');
    xlabel(ax2, '\sigma (mm)'); ylabel(ax2, 'detections');
    title(ax2, 'selected scale', 'FontWeight','normal');
    if isfield(o.Truth,'sigma') && ~isempty(o.Truth.sigma)
        hold(ax2,'on');
        xline(ax2, 1000*o.Truth.sigma, 'r--', 'LineWidth', 2, 'Label', 'true \sigma');
    end
    grid(ax2,'on');

    % ---- panel 3: ridge lifetime ----
    ax3 = axes('Parent', hFig, 'Position', [0.615 0.15 0.17 0.70]);
    histogram(ax3, 1000*life, 'FaceColor', [0.35 0.60 0.40], 'EdgeColor','none');
    xlabel(ax3, 'lifetime (ms)'); ylabel(ax3, 'ridges');
    title(ax3, sprintf('lifetime (median %.1f ms, max %.0f)', 1000*median(life), 1000*max(life)), ...
        'FontWeight','normal');
    grid(ax3,'on');

    % ---- panel 4: net speed vs lifetime ----
    ax4 = axes('Parent', hFig, 'Position', [0.845 0.15 0.13 0.70]);
    scatter(ax4, 1000*life, spd, 18, 1000*smod, 'filled', 'MarkerFaceAlpha', 0.65);
    xlabel(ax4, 'lifetime (ms)'); ylabel(ax4, 'net speed (m/s)');
    cb = colorbar(ax4); cb.Label.String = '\sigma (mm)';
    title(ax4, sprintf('net speed (median %.3f m/s)', median(spd)), 'FontWeight','normal');
    grid(ax4,'on');
    if isfield(o.Truth,'speed') && ~isempty(o.Truth.speed)
        hold(ax4,'on');
        yline(ax4, o.Truth.speed, 'r--', 'LineWidth', 2, 'Label', 'true speed');
    end

    sgtitle(hFig, o.Title, 'FontWeight','bold');
end

% Author: Diellor Basha, 2026
