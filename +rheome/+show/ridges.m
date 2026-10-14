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
%   'Title'      figure title: the figure's name (the legend carries it, not the artwork)
%   'Export'     file stem: write <stem>.pdf/.svg/.png at print size (nxr_save)
%   'Visible'    'on' (default) | 'off'
%
% OUTPUT: hFig figure handle.
%
% See also: rheome.detect.ridges, rheome.show.surface
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('MaxRidges', 60);
    p.addParameter('View', [-90 10]);
    p.addParameter('Background', []);
    p.addParameter('Truth', struct());
    p.addParameter('Title', 'scale-space ridges');
    p.addParameter('Lift', 0.012);   % radial offset for trajectories, fraction of the mesh extent
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});
    o = p.Results;

    R = out.ridges;
    if isempty(R), error('show:ridges:empty', 'No ridges to show.'); end
    sig  = out.summary.sigmas;
    life = [R.lifetime];
    spd  = [R.speedNet];
    smod = [R.sigmaMode];

    [~, ord] = sort([R.nFrames], 'descend');
    ord = ord(1:min(o.MaxRidges, numel(ord)));
    lo = min(sig); hi = max(sig);  if hi <= lo, hi = lo + eps; end
    F = nxr_figure('double', 65, {'ABCD'}, 'Widths', [2 1 1 1], 'Kinds', struct('A', 'cortex'), 'Titles', struct( ...
        'A', sprintf('%d ridges, longest %d frames', numel(R), R(ord(1)).nFrames), ...
        'B', 'Selected scale', ...
        'C', 'Lifetime', 'D', 'Net speed'));      % medians are statistics: the legend's
    hFig = F.fig;  hFig.Visible = o.Visible;
    ref = {'--', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.leader_pt, ...
           'FontSize', T.type.annotation.size_pt};

    % ---- trajectories on the cortex, coloured by scale (the key is panel D's colourbar) ----
    ax1 = F.ax.A;
    if isempty(o.Background)                          % flat: the neutral cortex, no colormap used
        rheome.show.surface(surf, [], 'Parent', ax1, 'View', o.View, 'Title', '');
    else
        rheome.show.surface(surf, o.Background, 'Parent', ax1, 'View', o.View, 'Title', '', 'Colorbar', false);
    end
    hold(ax1,'on');
    if isfield(surf,'VertNormals') && ~isempty(surf.VertNormals)
        Nrm = surf.VertNormals ./ max(vecnorm(surf.VertNormals,2,2), eps);
    else
        Nrm = surf.Vertices ./ max(vecnorm(surf.Vertices,2,2), eps);
    end
    span = max(max(surf.Vertices, [], 1) - min(surf.Vertices, [], 1));
    lift = @(P, idx) P + o.Lift * span * Nrm(idx, :);
    cmap = nxr_cmap('atlas-tile-code');
    for k = flip(ord)                                   % longest drawn last, so it sits on top
        f = (R(k).sigmaMode - lo) / (hi - lo);
        col = cmap(max(1, min(256, round(1 + f*255))), :);
        P = lift(R(k).pos, R(k).v);
        lw = T.line.data_pt;  if k == ord(1), lw = 2*T.line.event_pt; end
        plot3(ax1, P(:,1), P(:,2), P(:,3), '-', 'LineWidth', lw, 'Color', col);
        plot3(ax1, P(1,1), P(1,2), P(1,3), 'o', 'Color', col, 'MarkerSize', T.marker.size_pt, 'MarkerFaceColor', col);
    end
    if isfield(o.Truth,'path') && ~isempty(o.Truth.path)
        Tp = o.Truth.path;
        Tp = Tp .* (1 + 1.4*o.Lift*span./max(vecnorm(Tp,2,2),eps));
        plot3(ax1, Tp(:,1), Tp(:,2), Tp(:,3), '--', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.event_pt);
    end

    % ---- scale histogram ----
    ax2 = F.ax.B;
    bar(ax2, 1000*sig, out.summary.scaleHist, 'FaceColor', T.rgb.neutral_slate, 'EdgeColor','none');
    if isfield(o.Truth,'sigma') && ~isempty(o.Truth.sigma)
        hold(ax2,'on');  xline(ax2, 1000*o.Truth.sigma, ref{1}, 'true σ', ref{2:end});
    end
    nxr_axes_style(ax2, 'σ (mm)', 'detections (count)');

    % ---- lifetime ----
    ax3 = F.ax.C;
    histogram(ax3, 1000*life, 'FaceColor', T.rgb.neutral_slate, 'EdgeColor','none');
    nxr_axes_style(ax3, 'lifetime (ms)', 'ridges (count)');

    % ---- speed vs lifetime, coloured by scale ----
    ax4 = F.ax.D;
    scatter(ax4, 1000*life, spd, T.marker.size_pt^2, 1000*smod, 'filled', 'MarkerFaceAlpha', 0.65);
    if isfield(o.Truth,'speed') && ~isempty(o.Truth.speed)
        hold(ax4,'on');  yline(ax4, o.Truth.speed, ref{1}, 'true speed', ref{2:end});
    end
    nxr_axes_style(ax4, 'lifetime (ms)', 'net speed (m/s)');
    nicebar(nxr_colorbar(ax4, 'ridge scale σ', 'mm', 'code', 1000*[lo hi]), ax4);
    finish(F, o.Title, o.Export);
end

% Author: Diellor Basha, 2026
