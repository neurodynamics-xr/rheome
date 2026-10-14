function [hFig, info] = activity(surf, c, Psi, t, varargin)
% SHOW.ACTIVITY  Where and when vortical activity is strong.
%
%   [hFig, info] = rheome.show.activity(surf, c, Psi, t)
%
% Three panels: WHERE (time-averaged amplitude map on the cortex), WHEN (total amplitude over
% time), and the distribution of instantaneous amplitude. Use it to choose the regions and time
% windows worth looking at before committing to detection.
%
% The WHEN panel is free. Psi is M-orthonormal, so by Parseval the total vortical energy at time t
% is just the coefficient norm:  sum_v |Psi c(:,t)|^2  ==  ||c(:,t)||^2 -- no synthesis needed.
% Only the WHERE panel requires the vertex-space multiply.
%
% INPUTS:
%   surf  surface struct   c [K x nT] complex coefficients   Psi [N x K]   t [1 x nT] times (s)
% OPTIONS:
%   'Percentile'  mark windows above this percentile of amplitude (default 90)
%   'Unit'        unit of the amplitude (default 'a.u.')
%   'View','Visible'; 'Title' the figure's name (the legend carries it); 'Export' file stem
%   Sizes and colours are the figure tokens; 190 x 65 mm. The peak is ringed.
%
% OUTPUT:
%   info .whenPower [1 x nT] energy per frame (Parseval)   .peakT time of maximum
%        .windows  [n x 2] start/stop of the runs above the percentile
%        .whereMap [N x 1] time-averaged amplitude
%
% See also: rheome.filters.frame_scalogram, rheome.filters.joint_scalogram
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('Percentile', 90);
    p.addParameter('View', [-90 10]);
    p.addParameter('Unit', 'a.u.');
    p.addParameter('Title','');
    p.addParameter('Export','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    t  = double(t(:)).';
    pw = sum(abs(c).^2, 1);                 % Parseval: energy per frame, no synthesis
    A  = abs(Psi * c);                       % [N x nT] amplitude map over time
    mp = mean(A, 2);

    thr  = prctile(pw, o.Percentile);
    above = pw >= thr;
    d = diff([false above false]);
    runs = [t(d==1).', t(find(d==-1)-1).'];

    info.whenPower = pw;
    info.whereMap  = mp;
    [~,ip] = max(pw);  info.peakT = t(ip);
    info.windows = runs;

    F = nxr_figure('double', 65, {'ABC'}, 'Kinds', struct('A', 'cortex'), 'Titles', struct( ...
        'A', 'Where: time-averaged amplitude', ...
        'B', sprintf('When: %d window(s) above p%g', size(runs,1), o.Percentile), ...
        'C', 'Amplitude distribution'));
    hFig = F.fig;  hFig.Visible = o.Visible;

    ax1 = F.ax.A;
    rheome.show.surface(surf, mp, 'Parent', ax1, 'View', o.View, 'Colorbar', false, 'Title', '');
    nicebar(nxr_colorbar(ax1, 'mean |ω|', o.Unit, 'magnitude', [0 max(max(mp), eps)]), ax1);

    ax2 = F.ax.B;  hold(ax2, 'on');
    plot(ax2, t, pw, 'Color', T.rgb.trace_neutral, 'LineWidth', T.line.data_pt);
    y = pw;  y(~above) = NaN;                % the windows on the trace itself (§5c), no shaded spans
    plot(ax2, t, y, 'Color', T.rgb.trace_event, 'LineWidth', T.line.event_pt);
    yline(ax2, thr, '--', sprintf('p%g', o.Percentile), 'Color', T.rgb.neutral_slate, ...
        'LineWidth', T.line.axis_pt, 'FontSize', T.type.annotation.size_pt);
    plot(ax2, info.peakT, pw(ip), 'o', 'MarkerSize', 2*T.marker.size_pt, 'MarkerEdgeColor', T.rgb.neutral_ink, ...
        'LineWidth', T.marker.edge_pt);
    nxr_axes_style(ax2, 'time (s)', ['energy Σ|ω|² (' o.Unit '²)']);
    xlim(ax2, [t(1) t(end)]);

    ax3 = F.ax.C;
    histogram(ax3, A(:), 60, 'FaceColor', T.rgb.neutral_slate, 'EdgeColor', 'none');
    set(ax3, 'YScale', 'log');
    nxr_axes_style(ax3, ['|ω| per vertex and frame (' o.Unit ')'], 'count (vertex-frames)');
    plainlog(ax3);
    finish(F, o.Title, o.Export);
end

% Author: Diellor Basha, 2026
