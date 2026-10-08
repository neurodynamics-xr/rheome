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
%   'View','Title','Visible'
%
% OUTPUT:
%   info .whenPower [1 x nT] energy per frame (Parseval)   .peakT time of maximum
%        .windows  [n x 2] start/stop of the runs above the percentile
%        .whereMap [N x 1] time-averaged amplitude
%
% See also: rheome.filters.frame_scalogram, rheome.filters.joint_scalogram
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Percentile', 90);
    p.addParameter('View', [-90 10]);
    p.addParameter('Title','');
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

    hFig = figure('Color','w','Position',[50 50 1320 430],'Visible',o.Visible);

    ax1 = subplot(1,3,1);
    rheome.show.surface(surf, mp, 'Parent', ax1, 'View', o.View);
    title(ax1, 'WHERE: time-averaged |\omega| amplitude', 'FontWeight','normal');

    ax2 = subplot(1,3,2);
    plot(ax2, t, pw, 'LineWidth', 1.1, 'Color', [0.15 0.35 0.65]); hold(ax2,'on');
    yline(ax2, thr, 'r--', sprintf('p%g', o.Percentile), 'LineWidth', 1.2);
    yl = ylim(ax2);
    for i = 1:size(runs,1)
        patch(ax2, [runs(i,1) runs(i,2) runs(i,2) runs(i,1)], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.9 0.3 0.2], 'FaceAlpha', 0.18, 'EdgeColor','none');
    end
    plot(ax2, info.peakT, pw(ip), 'r*', 'MarkerSize', 10);
    xlabel(ax2,'time (s)'); ylabel(ax2,'\Sigma_v |\omega|^2  (= ||c||^2)'); grid(ax2,'on');
    xlim(ax2,[t(1) t(end)]);
    title(ax2, sprintf('WHEN: %d window(s) above p%g; peak at %.2f s', size(runs,1), o.Percentile, info.peakT), ...
        'FontWeight','normal');

    ax3 = subplot(1,3,3);
    histogram(ax3, A(:), 60, 'FaceColor', [0.55 0.62 0.72], 'EdgeColor','none');
    set(ax3,'YScale','log'); xlabel(ax3,'|\omega| per vertex per frame'); ylabel(ax3,'count');
    grid(ax3,'on'); title(ax3,'amplitude distribution','FontWeight','normal');

    if ~isempty(o.Title), sgtitle(hFig, o.Title, 'FontWeight','bold','FontSize',13); end
end

% Author: Diellor Basha, 2026
