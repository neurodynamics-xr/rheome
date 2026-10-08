function hFig = specparam(P, f, ap, varargin)
% SHOW.SPECPARAM  Aperiodic fits and the periodic/aperiodic split, per spatial mode.
%
%   rheome.show.specparam(P, f, ap)
%   hFig = rheome.show.specparam(P, f, ap, name, value, ...)
%
% Left: a few modes' spectra with their fitted aperiodic background, log-log.
% Middle: the split gains h_per and h_ap for one mode, which square to one.
% Right: chi(lambda) -- the fitted exponent against spatial frequency. That panel is the point of
% fitting per mode rather than per channel: it asks whether the aperiodic slope depends on
% spatial scale, which a per-channel fit cannot see.
%
% INPUTS: P [nF x nS] linear power, f [nF x 1], ap from rheome.spectral.aperiodic
% OPTIONS: 'Modes' which columns to draw (default: 5 spread across), 'Lambda', 'Title', 'Visible'
%
% See also: rheome.spectral.aperiodic, rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Modes', []);
    p.addParameter('Lambda', []);
    p.addParameter('Title','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    f = double(f(:));  nS = size(P,2);
    if isempty(o.Modes), o.Modes = unique(round(linspace(1, nS, min(5,nS)))); end

    hFig = figure('Color','w','Position',[60 60 1300 400],'Visible',o.Visible);
    cmap = parula(numel(o.Modes)+1);

    ax1 = subplot(1,3,1); hold(ax1,'on');
    for i = 1:numel(o.Modes)
        k = o.Modes(i);
        plot(ax1, f, P(:,k), '-',  'Color', cmap(i,:), 'LineWidth', 1.0);
        plot(ax1, f, ap.Pap(:,k), '--', 'Color', cmap(i,:), 'LineWidth', 1.6);
    end
    set(ax1,'XScale','log','YScale','log'); grid(ax1,'on');
    xlabel(ax1,'frequency (Hz)'); ylabel(ax1,'power');
    title(ax1,'spectra (—) and aperiodic fits (--)','FontWeight','normal');

    ax2 = subplot(1,3,2);
    k = o.Modes(ceil(end/2));
    [hp, ha] = rheome.spectral.gains(P(:,k), ap.Pap(:,k));
    plot(ax2, f, hp, 'LineWidth', 1.6); hold(ax2,'on');
    plot(ax2, f, ha, 'LineWidth', 1.6);
    plot(ax2, f, hp.^2 + ha.^2, 'k:', 'LineWidth', 1.2);
    legend(ax2, {'h_{per}','h_{ap}','h^2 sum = 1'}, 'Location','east'); grid(ax2,'on');
    xlabel(ax2,'frequency (Hz)'); ylabel(ax2,'gain'); ylim(ax2,[0 1.1]);
    title(ax2, sprintf('split gains, mode %d', k), 'FontWeight','normal');

    ax3 = subplot(1,3,3);
    if isempty(o.Lambda), x = 1:nS; xl = 'mode index';
    else,                 x = sqrt(double(o.Lambda(:))).'; xl = '\surd\lambda (rad m^{-1})'; end
    plot(ax3, x, ap.exponent, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor','w'); grid(ax3,'on');
    xlabel(ax3, xl); ylabel(ax3,'aperiodic exponent \chi');
    title(ax3,'\chi(\lambda) — does the slope depend on scale?','FontWeight','normal');

    if ~isempty(o.Title), sgtitle(hFig, o.Title, 'FontWeight','bold','FontSize',13); end
end

% Author: Diellor Basha, 2026
