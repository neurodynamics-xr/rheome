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
% OPTIONS: 'Modes' which columns to draw (default: 5 spread across), 'Lambda', 'Unit' of P ('a.u.'),
%          'Visible', 'Title' (the figure's name), 'Export' file stem (<stem>.pdf/.svg/.png).
% Sizes and colours are the figure tokens; 190 x 60 mm; modes coloured low -> high.
%
% See also: rheome.spectral.aperiodic, rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('Modes', []);
    p.addParameter('Lambda', []);
    p.addParameter('Unit', 'a.u.');
    p.addParameter('Title','');
    p.addParameter('Export','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    f = double(f(:));  nS = size(P,2);
    if isempty(o.Modes), o.Modes = unique(round(linspace(1, nS, min(5,nS)))); end
    k = o.Modes(ceil(end/2));

    F = nxr_figure('double', 60, {'ABC'}, 'Titles', struct( ...
        'A', 'Spectra (solid) and aperiodic fits (dashed)', ...
        'B', sprintf('Split gains, mode %d', k), ...
        'C', 'Aperiodic exponent across scale'));
    hFig = F.fig;  hFig.Visible = o.Visible;
    cmap = ramp(numel(o.Modes));

    ax1 = F.ax.A;  hold(ax1,'on');
    for i = 1:numel(o.Modes)
        m = o.Modes(i);
        plot(ax1, f, P(:,m), '-',  'Color', cmap(i,:), 'LineWidth', T.line.data_pt);
        plot(ax1, f, ap.Pap(:,m), '--', 'Color', cmap(i,:), 'LineWidth', T.line.data_pt);
    end
    set(ax1,'XScale','log','YScale','log');
    nxr_axes_style(ax1, 'frequency (Hz)', ['power (' o.Unit ')']);
    plainlog(ax1);

    ax2 = F.ax.B;
    [hp, ha] = rheome.spectral.gains(P(:,k), ap.Pap(:,k));
    plot(ax2, f, hp, '-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt); hold(ax2,'on');
    plot(ax2, f, ha, '--', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt);
    plot(ax2, f, hp.^2 + ha.^2, ':', 'Color', T.rgb.neutral_slate, 'LineWidth', T.line.data_pt);
    ylim(ax2,[0 1.1]);
    nxr_axes_style(ax2, 'frequency (Hz)', 'gain (ratio)');
    legend(ax2, {'periodic h_per','aperiodic h_ap','h_per² + h_ap² = 1'}, 'Location','east', 'Box','off', ...
        'Interpreter','none', 'FontSize', T.type.legendText.size_pt);

    ax3 = F.ax.C;
    if isempty(o.Lambda), x = 1:nS; xl = 'mode (index)';
    else,                 x = sqrt(double(o.Lambda(:))).'; xl = '√λ (rad/m)'; end
    plot(ax3, x, ap.exponent, 'o-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt, ...
        'MarkerSize', T.marker.size_pt, 'MarkerFaceColor', T.rgb.neutral_paper);
    nxr_axes_style(ax3, xl, 'aperiodic exponent χ (ratio)');
    finish(F, o.Title, o.Export);
end

% Author: Diellor Basha, 2026
