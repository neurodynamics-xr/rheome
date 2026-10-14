function hFig = jointscalogram(scal, varargin)
% SHOW.JOINTSCALOGRAM  Energy over spatial scale x temporal frequency, with the speed ridge.
%
%   rheome.show.jointscalogram(scal)                 % scal from rheome.filters.joint_scalogram
%   hFig = rheome.show.jointscalogram(scal, name, value, ...)
%
% Left: E(m,omega) -- which spatial scale carries the field at each temporal frequency.
% Right: the ridge, the energy-weighted mean wavenumber per frequency, with the fitted line
% omega = c*sqrt(lambda). The fitted speed and its R^2 are printed, together with .fSupport, the
% frequency interval carrying 90%% of the energy -- if the ridge does not span it, the slope is
% clipped and should not be believed.
%
% OPTIONS: 'Visible', 'Title' (the figure's name), 'Export' file stem (<stem>.pdf/.svg/.png),
%          'Theme' ('' = the tokens' default, Pergammon; 'neuroimage' for MS1).
% Sizes and colours are the figure tokens; double column x 58 mm; colourbars in the key bands.
%
% See also: rheome.filters.joint_scalogram, rheome.show.jointspectrum
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Title','');
    p.addParameter('Export','');
    p.addParameter('Visible','on');
    p.addParameter('Theme', '');
    p.parse(varargin{:});  o = p.Results;
    T = tokens(o.Theme);

    hasRidge = ~isempty(scal.kbar);
    tB = 'No ridge';
    if hasRidge
        tB = sprintf('Speed ridge %.3f m/s over %.1f–%.1f Hz', scal.slope, scal.fSupport(1), scal.fSupport(2));
    end
    F = nxr_figure('double', 58, {'AB'}, 'Theme', o.Theme, 'Keys', struct('A', 's', 'B', 's'), 'Titles', struct('A', 'Joint scalogram', 'B', tB));
    hFig = F.fig;  hFig.Visible = o.Visible;

    ax1 = F.ax.A;
    lab = arrayfun(@(s) ternary(isnan(s),'LP',sprintf('%.0f',1000*s)), scal.sigma, 'UniformOutput', false);
    Z = log10(max(scal.energy, max(scal.energy(:))*1e-8));
    imagesc(ax1, scal.f, 1:numel(scal.sigma), Z);
    axis(ax1,'xy');
    set(ax1,'YTick',1:numel(lab),'YTickLabel',lab);
    nxr_axes_style(ax1, 'temporal frequency (Hz)', 'member scale $\sigma$ (mm)');
    nicebar(nxr_colorbar(ax1, 'energy $\log_{10} E$', 'a.u.', 'power', [min(Z(:)) max(Z(:)) + eps], 'F', F), ax1);

    ax2 = F.ax.B;
    if hasRidge
        scatter(ax2, scal.kbar, 2*pi*scal.fRidge, T.marker.size_pt^2, scal.wRidge, 'filled'); hold(ax2,'on');
        kfit = linspace(min(scal.kbar), max(scal.kbar), 50);
        plot(ax2, kfit, scal.slope*kfit + (2*pi*mean(scal.fRidge) - scal.slope*mean(scal.kbar)), ...
            '-', 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt);
        text(ax2, 0.03, 0.97, sprintf('$\\omega = c\\sqrt{\\lambda}$, $R^2 = %.3f$', scal.slopeR2), 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', T.type.annotation.size_pt);
        nxr_axes_style(ax2, 'energy-weighted $\sqrt{\lambda}$ (rad/m)', 'angular frequency $\omega$ (rad/s)');
        nicebar(nxr_colorbar(ax2, 'energy weight', 'ratio', 'power', [0 max(max(scal.wRidge), eps)], 'F', F), ax2);
    else
        axis(ax2,'off');
    end
    finish(F, o.Title, o.Export);
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
