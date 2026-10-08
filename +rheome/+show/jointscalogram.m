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
% OPTIONS: 'Title', 'Visible'
%
% See also: rheome.filters.joint_scalogram, rheome.show.jointspectrum
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Title','');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

    hFig = figure('Color','w','Position',[60 60 1150 420],'Visible',o.Visible);

    ax1 = subplot(1,2,1);
    lab = arrayfun(@(s) ternary(isnan(s),'LP',sprintf('%.0f',1000*s)), scal.sigma, 'UniformOutput', false);
    imagesc(ax1, scal.f, 1:numel(scal.sigma), log10(max(scal.energy, max(scal.energy(:))*1e-8)));
    axis(ax1,'xy');
    set(ax1,'YTick',1:numel(lab),'YTickLabel',lab);
    xlabel(ax1,'temporal frequency (Hz)'); ylabel(ax1,'member \sigma_m (mm)');
    cb = colorbar(ax1); cb.Label.String = 'log_{10} energy';
    title(ax1,'joint scalogram','FontWeight','normal');

    ax2 = subplot(1,2,2);
    if ~isempty(scal.kbar)
        scatter(ax2, scal.kbar, 2*pi*scal.fRidge, 14, scal.wRidge, 'filled'); hold(ax2,'on');
        kfit = linspace(min(scal.kbar), max(scal.kbar), 50);
        plot(ax2, kfit, scal.slope*kfit + (2*pi*mean(scal.fRidge) - scal.slope*mean(scal.kbar)), ...
            'r-', 'LineWidth', 1.8);
        xlabel(ax2,'energy-weighted \surd\lambda (rad m^{-1})'); ylabel(ax2,'\omega (rad s^{-1})');
        grid(ax2,'on'); cb2 = colorbar(ax2); cb2.Label.String = 'energy weight';
        title(ax2, sprintf('c = %.3f m/s   (R^2 = %.3f)   support %.1f–%.1f Hz', ...
            scal.slope, scal.slopeR2, scal.fSupport(1), scal.fSupport(2)), 'FontWeight','normal');
    else
        text(ax2,0.5,0.5,'no ridge','HorizontalAlignment','center'); axis(ax2,'off');
    end
    if ~isempty(o.Title), sgtitle(hFig, o.Title, 'FontWeight','bold','FontSize',13); end
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
