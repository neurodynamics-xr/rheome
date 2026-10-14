function hFig = frame(frameStruct, Lambda, varargin)
% SHOW.FRAME  A spectral-graph wavelet frame: member gains, coverage and scales.
%
%   rheome.show.frame(frame, Lambda)
%   hFig = rheome.show.frame(frame, Lambda, name, value, ...)
%
% Three panels that together say whether a filterbank is fit to use: the member gains g_m(lambda);
% the coverage curve sum_m g_m^2 with the frame bounds A and B marked (A -> 0 means part of the
% spectrum is uncovered, B/A = 1 means tight); and the member scales sigma_m = sqrt(2 t_m) against
% the wavenumbers the basis actually spans.
%
% This is the figure that replaces "we chose sigma = [10 15 20 30 40] mm" with a designed bank and
% a stated coverage.
%
% OPTIONS: 'Visible', 'Title' (the figure's name), 'Export' file stem (<stem>.pdf/.svg/.png),
%          'Theme' ('' = the tokens' default, Pergammon; 'neuroimage' for MS1).
% Sizes and colours are the figure tokens; double column x 60 mm; members coloured fine -> coarse.
%
% See also: rheome.filters.frame, rheome.filters.frame_bounds, rheome.filters.frame_gains
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Title', '');
    p.addParameter('Export', '');
    p.addParameter('Visible','on');
    p.addParameter('Theme', '');
    p.parse(varargin{:});  o = p.Results;
    T = tokens(o.Theme);

    lam = sort(double(Lambda(:)));
    H   = rheome.filters.frame_gains(frameStruct, lam);       % [K x M]
    b   = rheome.filters.frame_bounds(frameStruct, lam);
    M   = size(H,2);
    kk  = sqrt(lam);
    % EXACT member scale sqrt(2t), not the gain-weighted centroid -- the centroid compresses the
    % fine end and can give two distinct members the same label.
    if isfield(frameStruct,'Sigma') && ~isempty(frameStruct.Sigma)
        sig = frameStruct.Sigma;
    else
        sig = sqrt(2) ./ max(frameStruct.Centers, eps);
    end
    isScaling = strcmpi(frameStruct.Family,'mexhat');

    finest = 1000*pi/sqrt(max(lam));
    F = nxr_figure('double', 60, {'ABC'}, 'Theme', o.Theme, 'Titles', struct( ...
        'A', sprintf('%s frame, %d members', frameStruct.Family, M), ...
        'B', sprintf('Coverage, tightness %.3f', b.Tightness), ...
        'C', 'Member scales (1 = scaling function)'));
    hFig = F.fig;  hFig.Visible = o.Visible;
    ref = {'--', 'Color', T.rgb.neutral_slate, 'LineWidth', T.line.axis_pt, ...
           'FontSize', T.type.annotation.size_pt};

    ax1 = F.ax.A;  hold(ax1,'on');  cmap = ramp(M);
    for m = 1:M, plot(ax1, kk, H(:,m), 'LineWidth', T.line.data_pt, 'Color', cmap(m,:)); end
    nxr_axes_style(ax1, '$\sqrt{\lambda}$ (rad/m)', 'member gain $g_m$ (ratio)');

    ax2 = F.ax.B;
    plot(ax2, kk, sum(H.^2,2), 'Color', T.rgb.neutral_ink, 'LineWidth', T.line.data_pt); hold(ax2,'on');
    yline(ax2, b.A, ref{1}, sprintf('$A = %.3g$', b.A), ref{2:end});
    yline(ax2, b.B, ref{1}, sprintf('$B = %.3g$, $B/A = %.3f$', b.B, b.Tightness), ref{2:end});
    nxr_axes_style(ax2, '$\sqrt{\lambda}$ (rad/m)', 'coverage $\sum_m g_m^2$ (ratio)');

    ax3 = F.ax.C;
    idx = 1:M;  lab = 1000*sig;
    if isScaling, lab(1) = NaN; end
    stem(ax3, idx, lab, 'filled', 'LineWidth', T.line.data_pt, 'Color', T.rgb.neutral_ink, 'MarkerSize', T.marker.size_pt);
    set(ax3,'YScale','log');
    yline(ax3, finest, ref{1}, sprintf('basis limit %.1f mm', finest), ref{2:end});
    nxr_axes_style(ax3, 'member (index)', 'scale $\sigma$ (mm)');
    plainlog(ax3);
    finish(F, o.Title, o.Export);
end

% Author: Diellor Basha, 2026
