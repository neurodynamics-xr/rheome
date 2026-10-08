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
% OPTIONS: 'Title', 'Visible'
%
% See also: rheome.filters.frame, rheome.filters.frame_bounds, rheome.filters.frame_gains
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Title', '');
    p.addParameter('Visible','on');
    p.parse(varargin{:});  o = p.Results;

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

    hFig = figure('Color','w','Position',[60 60 1300 400],'Visible',o.Visible);
    cmap = parula(M);

    ax1 = subplot(1,3,1); hold(ax1,'on');
    for m = 1:M, plot(ax1, kk, H(:,m), 'LineWidth', 1.5, 'Color', cmap(m,:)); end
    xlabel(ax1,'\surd\lambda (rad m^{-1})'); ylabel(ax1,'g_m(\lambda)'); grid(ax1,'on');
    title(ax1, sprintf('%s frame, %d members', frameStruct.Family, M), 'FontWeight','normal');

    ax2 = subplot(1,3,2);
    plot(ax2, kk, sum(H.^2,2), 'k', 'LineWidth', 1.8); hold(ax2,'on');
    yline(ax2, b.A, 'r--', sprintf('A = %.3g', b.A), 'LineWidth', 1.2);
    yline(ax2, b.B, 'b--', sprintf('B = %.3g', b.B), 'LineWidth', 1.2);
    xlabel(ax2,'\surd\lambda (rad m^{-1})'); ylabel(ax2,'\Sigma_m g_m^2'); grid(ax2,'on');
    title(ax2, sprintf('coverage — tightness B/A = %.3f', b.Tightness), 'FontWeight','normal');

    ax3 = subplot(1,3,3);
    idx = 1:M;  lab = 1000*sig;
    if isScaling, lab(1) = NaN; end
    stem(ax3, idx, lab, 'filled', 'LineWidth', 1.4, 'Color', [0.15 0.35 0.65]);
    xlabel(ax3,'member'); ylabel(ax3,'\sigma_m (mm)'); grid(ax3,'on');
    set(ax3,'YScale','log');
    finest = 1000*pi/sqrt(max(lam));
    yline(ax3, finest, 'r--', sprintf('basis limit %.1f mm', finest), 'LineWidth', 1.2);
    title(ax3,'member scales (member 1 = scaling function)','FontWeight','normal');

    if ~isempty(o.Title), sgtitle(hFig, o.Title, 'FontWeight','bold','FontSize',13); end
end

% Author: Diellor Basha, 2026
