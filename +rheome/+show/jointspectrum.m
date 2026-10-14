function hFig = jointspectrum(C, f, Lambda, varargin)
% SHOW.JOINTSPECTRUM  The joint spectrum |C(lambda,omega)|^2 -- spatial vs temporal frequency.
%
%   rheome.show.jointspectrum(C, f, Lambda)
%   hFig = rheome.show.jointspectrum(C, f, Lambda, name, value, ...)
%
% Renders the endpoint of flow mapping: rows are spatial frequency sqrt(lambda), columns temporal
% frequency. A travelling structure satisfies omega = c*sqrt(lambda) and therefore appears as a
% DIAGONAL; a standing one as a horizontal band. Marginals are drawn alongside so the spatial and
% temporal spectra can be read off the same figure.
%
% OPTIONS:
%   'Speed'   overlay omega = c*sqrt(lambda) lines for these speeds (m/s), e.g. [0.1 0.5 1]
%   'Log'     true (default) -- log10 colour scale
%   'Title'   the joint spectrum's panel title    'Unit' of |C|^2 ('a.u.')
%   'Visible' 'on' (default) | 'off'          'Export' file stem (<stem>.pdf/.svg/.png)
% Sizes and colours are the figure tokens; 190 x 110 mm: spatial marginal in the space colour,
% temporal marginal in the time colour (atlas-space / atlas-time).
%
% See also: rheome.flow.joint, rheome.filters.joint_scalogram, rheome.show.jointscalogram
%
% Author: Diellor Basha, 2026

    T = tokens();
    p = inputParser;
    p.addParameter('Speed', []);
    p.addParameter('Log', true);
    p.addParameter('Title', 'Joint spectrum of cortical vorticity');
    p.addParameter('Unit', 'a.u.');
    p.addParameter('Export', '');
    p.addParameter('Visible', 'on');
    p.parse(varargin{:});  o = p.Results;

    f  = double(f(:)).';
    lamv = double(Lambda(:));
    % Sort by spatial frequency. The basis is assembled PER HEMISPHERE, so Lambda restarts partway
    % through and is not monotonic; plotting it directly as a y coordinate draws rows in the wrong
    % places and puts a spurious discontinuity at the hemisphere boundary.
    [ls, ord] = sort(sqrt(lamv));
    P  = abs(C(ord,:)).^2;
    K  = numel(ls);
    Z  = P;  if o.Log, Z = log10(max(P, max(P(:))*1e-8)); end

    F = nxr_figure('double', 110, {'AAAB'; 'AAAB'; 'CCC.'}, 'Titles', struct( ...
        'A', o.Title, 'B', 'Spatial marginal', 'C', 'Temporal marginal'));
    hFig = F.fig;  hFig.Visible = o.Visible;

    ax = F.ax.A;
    imagesc(ax, f, 1:K, Z);  axis(ax,'xy');
    colormap(ax, nxr_cmap('atlas-power'));
    nxr_axes_style(ax, '', 'spatial frequency √λ (rad/m)');
    rheome.show.lambdaaxis(ax, lamv);
    set(ax,'XTickLabel',[]);                        % the marginal below carries the x label
    cb = nxr_colorbar(ax, 'power |C|²', ternary(o.Log, ['log10 ' o.Unit], o.Unit), 'power', ...
        [min(Z(:)) max(Z(:)) + eps], '', 'keep');
    nicebar(cb, ax);
    rb = F.rect_mm.B;  rc = F.rect_mm.C;  W = F.size_mm(1);  H = F.size_mm(2);   % the bar in the free cell
    L = T.spacing.colorbarLength_frac * rb(3);
    cb.Position = [(rb(1) + (rb(3) - L)/2)/W, 1 - (rc(2) + rc(4)/2)/H, L/W, T.spacing.colorbarThickness_mm/H];
    if ~isempty(o.Speed)
        hold(ax,'on');
        for c = o.Speed(:)'
            fc = c*ls/(2*pi);                        % omega = c*sqrt(lambda) -> f = c*k/2pi
            inR = fc >= min(f) & fc <= max(f);
            if any(inR)
                plot(ax, fc(inR), find(inR), '--', 'Color', T.rgb.neutral_paper, 'LineWidth', T.line.leader_pt);
                j = find(inR, 1, 'last');
                text(ax, fc(j), j, sprintf('%.2g m/s ', c), 'Color', T.rgb.neutral_paper, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            end
        end
        xlim(ax, [min(f) max(f)]);
    end

    ax2 = F.ax.C;
    plot(ax2, f, sum(P,1), 'LineWidth', T.line.data_pt, 'Color', T.rgb.atlas_time);
    nxr_axes_style(ax2, 'temporal frequency (Hz)', ['power (' o.Unit ')']);
    xlim(ax2,[min(f) max(f)]);

    ax3 = F.ax.B;
    plot(ax3, sum(P,2), 1:K, 'LineWidth', T.line.data_pt, 'Color', T.rgb.atlas_space);
    rheome.show.lambdaaxis(ax3, lamv);
    nxr_axes_style(ax3, ['power (' o.Unit ')']);  ylim(ax3,[1 K]);  set(ax3, 'YTickLabel', []);
    finish(F, '', o.Export);
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
