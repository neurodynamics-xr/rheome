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
%   'Title'   overall title              'Visible' 'on' (default) | 'off'
%
% See also: rheome.flow.joint, rheome.filters.joint_scalogram, rheome.show.jointscalogram
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Speed', []);
    p.addParameter('Log', true);
    p.addParameter('Title', 'Joint spectrum of cortical vorticity');
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

    hFig = figure('Color','w','Position',[60 60 1240 560],'Visible',o.Visible);

    ax = subplot(4,4,[1 2 3 5 6 7 9 10 11]);
    imagesc(ax, f, 1:K, Z);  axis(ax,'xy');
    rheome.show.lambdaaxis(ax, lamv);
    ylabel(ax,'spatial frequency \surd\lambda (rad m^{-1})');
    set(ax,'XTickLabel',[]);                        % the marginal below carries the x label
    cb = colorbar(ax);  cb.Label.String = ternary(o.Log,'log_{10} |C|^2','|C|^2');
    title(ax, o.Title, 'FontWeight','bold');
    if ~isempty(o.Speed)
        hold(ax,'on');
        for c = o.Speed(:)'
            fc = c*ls/(2*pi);                        % omega = c*sqrt(lambda) -> f = c*k/2pi
            inR = fc >= min(f) & fc <= max(f);
            if any(inR)
                plot(ax, fc(inR), find(inR), 'w--', 'LineWidth', 1.4);
                j = find(inR, 1, 'last');
                text(ax, fc(j), j, sprintf(' %.2g m/s', c), 'Color','w', ...
                    'VerticalAlignment','bottom','FontSize',9);
            end
        end
        xlim(ax, [min(f) max(f)]);
    end

    ax2 = subplot(4,4,[13 14 15]);
    plot(ax2, f, sum(P,1), 'LineWidth', 1.2, 'Color', [0.15 0.35 0.65]);
    xlabel(ax2,'temporal frequency (Hz)'); ylabel(ax2,'power'); grid(ax2,'on');
    xlim(ax2,[min(f) max(f)]);

    ax3 = subplot(4,4,[4 8 12]);
    plot(ax3, sum(P,2), 1:K, 'LineWidth', 1.2, 'Color', [0.65 0.35 0.15]);
    rheome.show.lambdaaxis(ax3, lamv);
    xlabel(ax3,'power'); grid(ax3,'on');  ylim(ax3,[1 K]);
    title(ax3,'spatial marginal','FontWeight','normal');
end

function y = ternary(c,a,b), if c, y=a; else, y=b; end, end

% Author: Diellor Basha, 2026
