function h = zstack(ax, n, current)
% SHOW.ZSTACK  A z-stack indicator: n stacked planes (the dyadic levels), the current slice filled in ink.
%
%   h = rheome.show.zstack(ax, n, current)     current: 1 (finest, bottom) .. n (coarsest, top)
%
% Beside a cortical map or a scalogram slice it says WHICH level of the constant-Q stack is shown, without a colour
% (a level is the grid, not a field). The planes are drawn in oblique projection, outlined in slate; the current one
% filled in ink. ax is made equal-aspect and its axes hidden.
%
% See also: rheome.show.levelstack
%
% Author: Diellor Basha, 2026
    T = tokens();
    if ~(isscalar(n) && n >= 1 && n == round(n)) || ~(isscalar(current) && current >= 1 && current <= n)
        error('rheome:show:zstack', 'n planes (>= 1) and current in 1..n.');
    end
    hold(ax, 'on');
    sh = 0.45;  dy = 0.35;   % the oblique shear and the planes' spacing, in plane widths
    h = gobjects(n, 1);
    for k = 1:n
        y = (k - 1) * dy;
        x = [0 1 1 + sh sh];  yy = y + [0 0 0.3 0.3];
        if k == current
            h(k) = patch(ax, x, yy, T.rgb.neutral_ink, 'EdgeColor', T.rgb.neutral_ink, 'LineWidth', T.line.outline_pt);
        else
            h(k) = patch(ax, x, yy, T.rgb.neutral_paper, 'EdgeColor', T.rgb.neutral_slate, 'LineWidth', T.line.outline_pt, 'FaceAlpha', 0.9);
        end
    end
    axis(ax, 'equal', 'off');  hold(ax, 'off');
end
