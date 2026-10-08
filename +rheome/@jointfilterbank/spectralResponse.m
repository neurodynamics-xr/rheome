function [W, k, f] = spectralResponse(obj, m)
% SPECTRALRESPONSE  One member as an image on the (wavenumber, frequency) plane.
%
%   [W, k, f] = spectralResponse(jfb, m)
%   spectralResponse(jfb, m)                 plots
%
% Not called freqz: one axis is a WAVENUMBER, and freqz's 'z' is from the z-transform,
% which has no analogue here.
%
% A constant-speed structure lies on the DIAGONAL omega = c*sqrt(lambda), so a
% speed-selective member appears as a diagonal band in this picture.
%
% ⚠ THE WAVENUMBER AXIS IS PLOTTED BY ROW INDEX, NOT BY VALUE. sqrt(lambda) is strongly
% non-uniformly spaced, and imagesc assumes a UNIFORM y vector -- given one that is not,
% it maps row index linearly onto [min k, max k] and silently misplaces every row.
% Measured on an ico4/K=100 sphere: the degree-3 harmonic sits at k = 34.6 rad/m but was
% drawn at 8.8, an error of 25.9 rad/m. Rows are therefore indices and the tick LABELS
% carry the wavenumbers. To overlay a curve, map k to a row with
% interp1(k, 1:numel(k), kWanted).
%
% See also: dispersion, jointfilters
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(m), m = 1; end
    W = jointfilters(obj, m);
    k = sqrt(obj.Lambda);
    f = obj.Frequencies;
    if nargout == 0
        imagesc(f, 1:numel(k), abs(W));
        set(gca, 'YDir', 'normal');
        yt = round(linspace(1, numel(k), min(6, numel(k))));
        set(gca, 'YTick', yt, 'YTickLabel', compose('%.0f', k(yt)));
        xlabel('temporal frequency  f  (Hz)');
        ylabel('wavenumber  k = \surd\lambda  (rad/m)');
        title(sprintf('member %d of %d  (%s)', m, obj.NumMembers, obj.Labels{m}));
        colorbar;
        clear W k f
    end
end

% Author: Diellor Basha, 2026
