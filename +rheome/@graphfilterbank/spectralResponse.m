function [H, k] = spectralResponse(obj, n)
% SPECTRALRESPONSE  Member gains against wavenumber. The freqz analogue.
%
%   [H, k] = spectralResponse(gfb)      H [n x M], k [n x 1] = sqrt(lambda), rad/m
%   spectralResponse(gfb)               plots
%
% NOT called freqz, for two reasons: a graph has no frequencies (sqrt(lambda) is a
% WAVENUMBER, rad/m, for a geometric Laplacian), and freqz's 'z' is from the
% z-transform, which has no graph analogue either.
%
% See also: graphfilters, framebounds
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(n), n = 512; end
    lam = linspace(0, obj.Lmax_, n)';
    H   = graphfilters(obj, 'Lambda', lam);
    k   = sqrt(lam);
    if nargout == 0
        plot(k, H, 'LineWidth', 1);
        xlabel('wavenumber k = \surd\lambda  (rad/m)');
        ylabel('gain  g_m(\lambda)');
        title(sprintf('%s, %d members, B/A = %.3g', ...
              obj.Wavelet, obj.NumMembers, framebounds(obj).Tightness));
        grid on;
        clear H k
    end
end

% Author: Diellor Basha, 2026
