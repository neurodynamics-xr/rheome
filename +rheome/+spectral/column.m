function [x, xa] = column(S, c, H)
% SPECTRAL.COLUMN  The sensor column at one time, cut from the spectral root with a gain on frequency: real and analytic.
%
%   x = rheome.spectral.column(S, c)            % S = rheome.spectral.root(...); c the block's column (0-based); no gain
%   [x, xa] = rheome.spectral.column(S, c, H)   % H [N/2+1 x 1] the gain of each bin (e.g. rheome.spectral.meyerband)
%
%   x   [C x 1] the real column: (1/N) Re sum_k g_k H_k X_k e^{2 pi i k n / N}, g_k = 1 at DC and Nyquist, else 2; n = halo + c
%   xa  [C x 1] complex, the analytic column x + i*Hilbert(x): its imaginary part (1/N) 2 sum_{0<k<N/2} H_k (Xre sin + Xim cos),
%       nothing at DC or Nyquist. abs(xa) the amplitude, angle(xa) the phase of the band at that time.
%
% Then the modes' synthesis weights are w = B * x (rheome.inverse.hodge_weights), and the analytic weights B * xa.
% The phase k*n is reduced mod N in integers before the angle, as the editor does (k*n reaches 10^9).
%
% Port of the editor's spectralColumnCpu and spectralAnalyticCpu (modeView.ts), the CPU twins of its per-frame device pass.
%
% See also: rheome.spectral.root, rheome.spectral.meyerband
%
% Author: Diellor Basha, 2026
arguments, S (1,1) struct, c (1,:) double, H (:,1) double = ones(S.N/2+1, 1), end
N = S.N; k = (0:N/2)';
n = S.halo + c;
th = 2 * pi * mod(k * n, N) / N;                                 % [bins x numel(c)]
g = H * 2; g([1 end]) = H([1 end]);
x = (real(S.X) .* g' * cos(th) - imag(S.X) .* g' * sin(th)) / N;
if nargout > 1
    h = 2 * H; h([1 end]) = 0;
    xa = complex(x, (real(S.X) .* h' * sin(th) + imag(S.X) .* h' * cos(th)) / N);
end
end
