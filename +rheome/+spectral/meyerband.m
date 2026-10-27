function H = meyerband(f, lo, hi, halfWidth)
% SPECTRAL.MEYERBAND  A band's gain at frequencies f: smooth Meyer transitions at its edges, each 2h octaves wide (h = 1/3: from e*2^-h to e*2^h).
%
%   H = rheome.spectral.meyerband(f, lo, hi)                  % f in Hz, the band [lo, hi] Hz
%   H = rheome.spectral.meyerband((0:N/2)' * fs / N, 8, 13)  % one gain a bin of rheome.spectral.root
%
% rise(f) = sin(pi/2 nu((log2(f/lo) + h)/(2h))), fall(f) = cos(pi/2 nu((log2(f/hi) + h)/(2h))), H = rise * fall; nu Meyer's
% auxiliary polynomial (0 to 1 on [0, 1], nu(x) + nu(1-x) = 1). lo = 0 has no rise, hi = Inf no fall; H(0) = 0 unless lo = 0.
% Neighbouring bands sharing an edge partition power: H1^2 + H2^2 = 1 across it.
% Port of the editor's meyerBand and meyerNu (dyadicFrame.ts), the gain its per-frame view applies.
%
% See also: rheome.spectral.column
%
% Author: Diellor Basha, 2026
arguments, f double, lo (1,1) double, hi (1,1) double, halfWidth (1,1) double = 1/3, end
nu = @(x) (x > 0 & x < 1) .* x.^4 .* (35 - 84*x + 70*x.^2 - 20*x.^3) + (x >= 1);
ang = @(x, e) pi / 2 * nu((log2(x / e) + halfWidth) / (2 * halfWidth));
rise = ones(size(f)); fall = ones(size(f));
if lo > 0, rise = sin(ang(f, lo)); end
if isfinite(hi), fall = cos(ang(f, hi)); end
H = rise .* fall;
H(~(f > 0)) = double(lo <= 0);
end
