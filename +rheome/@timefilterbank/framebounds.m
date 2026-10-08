function b = framebounds(obj)
% FRAMEBOUNDS  Frame bounds of S/2 over the whole axis (0, fs/2), excluding the DC and Nyquist bins.
%
%   b = framebounds(tfb)   .A .B .Tightness (B/A) .Uncovered (fraction of bins with S ~ 0)
%
% The half is the real-signal convention used by the store: for real x,
% sum_m ||w_m||^2 = A * sum x^2 = B * sum x^2 when tight, up to the two excluded bins.
% Tight by construction: A = B = 1 to round-off (tests/tTimeFilterBank.m).
%
% See also: isframetight, rheome.graphfilterbank/framebounds
%
% Author: Diellor Basha, 2026

    H = freqz(obj);
    S = sum(H.^2, 1) / 2;
    S = S(2:end-1);                                   % interior of the one-sided grid
    b = struct('A', min(S), 'B', max(S), 'Tightness', max(S) / max(min(S), eps), ...
               'Uncovered', mean(S <= 1e-12));
end
% Author: Diellor Basha, 2026
