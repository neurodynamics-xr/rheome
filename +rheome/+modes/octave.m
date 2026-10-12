function [z, rate, valid] = octave(Y, fs, edges, mask)
% MODES.OCTAVE  The time factor of the joint field: one jointfilterbank band member as an analytic, decimated signal.
%
%   [z, rate, valid] = rheome.modes.octave(Y, fs, [8 16])
%   [z, rate, valid] = rheome.modes.octave(Y, fs, [8 16], mask)
%
% The raised-cosine band member of rheome.jointfilterbank (rheome.jointfilterbank.band, exactly zero at both
% edges) applied to the one-sided spectrum of each row, x 2: the analytic signal of the band (the joint field's
% octave f, protocol section 14 J1). Decimated to the largest rate fs/m >= 3 x the upper edge; valid is false on
% mask (any sample of the block) and within 3 cycles of the lower edge of either end -- as rheome.modes.analytic,
% which is the Butterworth + Hilbert counterpart of the band-limited field.
%
%   z [rows x T] single complex   rate (Hz)   valid [1 x T] logical
%
% See also: rheome.modes.jointfields, rheome.modes.analytic, rheome.jointfilterbank
%
% Author: Diellor Basha, 2026

    N = size(Y, 2);  if nargin < 4 || isempty(mask), mask = false(1, N); end
    h = rheome.jointfilterbank.band(edges(1), edges(2));
    fk = (0:N-1) * fs / N;  pos = 2:ceil(N/2);  H = zeros(1, N);  H(pos) = 2 * h(2*pi*fk(pos));
    m = max(1, floor(fs / (3 * edges(2))));  rate = fs / m;
    z = ifft(fft(Y, [], 2) .* H, [], 2);  z = single(z(:, 1:m:end));  T = size(z, 2);
    bad = any(reshape([logical(mask(:)'), false(1, m*T - N)], m, T), 1);
    e = ceil(3 / edges(1) * rate);  bad([1:min(e, T), max(1, T-e+1):T]) = true;
    valid = ~bad;
end

% Author: Diellor Basha, 2026
