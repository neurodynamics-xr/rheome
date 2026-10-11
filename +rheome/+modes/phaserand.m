function Y = phaserand(Y, rs)
% MODES.PHASERAND  Multivariate Fourier phase randomisation: one random phase per frequency bin, shared by rows.
%
%   Y = rheome.modes.phaserand(Y)            % Y [rows x N] real
%   Y = rheome.modes.phaserand(Y, RandStream('mt19937ar', 'Seed', 6))
%
% Every auto- and cross-spectrum is kept, so every stationary second-order pattern (a steady travelling or
% standing phase gradient included) survives; what is destroyed is structure beyond second order -- the
% surrogate of MS1 G5 and of rheome.scale.measure_atlas's atlas events.
%
% Author: Diellor Basha, 2026

    if nargin < 2, rs = RandStream.getGlobalStream; end
    N = size(Y, 2);  Zf = fft(Y, [], 2);  h = floor((N - 1) / 2);
    ph = exp(2i*pi*rand(rs, 1, h));
    Zf(:, 2:h+1) = Zf(:, 2:h+1) .* ph;  Zf(:, N-h+1:N) = conj(fliplr(Zf(:, 2:h+1)));
    Y = real(ifft(Zf, [], 2));
end

% Author: Diellor Basha, 2026
