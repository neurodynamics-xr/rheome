function S = root(x, halo)
% SPECTRAL.ROOT  The spectral root of a chunk: its rows' one-sided spectra, from which every view at any time is cut.
%
%   S = rheome.spectral.root(x)              % x [C x w]: one chunk of sensor rows
%   S = rheome.spectral.root(x, halo)        % x [C x (halo + width + halo)]: a block with its context on each side
%
%   S.X      [C x N/2+1] complex, the one-sided spectrum of each row
%   S.N      the transform length, nextpow2(w) (a row shorter than N is mirror-padded)
%   S.halo   the context columns before the block's first own column (n = halo + c reads block column c, 0-based)
%   S.width  the block's own columns (w - 2*halo)
%
% THE LOAD (the editor's chunkFft): each row less its mean over the w loaded columns, mirror-padded to N (sample i >= w reads
% 2(w-1) - i), transformed; the mean put back as N*mean in bin 0, so the views keep the DC. The mean is taken out because a
% float32 transform rounds every butterfly at the offset's size; in double it changes nothing measurable, and the bins equal
% the FFT of the mirror-padded row.
% The transform is circular: a filter's needed context (about 3 cycles of its lowest frequency) must fit in the halo.
%
% Port of the Cortical Flow editor's chunkFft.ts (its spectral root); the views: rheome.spectral.column.
%
% See also: rheome.spectral.column, rheome.spectral.meyerband
%
% Author: Diellor Basha, 2026
arguments, x (:,:) {mustBeNumeric}, halo (1,1) double {mustBeNonnegative, mustBeInteger} = 0, end
x = double(x); [C, w] = size(x);
N = 2 ^ nextpow2(w);
mu = mean(x, 2);
i = 0:N-1; src = i; src(i >= w) = max(0, 2 * (w - 1) - i(i >= w));
F = fft(x(:, src + 1) - mu, N, 2);
F(:, 1) = F(:, 1) + N * mu;
S = struct('X', F(:, 1:N/2+1), 'N', N, 'halo', halo, 'width', w - 2 * halo, 'C', C);
end
