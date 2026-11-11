function fs = signedfrequencies(f, half)
% FILTERS.SIGNEDFREQUENCIES  Read a frequency grid (Hz) as signed frequencies.
%
%   fs = rheome.filters.signedfrequencies(f)             % f is a FULL FFT grid, 0 ... fs - df
%   fs = rheome.filters.signedfrequencies(f, 'full')
%   fs = rheome.filters.signedfrequencies(f, 'positive') % f is a ONE-SIDED grid, 0 ... fs/2
%
% The (lambda, f) kernels of this package take SIGNED frequencies in HERTZ. This turns the grid of a
% transform into them:
%
%   'full'      bins 0 ... N-1 at spacing df (as rheome.filters.jspectrum returns): the upper half is
%               negative, f - Fs for f >= Fs/2 with Fs = N*df, the order of MATLAB's fft (+Fs/2 maps to -Fs/2).
%   'positive'  bins 0 ... fs/2 (the one-sided half of a real signal, or the 'positive' Half of a
%               jointfilterbank): every bin is non-negative, so nothing wraps. Wrapping a one-sided
%               grid by Fs = N*df is wrong (Fs is then about fs/2 + df and the upper half of the grid
%               comes out negative), which is why the half is stated and never guessed.
%
% INPUT:   f [1 x N] or [N x 1] frequencies (Hz), uniformly spaced from 0
% OUTPUT:  fs [1 x N] signed frequencies (Hz)
%
% See also: rheome.filters.jspectrum, rheome.filters.travwave, rheome.filters.gabor
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(half), half = 'full'; end
    f = double(f(:)');
    switch lower(char(half))
        case 'positive', fs = f;
        case 'full'
            if numel(f) < 2, fs = f; return; end
            Fs = numel(f) * (f(2) - f(1));
            fs = f - Fs .* (f >= Fs/2);
        otherwise
            error('filters:signedfrequencies:half', 'half must be ''full'' or ''positive'', got ''%s''.', char(half));
    end
end

% Author: Diellor Basha, 2026
