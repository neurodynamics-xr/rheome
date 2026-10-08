function H = frame_gains(frame, lambda)
% FILTERS.FRAME_GAINS  Evaluate a frame's member gains on eigenvalues -> H [K x M].
%
%   H = rheome.filters.frame_gains(frame, lambda)   % lambda = basis.Lambda [K x 1]
%
% Column m is g_m(lambda), the spectral gain of frame member m. This is the frame
% analogue of a single filter kernel evaluated on the spectrum. Port of
% bst_eigenwavelet('Evaluate').
%
% See also: rheome.filters.frame, rheome.filters.frame_bounds
%
% Author: Diellor Basha, 2026

    lambda = double(lambda(:));
    M = numel(frame.g);
    H = zeros(numel(lambda), M);
    for m = 1:M
        v = frame.g{m}(lambda);
        H(:, m) = v(:);
    end
end

% Author: Diellor Basha, 2026
