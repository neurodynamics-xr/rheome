function b = frame_bounds(frame, lambda)
% FILTERS.FRAME_BOUNDS  Frame bounds A, B and tightness B/A over the eigenvalues.
%
%   b = rheome.filters.frame_bounds(frame, lambda)   % lambda = basis.Lambda [K x 1]
%
% For the per-mode frame response S(lambda) = sum_m g_m(lambda)^2, returns the lower/upper
% frame bounds A = min S, B = max S and Tightness = B/A. A frame with A = B is TIGHT
% (Tightness = 1) and reconstructs exactly: Synthesis(Analysis(F)) = A*F (so exactly F for
% A = 1). A near 0 means the frame leaves part of the spectrum uncovered (widen lrange).
% Port of bst_eigenwavelet('Bounds').
%
% OUTPUT (struct b): .A .B .Tightness
%
% See also: rheome.filters.frame, rheome.filters.frame_gains, rheome.filters.frame_synthesis
%
% Author: Diellor Basha, 2026

    H = rheome.filters.frame_gains(frame, lambda);
    S = sum(H.^2, 2);
    A = min(S);  B = max(S);
    b = struct('A', A, 'B', B, 'Tightness', B / max(A, eps));
end

% Author: Diellor Basha, 2026
