function R = pairsqrt(A)
% CORES.PAIRSQRT  Symmetric square root of a PSD matrix, negative eigenvalues clipped to 0 (nxr sqrt_psd): E diag(sqrt(e)) E'.
%
%   R = rheome.cores.pairsqrt(A)
%
% Ports nxr-analysis src/nxr/vortex.py sqrt_psd (commit 3da1320).
%
% Author: Diellor Basha, 2026
    arguments
        A (:,:) double
    end
    [E, e] = eig((A + A') / 2, 'vector');
    R = (E .* sqrt(max(e(:).', 0))) * E';
end
