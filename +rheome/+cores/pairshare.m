function r = pairshare(L, Sx)
% CORES.PAIRSHARE  Lagged share of operator L under Sx: ||L Im(Sx) L'||_F / ||L Sx L'||_F, in 0 ... 1 (nxr share).
%
%   r = rheome.cores.pairshare(L, Sx)
%
% = sqrt(tr(Q A Q A') / Re tr(Q Sx Q Sx')), Q = L'L, A = Im(Sx); the denominator is the complex Frobenius norm.
%
% Ports nxr-analysis src/nxr/vortex.py share (:143, commit 3da1320).
%
% See also: rheome.cores.pair, rheome.cores.pairs3
%
% Author: Diellor Basha, 2026
    arguments
        L (:,:) double
        Sx (:,:) {mustBeNumeric}
    end
    Q = L' * L; A = imag(Sx);
    r = sqrt(trace(Q * A * Q * A.') / real(trace(Q * Sx * Q * Sx')));
end
