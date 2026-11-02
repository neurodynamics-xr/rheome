function [fa, fb, u] = pair(L, Sx)
% CORES.PAIR  The lagged pair of operator L under the cross-spectrum Sx: sensor filters fa, fb (q a quarter cycle after p) and u = p + i q.
%
%   [fa, fb, u] = rheome.cores.pair(L, Sx)
%
% L [K x C] real (the reading in mode coordinates), Sx [C x C] Hermitian. Reduced SVD L = U diag(s) V'; SV = diag(s) V';
% M = i SV Im(Sx) SV' (Hermitian); w = the eigenvector of the SMALLEST eigenvalue; fa = Re(w) SV, fb = Im(w) SV (row vectors,
% length C); sg = +1 if angle(fa Sx fb.') >= 0 else -1, fb = sg fb (so Im(fa Sx fb.') >= 0: q after p); u = U (Re w + i sg Im w).
%
% The eigenvector's complex phase is FREE (w -> e^{i theta} w rotates fa, fb into each other and u by e^{-i theta}); MATLAB's
% eig and numpy's eigh differ in it, as nxr has no gauge fixing. Compare raw fa, fb, u only after rheome.cores.pairgauge, or
% compare the plane fa'fa + fb'fb and the form fa'fb - fb'fa, share, s3 and the null statistics, which are gauge-free. If
% Im(Sx) = 0 the eigenvalues are all 0 and the pair is arbitrary (only share and s3 ~ 0 mean anything).
%
% Ports nxr-analysis src/nxr/vortex.py pair (:135, commit 3da1320; goldens tests/goldens/vortex-pair).
%
% See also: rheome.cores.pairgauge, rheome.cores.pairshare, rheome.cores.pairs3
%
% Author: Diellor Basha, 2026
    arguments
        L (:,:) double
        Sx (:,:) {mustBeNumeric}
    end
    [U, s, V] = svd(L, 'econ');
    SV = diag(diag(s)) * V';
    M = 1i * (SV * imag(Sx) * SV');
    [W, d] = eig((M + M') / 2, 'vector');
    [~, ix] = min(d);
    w = W(:, ix);
    fa = real(w).' * SV;
    fb = imag(w).' * SV;
    if angle(fa * Sx * fb.') >= 0, sg = 1; else, sg = -1; end
    fb = sg * fb;
    u = U * (real(w) + 1i * sg * imag(w));
end
