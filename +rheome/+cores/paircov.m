function S = paircov(X, m)
% CORES.PAIRCOV  Cross-spectrum of analytic signals over a mask: X(:,m) X(:,m)' / max(sum(m),1) (nxr cov).
%
%   S = rheome.cores.paircov(X, m)
%
% Plain mean over the masked samples: no mean removal, no ddof. X [C x T] complex, m logical or 0/1 [1 x T]; S [C x C] Hermitian.
% Note ' is the conjugate transpose.
%
% Ports nxr-analysis src/nxr/vortex.py cov (commit 3da1320).
%
% Author: Diellor Basha, 2026
    arguments
        X (:,:) {mustBeNumeric}
        m (1,:) {mustBeNumericOrLogical}
    end
    m = logical(m);
    Xm = X(:, m);
    S = Xm * Xm' / max(nnz(m), 1);
end
