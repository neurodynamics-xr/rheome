function f = pairdof(X)
% CORES.PAIRDOF  Effective independent samples per sample of an analytic signal: (sum P)^2 / (N sum P^2) (nxr dof_fraction).
%
%   f = rheome.cores.pairdof(X)
%
% X [C x N] complex, already restricted to the clean columns (X(:, ok), so N = nnz(ok)). P is the mean over channels of
% |fft(X, [], 2)|^2. f is 1 for white noise and about 1/N for a pure line; the Wishart dof is (samples used) * f.
%
% Ports nxr-analysis src/nxr/vortex.py dof_fraction (commit 3da1320).
%
% See also: rheome.cores.pairnull
%
% Author: Diellor Basha, 2026
    arguments
        X (:,:) {mustBeNumeric}
    end
    P = mean(abs(fft(X, [], 2)).^2, 1);
    f = sum(P)^2 / (numel(P) * sum(P.^2));
end
