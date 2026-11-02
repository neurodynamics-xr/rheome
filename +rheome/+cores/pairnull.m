function S = pairnull(Rh, dof, draws)
% CORES.PAIRNULL  A lag-free Wishart sample covariance: k = max(round(dof), 2) complex Gaussian samples of covariance Rh^2 (nxr wishart).
%
%   S = rheome.cores.pairnull(Rh, dof, draws)
%
% Rh [C x C] = rheome.cores.pairsqrt(Re S_all) (so the null's covariance is Re S_all, no lag). k = evenround(dof), half to even
% as Python's round, but at least 2. X = Rh ((N1 + i N2)/sqrt 2), S = X X'/k.
%
% draws is either a [C x k x 2] array (draws(:,:,1) = N1 the real draw, draws(:,:,2) = N2 the imaginary one: nxr's two
% standard_normal((C, k)) calls, real first), or a RandStream (default: the global stream) from which N1 then N2 are drawn.
% numpy's generator cannot be reproduced; tests feed the stored draws.
%
% Ports nxr-analysis src/nxr/vortex.py wishart (:170, commit 3da1320).
%
% See also: rheome.cores.pairdof, rheome.cores.pairsqrt, rheome.cores.pairshare
%
% Author: Diellor Basha, 2026
    arguments
        Rh (:,:) double
        dof (1,1) double
        draws = RandStream.getGlobalStream
    end
    C = size(Rh, 1);
    k = max(evenround(dof), 2);
    if isa(draws, 'RandStream')
        N1 = randn(draws, C, k); N2 = randn(draws, C, k);
    else
        assert(isequal(size(draws, [1 2 3]), [C k 2]), 'rheome:cores:pairnull', 'draws must be [%d x %d x 2]', C, k);
        N1 = draws(:, :, 1); N2 = draws(:, :, 2);
    end
    X = Rh * ((N1 + 1i * N2) / sqrt(2));
    S = X * X' / k;
end
