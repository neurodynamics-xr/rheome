function [V, D] = smallest(A, B, k, opts)
% EIGEN.SMALLEST  k smallest generalized eigenpairs of a (near-)singular symmetric/Hermitian pencil.
%
%   [V, D] = rheome.eigen.smallest(A, B, k)
%   [V, D] = rheome.eigen.smallest(A, B, k, opts)     % opts.tol, opts.maxit, opts.disp
%
% Faithful port of Brainstorm's bst_eigs_smallest. Returns the k smallest generalized eigenpairs
% of (A, B) with B SPD, sorted ascending and real-valued (the pencil is symmetric/Hermitian).
% AVOIDS the sigma=0 shift-invert of eigs(A,B,k,'smallestabs'): at sigma=0 eigs factorizes the
% (near-)singular A, which is ill-posed and returns WRONG low modes for pencils with a near-null
% mode (LBO's constant, Dirac's quaternion null, LB-Connectome's ~1e-9) -- verified: sigma=0 gives
% ~40% eigenpair residual even brute-forced (tol 1e-14, 5000 iters). Instead: estimate lambda_max
% (factorizing only the SPD mass B, so it cannot warn), then use a small NEGATIVE sigma shift so
% (A - sigma*B) = A + |sigma|*B is SPD/well-conditioned while 'nearest sigma' still returns the
% smallest modes (including the near-zero kernel). The -1e-7 -> -1e-4 escalation lifts the shift
% if it lands too close to the kernel; 'smallestabs' is the last-resort fallback.
%
% INPUTS:  A [n x n] symmetric/Hermitian (possibly singular);  B [n x n] SPD mass;  k modes
% OPTIONS: opts.tol (1e-6) .maxit (1000) .disp (0)
% OUTPUTS: V [n x k] eigenvectors (ascending);  D [k x k] diagonal eigenvalues (ascending, real)
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);
%   [L, M] = rheome.operators.laplace_beltrami(V, F);
%   [Phi, D] = rheome.eigen.smallest(L, M, 9);          % the constant, then the l = 1 and l = 2 shells
%   assert(issorted(diag(D)) && norm(Phi' * M * Phi - eye(9)) < 1e-6)
%
% See also: rheome.eigen.modes, rheome.operators.laplace_beltrami   (source: Brainstorm bst_eigs_smallest)
%
% Author: Diellor Basha, 2026 (port of bst_eigs_smallest)

    if nargin < 4 || isempty(opts)
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
    end
    A = (A + A') / 2;                    % symmetrize -> eigs uses the Lanczos (real-spectrum) path
    B = (B + B') / 2;

    % factorization-free spectrum-scale estimate (only an order of magnitude is needed to place
    % the shift); 'largestabs' factorizes only the SPD, well-conditioned mass B, never A.
    estOpts = struct('tol', 1e-2, 'maxit', 300, 'disp', 0);
    try
        lmax = abs(eigs(A, B, 1, 'largestabs', estOpts));
    catch
        lmax = 0;
    end

    if ~isfinite(lmax) || (lmax <= 0)
        [V, D] = eigs(A, B, k, 'smallestabs', opts);          % degenerate estimate: legacy path
    else
        try
            [V, D] = eigs(A, B, k, -1e-7 * lmax, opts);       % small negative shift (SPD factor)
        catch
            try
                [V, D] = eigs(A, B, k, -1e-4 * lmax, opts);   % escalate the shift
            catch
                [V, D] = eigs(A, B, k, 'smallestabs', opts);  % last resort
            end
        end
    end

    d = real(diag(D));
    [d, idx] = sort(d, 'ascend');
    V = V(:, idx);
    D = diag(d);
end

% Author: Diellor Basha, 2026
