function basis = modes(L, M, K, varargin)
% EIGEN.MODES  Smallest Laplace-Beltrami eigenmodes as an M-orthonormal basis.
%
%   basis = rheome.eigen.modes(L, M, K)
%   basis = rheome.eigen.modes(L, M, K, 'tol', 1e-6)
%
% Solves the generalized eigenproblem  L phi = lambda M phi  for the K smallest
% eigenvalues and returns an M-ORTHONORMAL basis (phi_i' M phi_j = delta_ij) sorted
% by increasing eigenvalue. This is the spectral basis every atom is projected onto
% (via rheome.filters.apply / rheome.filters.localize) and reconstructed from.
%
% INPUTS:
%   L   [nV x nV] symmetric PSD stiffness (from rheome.operators.laplace_beltrami)
%   M   [nV x nV] symmetric PD mass matrix (from rheome.operators.mass)
%   K   number of modes to compute (clamped to nV-1)
%
% OPTIONS:
%   'tol'    eigs convergence tolerance (default 1e-6)
%
% The eigensolve goes through rheome.eigen.smallest (a port of Brainstorm's bst_eigs_smallest): a small
% NEGATIVE sigma shift so (L - sigma*M) = L + |sigma|*M is SPD/well-conditioned. This is REQUIRED
% -- the naive eigs(L,M,'smallestabs') factorizes the (near-)singular L at sigma=0 and returns
% WRONG low eigenvectors (~40% residual, non-reproducible) for LBO/Dirac/connectome pencils.
%
% OUTPUT (struct basis):
%   .Phi     [nV x K] eigenvectors, M-orthonormal, columns in ascending lambda
%   .Lambda  [K x 1]  eigenvalues (ascending; the smallest ~0)
%   .Mass    [nV x nV] the mass matrix (kept so downstream ops need only 'basis')
%   .nV, .K  counts
%
% Notes:
%   - A whole cortex is two disconnected hemispheres, so the first TWO eigenvalues
%     are ~0 (one constant mode per hemisphere). Both are kept; they are the DC of
%     each hemisphere and pass through low-pass kernels unchanged.
%   - Within a degenerate multiplet the individual eigenvectors are only defined up
%     to rotation; any atom built from a full multiplet is nonetheless well defined.
%
% See also: rheome.operators.laplace_beltrami, rheome.filters.localize, rheome.show.surface
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('tol', 1e-6);
    p.parse(varargin{:});

    nV = size(L, 1);
    K  = min(K, nV - 1);

    opts = struct('tol', p.Results.tol, 'maxit', 1000, 'disp', 0);

    % Robust smallest eigenpairs via the negative-sigma shift (rheome.eigen.smallest = bst_eigs_smallest):
    % NEVER the naive sigma=0 factorization of the singular L, which returns wrong low modes.
    [V, D] = rheome.eigen.smallest(L, M, K, opts);

    lam = real(diag(D));
    [lam, ord] = sort(lam, 'ascend');
    V   = real(V(:, ord));
    lam = max(lam, 0);   % clamp tiny negative numerical zeros (LBO is PSD; the null modes)

    % M-orthonormalize. eigs is NOT guaranteed B-orthonormal across (near-)degenerate
    % multiplets -- e.g. a symmetric sphere, where each degree l has 2l+1 identical
    % eigenvalues and eigs returns an arbitrary, non-orthogonal basis of that eigenspace.
    % Cholesky whitening in the M-inner product restores Phi'*M*Phi = I exactly (mixing
    % only within multiplets, which is allowed), which the spectral reconstruction needs.
    MV = M * V;
    G  = (V' * MV + (V' * MV)') / 2;
    R  = chol(G);
    V  = V / R;

    % ---- how well the mesh resolves the TOP of the basis ----
    % ⚠ COTAN-FEM EIGENVALUES ARE BIASED HIGH, AND THE BIAS GROWS AS (l*h)^2. Measured against the
    % exact spherical spectrum on ico5 (the sphere validation): +7.3e-4 at degree 2 rising to
    % +3.80% at degree 19, the top of a 400-mode basis -- and the per-degree errors collapse onto a
    % single constant once divided by (l*h/R)^2, which is the signature of discretisation rather
    % than of the solver. It therefore cannot be corrected here. But it CAN be declared, and the
    % quantity controlling it is computable on any mesh: how many vertices span the finest
    % half-wavelength the basis carries.
    %
    %   area = sum(M(:))            total surface area (the Galerkin mass matrix sums to it)
    %   h    = sqrt(area/nV)        mean vertex spacing
    %   L_2  = pi/sqrt(lambda_max)  finest half-wavelength in the basis
    %   vpw  = L_2/h                vertices per half-wavelength
    %
    % On ico5 at K = 400 this is 4.5 and the top-of-basis bias is 3.8%. Below ~3 the highest modes
    % are shape-unreliable, not merely inaccurate.
    area  = full(sum(M(:)));
    hMesh = sqrt(area / max(nV,1));
    lmaxK = max(lam);
    vpw   = (pi / sqrt(max(lmaxK, eps))) / hMesh;
    if isfinite(vpw) && vpw < 3
        warning('eigen:modes:resolution', ...
            ['K = %d puts lambda_max = %.4g at %.1f vertices per half-wavelength (mesh h = %.2f mm). ' ...
             'Cotan-FEM eigenvalues are biased HIGH as (l*h)^2, so the top of this basis is ' ...
             'unreliable -- reduce K or refine the mesh. For reference, 4.5 v/hw on ico5 gives +3.8%%.'], ...
            numel(lam), lmaxK, vpw, 1000*hMesh);
    end

    basis = struct('Phi', V, 'Lambda', lam, 'Mass', M, 'nV', nV, 'K', numel(lam), ...
                   'h', hMesh, 'vertsPerHalfWave', vpw);
end

% Author: Diellor Basha, 2026
