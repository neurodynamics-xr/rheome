function laplace_beltrami(SurfaceFile)
% DEMOS.LAPLACE_BELTRAMI  Read a Brainstorm cortex and build/validate the LBO.
%
%   rheome.demos.laplace_beltrami(SurfaceFile)
%
% Step 1 of the atom-design demo. Loads a Brainstorm cortical surface .mat, builds
% the mass matrix and the cotangent Laplace-Beltrami stiffness in pure MATLAB, and
% runs a set of algebraic sanity checks that any correct LBO must pass. Nothing here
% touches the Brainstorm environment beyond reading the surface file from disk.
%
% Needs the toolbox on the path (installed .mltbx, or addpath(genpath(<clone>))); it does
% not edit the path itself.
%
% INPUT:
%   SurfaceFile : path to a Brainstorm surface .mat (e.g. a tess_cortex_*_low.mat).
%
% Uses: rheome.io.read.surface, rheome.operators.mass, rheome.operators.laplace_beltrami
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(SurfaceFile)
        error('demos:laplace_beltrami:args', ['Pass the path to a Brainstorm cortex surface .mat, e.g.\n' ...
            '   rheome.demos.laplace_beltrami(''/path/to/tess_cortex_pial_low.mat'')']);
    end

    fprintf('\n=== Atom-design demo, step 1: Laplace-Beltrami ===\n');

    % --- 1. read the surface ---
    S = rheome.io.read.surface(SurfaceFile);
    fprintf('Surface: %s\n', S.Comment);
    fprintf('  vertices: %d   faces: %d\n', S.nV, S.nF);

    % --- 2. build the operators ---
    M      = rheome.operators.mass(S.Vertices, S.Faces, 'galerkin');
    [L, ~] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces, 'galerkin');

    % --- 3. validation ---
    tol = 1e-9;
    pass = true;

    % surface area from face areas (independent of the mass assembly)
    e1   = S.Vertices(S.Faces(:,2),:) - S.Vertices(S.Faces(:,1),:);
    e2   = S.Vertices(S.Faces(:,3),:) - S.Vertices(S.Faces(:,1),:);
    area = sum(0.5 * sqrt(sum(cross(e1, e2, 2).^2, 2)));

    massTotal = full(sum(M(:)));
    massSym   = i_maxabs(M - M') / max(i_maxabs(M), eps);
    fprintf('\n[mass]\n');
    fprintf('  symmetry (rel)          : %.2e\n', massSym);
    fprintf('  total mass vs area      : %.6g  vs  %.6g  (rel err %.2e)\n', ...
        massTotal, area, abs(massTotal - area)/area);
    pass = pass && massSym < tol && abs(massTotal - area)/area < 1e-6;

    Lsym    = i_maxabs(L - L') / max(i_maxabs(L), eps);
    rowSum  = full(max(abs(sum(L, 2))));
    diagMin = full(min(diag(L)));
    nullRes = full(max(abs(L * ones(S.nV, 1))));
    fprintf('\n[stiffness]\n');
    fprintf('  symmetry (rel)          : %.2e\n', Lsym);
    fprintf('  max |row sum|           : %.2e   (should be ~0)\n', rowSum);
    fprintf('  min diagonal            : %.6g   (should be >= 0)\n', diagMin);
    fprintf('  ||L * 1|| (null mode)   : %.2e   (should be ~0)\n', nullRes);
    pass = pass && Lsym < tol && rowSum < 1e-6 && diagMin >= -tol && nullRes < 1e-6;

    % smallest generalized eigenvalues: PSD spectrum, lambda_0 ~ 0
    k = min(6, S.nV - 2);
    opts.tol = 1e-6;
    lam = sort(real(eigs(L, M, k, 'smallestabs', opts)));
    fprintf('\n[spectrum]  smallest %d generalized eigenvalues L*phi = lambda*M*phi\n', k);
    fprintf('  %s\n', sprintf('%.4g  ', lam));
    fprintf('  lambda_0                : %.2e   (should be ~0)\n', lam(1));
    fprintf('  all eigenvalues >= 0    : %d\n', all(lam > -1e-6));
    pass = pass && abs(lam(1)) < 1e-6 && all(lam > -1e-6);

    fprintf('\n=== %s ===\n\n', i_ternary(pass, 'ALL CHECKS PASSED', 'SOME CHECKS FAILED'));
end

function s = i_ternary(c, a, b)
    if c, s = a; else, s = b; end
end

function m = i_maxabs(A)
% Largest absolute entry of a (sparse) matrix, robust to the all-zero case.
    if nnz(A) == 0, m = 0; else, m = full(max(abs(nonzeros(A)))); end
end

% Author: Diellor Basha, 2026
