function phi = poisson(f, S)
% DIFFERENTIAL.POISSON  Solve the surface Poisson equation  K phi = M f  (mean-zero).
%
%   phi = rheome.differential.poisson(f, S)
%
% Solves the cortical Poisson problem PER HEMISPHERE: K phi = M f, where K is the cotan
% Laplace-Beltrami stiffness and M the Galerkin mass (rheome.operators.laplace_beltrami). The
% constant nullspace is removed by projecting the source to mean-zero (mass metric), pinning
% one vertex, and recentering -- faithful to Brainstorm's bst_poisson.
%
% INPUTS:
%   f  [nV x nT] per-vertex scalar source
%   S  surface struct: .Vertices, .Faces (and .Hemi for the per-hemisphere nullspace)
%
% OUTPUT:
%   phi  [nV x nT] per-vertex potential (mean-zero per hemisphere)
%
% See also: rheome.operators.laplace_beltrami, rheome.differential.helmholtz, rheome.differential.divergence
%
% Author: Diellor Basha, 2026

    [K, M] = rheome.operators.laplace_beltrami(S.Vertices, S.Faces, 'galerkin');
    nV = size(S.Vertices, 1);
    hemis = dif_hemis(S, nV);
    phi = zeros(nV, size(f, 2));
    for h = 1:numel(hemis)
        vH = hemis{h};  nh = numel(vH);
        Kh = K(vH, vH);  Mh = M(vH, vH);  fh = f(vH, :);
        fh  = fh - sum(Mh * fh, 1) / sum(Mh(:));   % project source to mean-zero (mass metric)
        rhs = Mh * fh;
        free = 2:nh;                               % pin vertex 1 (nullspace = constant)
        x = zeros(nh, size(f, 2));
        x(free, :) = Kh(free, free) \ rhs(free, :);
        phi(vH, :) = x - mean(x, 1);               % recenter to the mean-zero gauge
    end
end

% Author: Diellor Basha, 2026
