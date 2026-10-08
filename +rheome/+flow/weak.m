function wk = weak(ctx)
% FLOW.WEAK  Weak (Galerkin) divergence / rotated-source operators from face_gradient.
%   wk = rheome.flow.weak(ctx)
%
% Builds the [V x 3V] sparse operators that act on an interleaved ambient field J [3V x m]:
%   divw  = wk.Bdiv * J      weak divergence  (Galerkin RHS of the Phi-Poisson)
%   vortw = wk.Brot * J      weak vorticity   (Galerkin RHS of the Psi-Poisson)
% These are the exact assemblies used inside rheome.differential.helmholtz, re-expressed as matrices
% so the Helmholtz potentials can be fused into single kernels. The scalar-LBO stiffness that
% inverts them, K = Gx'AGx + Gy'AGy + Gz'AGz, is the same cotan Laplacian whose (L,M) eigenmodes
% the potentials expand in -- so K^+ is diagonal (1/lambda) in that basis.
%
% See also: rheome.flow.potential, rheome.flow.stream, rheome.differential.helmholtz, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026

    fg = ctx.fg;  nV = fg.nV;  nF = fg.nF;
    Gx = fg.Gx; Gy = fg.Gy; Gz = fg.Gz; Nf = fg.FaceNormal; A = fg.FaceArea;

    Fvf = sparse([(1:nF)';(1:nF)';(1:nF)'], fg.Faces(:), 1/3, nF, nV);   % vertex -> face average
    % rotated gradient  Srot = N x G  (weak vorticity source)
    Sx = Nf(:,2).*Gz - Nf(:,3).*Gy;
    Sy = Nf(:,3).*Gx - Nf(:,1).*Gz;
    Sz = Nf(:,1).*Gy - Nf(:,2).*Gx;
    % component selectors: Jx = Ex*J, etc. (interleaved rows [x1 y1 z1 x2 ...])
    Ex = sparse(1:nV, 1:3:3*nV, 1, nV, 3*nV);
    Ey = sparse(1:nV, 2:3:3*nV, 1, nV, 3*nV);
    Ez = sparse(1:nV, 3:3:3*nV, 1, nV, 3*nV);
    AF = spdiags(A, 0, nF, nF);

    wk.Bdiv = Gx'*AF*Fvf*Ex + Gy'*AF*Fvf*Ey + Gz'*AF*Fvf*Ez;   % [V x 3V]
    wk.Brot = Sx'*AF*Fvf*Ex + Sy'*AF*Fvf*Ey + Sz'*AF*Fvf*Ez;   % [V x 3V]
    wk.Fvf  = Fvf;
end

% Author: Diellor Basha, 2026
