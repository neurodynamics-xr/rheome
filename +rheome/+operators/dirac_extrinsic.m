function E = dirac_extrinsic(V, F, N)
% OPERATORS.DIRAC_EXTRINSIC  Relative (extrinsic) Dirac energy block  E = D' *_F D.
%
%   E = rheome.operators.dirac_extrinsic(V, F)       % vertex normals from geometry
%   E = rheome.operators.dirac_extrinsic(V, F, N)    % supply the Gauss map (e.g. Brainstorm VertNormals)
%
% The Liu-Jacobson-Crane (SGP 2017) relative Dirac operator D_N, assembled from the
% Gauss map (vertex normals) in the Galerkin/energy form  E = D' (diag(A_f) (x) I4) D,
% a [4V x 4V] symmetric positive-semidefinite matrix. Per oriented face f = ijk with
% area A_f the first-order relative Dirac is
%   (D_N psi)_f = -(1/2A_f) * sum_{(p,q,r) in cyc(ijk)}  L_{(N_r - N_q)} * psi_p
% where L_v is the 4x4 left-multiplication matrix of the imaginary quaternion v=(x,y,z).
% Quaternion order [w,x,y,z], vertex-interleaved (index 4(v-1)+c). This is the extrinsic
% (E) block of rheome.operators.dirac_frame; it is pure geometry (no leadfield, no plugin).
%
% INPUTS:
%   V   [nV x 3] vertices        F  [nF x 3] triangle indices (consistent winding)
%   N   [nV x 3] unit vertex normals (optional; default = area-weighted from V,F)
%
% OUTPUT:
%   E   [4nV x 4nV] sparse, symmetric, PSD
%
% See also: rheome.operators.dirac_frame, rheome.operators.dirac_intrinsic_sq, rheome.eigen.dirac_frame
%
% Author: Diellor Basha, 2026

    nV = size(V,1);  nF = size(F,1);

    % face normals (unnormalized) and areas
    e1 = V(F(:,2),:) - V(F(:,1),:);
    e2 = V(F(:,3),:) - V(F(:,1),:);
    fn = cross(e1, e2, 2);
    A  = 0.5 * sqrt(sum(fn.^2, 2));                 % [nF x 1] face areas

    if nargin < 3 || isempty(N)
        fnu = fn ./ max(sqrt(sum(fn.^2,2)), eps);   % unit face normals
        N = zeros(nV,3);
        for c = 1:3
            N(:,c) = accumarray(F(:), repmat(A.*fnu(:,c),3,1), [nV 1]);
        end
        N = N ./ max(sqrt(sum(N.^2,2)), eps);       % area-weighted unit vertex normals
    else
        N = N ./ max(sqrt(sum(N.^2,2)), eps);       % enforce unit length
    end

    % assemble D [4nF x 4nV]: block-row f, block-col p gets -L_{(N_r-N_q)} / (2 A_f)
    nnzD = 48*nF;                                   % 3 cyclic terms x 16 per 4x4 block
    ii = zeros(nnzD,1);  jj = ii;  vv = ii;  ptr = 0;
    [rg, cg] = ndgrid(1:4, 1:4);  rg = rg(:);  cg = cg(:);
    for t = 1:nF
        f = F(t,:);  a2 = 2*A(t);
        cyc = [f(1) f(2) f(3);  f(2) f(3) f(1);  f(3) f(1) f(2)];   % rows = (p,q,r)
        r0 = 4*(t-1);
        for s = 1:3
            p = cyc(s,1);  q = cyc(s,2);  r = cyc(s,3);
            blk = -i_qleft(N(r,:) - N(q,:)) / a2;   % 4x4
            c0  = 4*(p-1);
            idx = ptr + (1:16);
            ii(idx) = r0 + rg;
            jj(idx) = c0 + cg;
            vv(idx) = blk(:);
            ptr = ptr + 16;
        end
    end
    D = sparse(ii, jj, vv, 4*nF, 4*nV);

    % star_F = diag(A_f) (x) I4  ;  E = D' star_F D
    W = spdiags(reshape(repmat(A(:)', 4, 1), [], 1), 0, 4*nF, 4*nF);
    E = D' * W * D;
    E = (E + E')/2;                                 % kill round-off asymmetry
end

% ---- 4x4 left-multiplication matrix of a PURE-IMAGINARY quaternion v=(x,y,z) ----
function Lv = i_qleft(v)
    x = v(1);  y = v(2);  z = v(3);
    Lv = [ 0 -x -y -z;
           x  0 -z  y;
           y  z  0 -x;
           z -y  x  0 ];
end

% Author: Diellor Basha, 2026
