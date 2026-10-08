function H = helmholtz(J, S)
% DIFFERENTIAL.HELMHOLTZ  Helmholtz-Hodge decomposition of an ambient source vector field.
%
%   H = rheome.differential.helmholtz(J, S)
%
% Splits the ambient (R^3) per-vertex field into orthogonal parts:
%   J  =  grad(Phi)   +   N x grad(Psi)   +   (J.n) n   +   harmonic
%         irrotational     solenoidal          normal        residual
% Phi is the SCALAR POTENTIAL (maxima = sources, minima = sinks); Psi is the STREAM FUNCTION
% (extrema = vortices / anti-vortices). The normal component is pulled out, so the source/sink
% and vortex structure of the TANGENTIAL flow can be read WITHOUT the curvature-coupled normal
% flux. Faithful port of Brainstorm's process_helmholtz Compute: weak divergence / vorticity
% sources -> two Poisson solves (rheome.differential.poisson nullspace handling) -> reconstruction.
%
% INPUTS:
%   J  [3nV x nT] ambient vectors, rows [x1,y1,z1, ...]
%   S  surface struct: .Vertices, .Faces, .VertNormals (and .Hemi for the Poisson nullspace)
%
% OUTPUT (struct H):
%   .Phi .Psi          [nV x nT]  scalar potential (sources/sinks) and stream function (vortices)
%   .Div .Curl         [nV x nT]  strong divergence / vorticity (rheome.differential.divergence/curl)
%   .Virr .Vsol .Hresid .Vtot  [3nV x nT]  irrotational / solenoidal / harmonic / total fields
%   .Fmag .Hmag        [nV x nT]  |J| and |harmonic residual| magnitudes
%   .HarmFrac          [1 x nT]   fraction of field energy in the harmonic residual
%
% See also: rheome.differential.divergence, rheome.differential.curl, rheome.differential.poisson, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026

    fg = rheome.operators.face_gradient(S.Vertices, S.Faces);
    nV = fg.nV;
    if size(J,1) ~= 3*nV
        error('differential:helmholtz:size', 'J has %d rows but expects 3*nV = %d.', size(J,1), 3*nV);
    end
    Nv = S.VertNormals ./ max(vecnorm(S.VertNormals, 2, 2), eps);      % unit vertex normals
    Gx = fg.Gx;  Gy = fg.Gy;  Gz = fg.Gz;  Nf = fg.FaceNormal;  A = fg.FaceArea;

    Jx = J(1:3:end,:);  Jy = J(2:3:end,:);  Jz = J(3:3:end,:);
    Fvf = sparse([(1:fg.nF)';(1:fg.nF)';(1:fg.nF)'], fg.Faces(:), 1/3, fg.nF, nV);   % vertex->face avg
    Jfx = Fvf*Jx;  Jfy = Fvf*Jy;  Jfz = Fvf*Jz;

    % rotated gradient  Srot = N x G  (weak vorticity source)
    Sx = Nf(:,2).*Gz - Nf(:,3).*Gy;
    Sy = Nf(:,3).*Gx - Nf(:,1).*Gz;
    Sz = Nf(:,1).*Gy - Nf(:,2).*Gx;

    % weak sources (Galerkin; already mean-zero per hemisphere) + cotan stiffness from SAME ops
    divw  = Gx'*(A.*Jfx) + Gy'*(A.*Jfy) + Gz'*(A.*Jfz);       % [nV x nT]
    vortw = Sx'*(A.*Jfx) + Sy'*(A.*Jfy) + Sz'*(A.*Jfz);
    K = Gx'*(A.*Gx) + Gy'*(A.*Gy) + Gz'*(A.*Gz);              % weak Laplace-Beltrami stiffness

    hemis = dif_hemis(S, nV);
    % ⭐ ONE factorisation for both potentials: K is shared, and backslash on [divw vortw] returns the
    % same bits as two solves (checked on a reference subject, max diff 0) in a third of the time.
    nT  = size(J, 2);
    PP  = i_hodge_solve(K, [divw vortw], hemis);
    Phi = PP(:, 1:nT);                                       % irrotational potential
    Psi = PP(:, nT+1:end);                                   % solenoidal stream function

    % reconstruct the component fields (per face -> area-weighted to vertices)
    VirrX = fg.W*(Gx*Phi);  VirrY = fg.W*(Gy*Phi);  VirrZ = fg.W*(Gz*Phi);         % grad Phi
    gsx = Gx*Psi;  gsy = Gy*Psi;  gsz = Gz*Psi;
    VsolX = fg.W*(Nf(:,2).*gsz - Nf(:,3).*gsy);                                     % (N x grad Psi)_x
    VsolY = fg.W*(Nf(:,3).*gsx - Nf(:,1).*gsz);
    VsolZ = fg.W*(Nf(:,1).*gsy - Nf(:,2).*gsx);
    Jn = Jx.*Nv(:,1) + Jy.*Nv(:,2) + Jz.*Nv(:,3);                                   % normal component
    HresX = Jx - VirrX - VsolX - Jn.*Nv(:,1);
    HresY = Jy - VirrY - VsolY - Jn.*Nv(:,2);
    HresZ = Jz - VirrZ - VsolZ - Jn.*Nv(:,3);

    % pack
    H.Phi = Phi;  H.Psi = Psi;
    H.Div  = rheome.differential.divergence(J, S);
    H.Curl = rheome.differential.curl(J, S);
    H.Virr   = i_interleave(VirrX, VirrY, VirrZ);
    H.Vsol   = i_interleave(VsolX, VsolY, VsolZ);
    H.Hresid = i_interleave(HresX, HresY, HresZ);
    H.Vtot   = J;
    H.Fmag = sqrt(Jx.^2 + Jy.^2 + Jz.^2);
    H.Hmag = sqrt(HresX.^2 + HresY.^2 + HresZ.^2);
    vertArea = accumarray(fg.Faces(:), repmat(A, 3, 1), [nV 1]) / 3;                % lumped vertex area
    H.HarmFrac = sum(vertArea .* H.Hmag.^2, 1) ./ max(sum(vertArea .* H.Fmag.^2, 1), eps);
end

% ---- solve K x = rhs per hemisphere, mean-zero (weak rhs is already mean-zero -> pin + recenter) ----
function x = i_hodge_solve(K, rhs, hemis)
    x = zeros(size(rhs));
    for h = 1:numel(hemis)
        vH = hemis{h};  nh = numel(vH);
        Kh = K(vH, vH);  rh = rhs(vH, :);
        free = 2:nh;
        xh = zeros(nh, size(rhs, 2));
        xh(free, :) = Kh(free, free) \ rh(free, :);
        x(vH, :) = xh - mean(xh, 1);
    end
end

function V = i_interleave(Vx, Vy, Vz)
    V = zeros(3*size(Vx,1), size(Vx,2));
    V(1:3:end,:) = Vx;  V(2:3:end,:) = Vy;  V(3:3:end,:) = Vz;
end

% Author: Diellor Basha, 2026
