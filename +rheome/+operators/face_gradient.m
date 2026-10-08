function fg = face_gradient(V, F)
% OPERATORS.FACE_GRADIENT  Per-face constant-gradient operator of a triangle mesh.
%
%   fg = rheome.operators.face_gradient(V, F)
%
% Builds the discrete surface-gradient primitives shared by the differential operators
% (divergence, curl, Helmholtz). Inside each triangle a per-vertex scalar is the linear
% interpolation of its three values, and the gradient of the vertex hat function is a
% constant tangent vector  grad phi_i = (N x e_i) / (2A). The three ambient components are
% returned as sparse operators Gx/Gy/Gz [nF x nV], so that for a per-vertex scalar field f,
% Gx*f is the per-FACE x-component of grad f (and Gx*Jx + Gy*Jy + Gz*Jz is the divergence).
%
% INPUTS:
%   V  [nV x 3] vertices        F  [nF x 3] triangle vertex indices
%
% OUTPUT (struct fg):
%   .Gx .Gy .Gz    [nF x nV] sparse per-face gradient operators (ambient x/y/z components)
%   .FaceNormal    [nF x 3]  unit face normals (from the vertex order)
%   .FaceArea      [nF x 1]  triangle areas
%   .Faces         [nF x 3]  the faces (as passed)
%   .W             [nV x nF] sparse area-weighted face->vertex average (row-normalized)
%   .nV .nF
%
% See also: rheome.differential.divergence, rheome.differential.curl, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    V = double(V);  F = double(F);
    nV = size(V, 1);  nF = size(F, 1);
    i1 = F(:,1);  i2 = F(:,2);  i3 = F(:,3);

    n    = cross(V(i2,:) - V(i1,:), V(i3,:) - V(i1,:), 2);   % 2A * unit normal
    dblA = max(sqrt(sum(n.^2, 2)), eps);                     % 2 * area
    N    = n ./ dblA;
    g1 = cross(N, V(i3,:) - V(i2,:), 2) ./ dblA;             % grad phi_1  (opposite edge p3-p2)
    g2 = cross(N, V(i1,:) - V(i3,:), 2) ./ dblA;             % grad phi_2  (p1-p3)
    g3 = cross(N, V(i2,:) - V(i1,:), 2) ./ dblA;             % grad phi_3  (p2-p1)

    rows = [(1:nF)'; (1:nF)'; (1:nF)'];  cols = [i1; i2; i3];
    fg.Gx = sparse(rows, cols, [g1(:,1); g2(:,1); g3(:,1)], nF, nV);
    fg.Gy = sparse(rows, cols, [g1(:,2); g2(:,2); g3(:,2)], nF, nV);
    fg.Gz = sparse(rows, cols, [g1(:,3); g2(:,3); g3(:,3)], nF, nV);
    fg.FaceNormal = N;
    fg.FaceArea   = dblA / 2;
    fg.Faces      = F;
    fg.nV = nV;  fg.nF = nF;

    Wi   = sparse(cols, rows, repmat(fg.FaceArea, 3, 1), nV, nF);   % incident face area  [nV x nF]
    fg.W = Wi ./ max(sum(Wi, 2), eps);                             % row-normalized area average
end

% Author: Diellor Basha, 2026
