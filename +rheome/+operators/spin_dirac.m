function S = spin_dirac(V, F)
% OPERATORS.SPIN_DIRAC  Spin Diracs on quaternion spinors per FACE of a closed triangle mesh.
%
%   S = rheome.operators.spin_dirac(V, F)
%       S.extrinsic    [4nF x 4nF]  D_f = 1/2 sum_j E_ij phi_j - H_i phi_i   (spin transformations; kills constants)
%       S.intrinsicHY  [4nF x 4nF]  D   = 1/2 sum_j E_ij phi_j              (Hoffmann & Ye, = D_f + H)
%       S.intrinsicYe  [4nF x 4nF]  D   = 1/2 sum_j cos(theta_ij/2) E_ij phi_j   (Ye et al. 2018; edge lengths only)
%       S.H, S.area    [nF x 1]     integrated mean curvature and area per face
%       S.mass         [4nF x 4nF]  area (x) I4: the eigenproblem is D phi = mu mass phi
%
% Port of nxr-cortical-flow editor/src/engine/toolbox/differential/spinDirac.ts (hyperedges, assemble).
% Hyperedge across the halfedge e (ccw in face i) to face j:  E_ij = |e| tan(theta/2) + e, quaternion [w x y z],
% sin(theta) = (n_i x n_j).e/|e|, cos(theta) = n_i.n_j. Rows 4(i-1)+r, interleaved [w x y z].
% Hoffmann & Ye arXiv 1802.06278; Ye, Diamanti, Tang, Guibas & Hoffmann, CGF 37(5) 2018.
%
% shortcut: the face-spinor eigenproblem has spurious small eigenvalues (the cell doubling); solve on vertex
% spinors through rheome.operators.spin_vertex_to_face, as Ye et al. do.
nF = size(F,1); nV = size(V,1);
e1 = V(F(:,2),:)-V(F(:,1),:); e2 = V(F(:,3),:)-V(F(:,1),:);
c = cross(e1, e2, 2); l = vecnorm(c, 2, 2);
if any(~(l > 0)), error('rheome:spin_dirac:degenerate', 'a zero-area face'); end
n = c./l; area = l/2;
h = (1:3*nF)'; f = ceil(h/3); k = mod(h-1,3);
a = F(sub2ind(size(F), f, k+1)); b = F(sub2ind(size(F), f, mod(k+1,3)+1));
T = sparse(a, b, h, nV, nV);
twin = full(T(sub2ind([nV nV], b, a)));
if any(twin == 0) || nnz(T) ~= 3*nF, error('rheome:spin_dirac:open', 'the mesh must be closed, oriented and manifold'); end
j = ceil(twin/3);
ev = V(b,:)-V(a,:); len = vecnorm(ev, 2, 2);
ni = n(f,:); nj = n(j,:);
sn = dot(cross(ni, nj, 2), ev, 2)./len; cs = dot(ni, nj, 2);
th = sn./(1+cs);                                   % tan(theta/2)
hc = 1./sqrt(1+th.^2);
Eq = [len.*th, ev];                                % [nH x 4]
H = accumarray(f, 0.5*len.*th, [nF 1]);
% left multiplication matrix of (w,x,y,z): rows r, cols c
% left multiplication matrix of (w,x,y,z), row-major
L = @(q) [q(:,1) -q(:,2) -q(:,3) -q(:,4)  q(:,2) q(:,1) -q(:,4) q(:,3)  q(:,3) q(:,4) q(:,1) -q(:,2)  q(:,4) -q(:,3) q(:,2) q(:,1)];
B = L(Eq);                                         % row-major 4x4 per halfedge
ri = 4*(f-1) + floor((0:15)/4) + 1;   ci = 4*(j-1) + mod(0:15,4) + 1;
build = @(s, d) sparse(ri(:), ci(:), reshape((0.5*s).*B, [], 1), 4*nF, 4*nF) + d;
I4 = @(x) spdiags(repelem(x, 4), 0, 4*nF, 4*nF);
S.extrinsic = build(ones(3*nF,1), I4(-H));
S.intrinsicHY = build(ones(3*nF,1), sparse(4*nF,4*nF));
S.intrinsicYe = build(hc, sparse(4*nF,4*nF));
S.H = H; S.area = area; S.mass = I4(area);
end
