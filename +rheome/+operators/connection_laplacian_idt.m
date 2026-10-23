function [L, M, S] = connection_laplacian_idt(V, F, nSym)
% OPERATORS.CONNECTION_LAPLACIAN_IDT  Connection (vector) Laplacian on the intrinsic Delaunay triangulation.
%
%   [L, M, S] = rheome.operators.connection_laplacian_idt(V, F, nSym)
%       L  [2nV x 2nV] real symmetric PSD (vertex i: rows 2i-1 = real, 2i = imaginary part)
%       M  [2nV x 1]   mass: the lumped mass of each vertex twice
%       S  the rheome.operators.intrinsic_delaunay struct (with signposts)
%
%   (L u)_i = sum_j w_ij (u_i - e^{i r_ji} u_j),  r_ji = psi_ij + pi - psi_ji,  psi the signpost angle
%   rescaled to 2*pi at its vertex (a boundary vertex keeps its own angle sum). nSym raises every
%   rotation to that power (n-RoSy: 1 vector, 2 line, 4 cross field). Same operator as
%   rheome.operators.connection_laplacian, but on the Delaunay triangulation (weights >= 0, no
%   sliver-driven spectrum top). Port of nxr-cortical-flow connectionLaplacian.ts
%   (Sharp, Soliman & Crane, "The Vector Heat Method", 2019, sec. 5).
%
% Author: Diellor Basha, 2026
arguments, V (:,3) double, F (:,3) double, nSym (1,1) double = 1, end
S = rheome.operators.intrinsic_delaunay(V, F, signposts=true);
nV = S.nV; nE = S.nE;
nxt = @(h) h - mod(h-1,3) + mod(mod(h-1,3)+1,3);
h = S.edgeHe; i = S.F(h); j = S.F(nxt(h));
w = S.star1;
cone = @(v) (~S.boundary(v))*2*pi + S.boundary(v).*S.angleSum(v);  % boundary keeps its own angle sum
psi = @(hh) cone(S.F(hh)) .* S.angle(hh) ./ S.angleSum(S.F(hh));
t = S.twin(h); rev = zeros(nE,1);
rev(t>0) = psi(t(t>0));
b = t == 0;                                                 % boundary: next halfedge's angle plus the corner between
if any(b)
    n = nxt(h(b)); v = S.F(n);
    x = S.len(S.edgeOf(n)); y = S.len(S.edgeOf(nxt(nxt(n)))); o = S.len(S.edgeOf(nxt(n)));
    c = acos(max(-1, min(1, (x.^2 + y.^2 - o.^2)./(2*x.*y))));
    rev(b) = cone(v) .* (S.angle(n) + c) ./ S.angleSum(v);
end
r = nSym * (psi(h) + pi - rev);
Z = @(p, q, A, B) sparse([2*p-1; 2*p-1; 2*p; 2*p], [2*q-1; 2*q; 2*q-1; 2*q], [A; -B; B; A], 2*nV, 2*nV);
L = Z(i, i, w, 0*w) + Z(j, j, w, 0*w) + Z(i, j, -w.*cos(r), -w.*sin(r)) + Z(j, i, -w.*cos(r), w.*sin(r));
M = repelem(S.mass, 2);
end
