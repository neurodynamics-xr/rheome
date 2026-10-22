function H = hodge_dirac(V, F, varargin)
% OPERATORS.HODGE_DIRAC  Hodge-Dirac d + delta on 0-, 1-, 2-forms over the intrinsic Delaunay triangulation.
%
%   H = rheome.operators.hodge_dirac(V, F)
%       H.K  [N x N] symmetric,  H.M [N x 1] diagonal mass, N = nV + nE + nF  (rows: vertices, edges, faces)
%       H.d0 [nE x nV], H.d1 [nF x nE], H.star0/star1/star2, H.mesh (rheome.operators.intrinsic_delaunay)
%
% D = M^-1 K is self-adjoint in M; its eigenvalues come in pairs +-mu and D^2 is block diagonal with
% the Hodge Laplacians. Solve the pencil (K, M). star0 = lumped mass, star1 = cotan weight,
% star2 = 1/area. Port of nxr-cortical-flow hodgeDirac.ts.
%
% Author: Diellor Basha, 2026
S = rheome.operators.intrinsic_delaunay(V, F, varargin{:});
nV = S.nV; nE = S.nE; nF = S.nF; N = nV + nE + nF;
nxt = @(h) h - mod(h-1,3) + mod(mod(h-1,3)+1,3);
eh = S.edgeHe; head = S.F(nxt(eh)); tail = S.F(eh);
d0 = sparse([1:nE 1:nE]', [head; tail], [ones(nE,1); -ones(nE,1)], nE, nV);
hs = (1:3*nF)'; sg = 1 - 2*(S.edgeHe(S.edgeOf) ~= hs);
d1 = sparse(ceil(hs/3), S.edgeOf, sg, nF, nE);
ln = reshape(S.len(S.edgeOf), 3, [])'; s = sum(ln,2)/2;
star2 = 1 ./ sqrt(max(0, s.*(s-ln(:,1)).*(s-ln(:,2)).*(s-ln(:,3))));
star1 = S.star1;
A = spdiags(star1, 0, nE, nE) * d0;  B = spdiags(star2, 0, nF, nF) * d1;
K = [sparse(nV,nV) A' sparse(nV,nF); A sparse(nE,nE) B'; sparse(nF,nV) B sparse(nF,nF)];
H = struct('K', K, 'M', [S.mass; star1; star2], 'd0', d0, 'd1', d1, 'star0', S.mass, 'star1', star1, 'star2', star2, 'mesh', S, 'N', N);
end

