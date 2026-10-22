function A = vertex_to_face(V, F, covariant)
% OPERATORS.VERTEX_TO_FACE  Vertex spinors to faces, A (4nF x 4nV): (A psi)_f = sum_{v in f} p_fv psi_v.
%
%   A = rheome.operators.vertex_to_face(V, F, true)    % covariant (Ye et al. adj_vf_cov)
%   A = rheome.operators.vertex_to_face(V, F, false)   % plain (adj_vf): p_fv = 1
%
% p_fv is the unit quaternion turning the area-weighted vertex normal onto the face normal,
% (1 + n_f.n_v, -n_f x n_v)/norm, applied by left multiplication. Unscaled (the quotient problem
% is homogeneous in A). Quaternions [w x y z], interleaved. Port of spinDirac.ts vertexToFace.
%
% Author: Diellor Basha, 2026
arguments, V (:,3) double, F (:,3) double, covariant (1,1) logical = true, end
nF = size(F,1); nV = size(V,1);
x = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
fn = x ./ vecnorm(x,2,2);
vn = zeros(nV,3); for c = 1:3, vn(:,c) = accumarray(F(:), repmat(x(:,c)/2, 3, 1), [nV 1]); end
vn = vn ./ max(vecnorm(vn,2,2), realmin);
I = zeros(48*nF,1); J = I; S = I; k = 0;
for s = 1:3
    v = F(:,s);
    if covariant
        nf = fn; nv = vn(v,:);
        q = [1 + dot(nf,nv,2), -cross(nf,nv,2)]; q = q ./ vecnorm(q,2,2);
    else
        q = repmat([1 0 0 0], nF, 1);
    end
    w = q(:,1); a = q(:,2); b = q(:,3); c = q(:,4);
    B = {w -a -b -c; a w -c b; b c w -a; c -b a w};          % left multiplication matrix, row r, col cc
    for r = 1:4, for cc = 1:4
        n = numel(w); idx = k+(1:n);
        I(idx) = 4*(0:nF-1)'+r; J(idx) = 4*(v-1)+cc; S(idx) = B{r,cc}; k = k + n;
    end, end
end
A = sparse(I(1:k), J(1:k), S(1:k), 4*nF, 4*nV);
end
