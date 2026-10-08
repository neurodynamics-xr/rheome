function [D, MF] = dirac_intrinsic(V, F)
% OPERATORS.DIRAC_INTRINSIC  Crane's intrinsic quaternionic Dirac itself: vertex quaternions -> faces.
%
%   [D, MF] = rheome.operators.dirac_intrinsic(V, F)
%       D   [4nF x 4nV] sparse    (D psi)_f = -(1/2A_f) sum_{i in f} e_i psi_i   (quaternion product)
%       MF  [4nF x 4nF] sparse    face-area mass, so rheome.operators.dirac_intrinsic_sq = D' * MF * D
%
% e_i is the edge OPPOSITE vertex i as a pure quaternion (Crane et al., "Spin Transformations of
% Discrete Surfaces"; gptoolbox dirac_operator.m; Brainstorm tess_operators). Quaternion order
% [w,x,y,z], vertex-interleaved (index 4(v-1)+c), face-interleaved likewise.
%
% ⭐ ON A VECTOR FIELD IT IS CURL, DIVERGENCE AND THE NORMAL COMPONENT'S GRADIENT, EXACTLY. For a pure
% quaternion field X (an ambient 3-vector per vertex, w = 0) and n_f the face normal, per face:
%
%       Re (D X)_f            = -curl X . n_f          the vorticity rheome.differential.curl computes
%       Im (D X)_f . n_f      = -div X                 the divergence rheome.differential.divergence computes
%       Im (D X)_f tangential = +grad (X . n_f)        the in-face gradient of the normal component
%
% because the edge vector IS the hat gradient turned by 90 degrees (e_i = +-2A n x grad(phi_i)), so the
% quaternion product e_i X_i = -e_i.X_i + e_i x X_i splits into the rotated-gradient pairing (curl) and
% the gradient pairing (div). MEASURED, both face orientations, a bumpy sphere: all three hold to
% 3e-16, and the overall sign is the same for either orientation (a derivation by hand got it wrong). ⚠ This is the SWAP of the smooth textbook D = sum e_k d_k (real part
% -div, imaginary normal part curl): this discrete operator is built on the rotated edge vectors.
% rheome.differential.diracparts reads the three parts off and tDiracParts checks the identity to 1e-12.
% ⚠⚠ THE THIRD PART IS NOT ZERO FOR A TANGENT FIELD. n_f is frozen per face, so grad(X . n_f) also
% holds the shape operator acting on X: on the unit sphere (shape operator = I) a tangent field's
% third-part energy converges to the integral of |X|^2 (ratio 0.9904 / 0.9976 / 0.9994 / 0.9998 at
% ico3-6), while its leaked curl falls at second order. On folded cortex this curvature term is large.
% ⚠ D'S SIGN FOLLOWS THE FACE WINDING, THE IDENTITY DOES NOT: reversing every face negates D and n_f
% together: a copy of a surface with every face reversed (same vertices, same face order) stores
% -D, to 3e-11 on a cortical hemisphere.
% ⚠ A real part w contributes a rotated gradient of w to the TANGENTIAL imaginary part only.
%
% See also: rheome.operators.dirac_intrinsic_sq, rheome.differential.diracparts, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026 (port of Crane / gptoolbox via Brainstorm tess_operators)

    V = double(V);  F = double(F);
    nF = size(F, 1);  nV = size(V, 1);
    e1 = V(F(:,2),:) - V(F(:,1),:);
    e2 = V(F(:,3),:) - V(F(:,1),:);
    dblA = sqrt(sum(cross(e1, e2, 2).^2, 2));               % doublearea [nF x 1]
    EV = [zeros(numel(F),1), V(F(:,[2 3 1]),:) - V(F(:,[3 1 2]),:)];   % opposite edge vectors (Im quaternion)
    Q = [1 0 0 0;0 -1  0 0;0 0 -1  0;0  0 0 -1; ...
         0 1 0 0;1  0  0 0;0 0  0 -1;0  0 1  0; ...
         0 0 1 0;0  0  0 1;1 0  0  0;0 -1 0  0; ...
         0 0 0 1;0  0 -1 0;0 1  0  0;1  0 0  0]';
    II = repmat(repmat((0:nF-1)'*4 + (1:4), 1, 4), 3, 1);
    JJ = (repmat(F(:),1,16)-1)*4 + reshape(repmat(1:4,4,1),1,[]);
    D  = sparse(II, JJ, -EV*Q ./ [dblA;dblA;dblA], 4*nF, 4*nV);
    MF = kron(spdiags(dblA/2, 0, nF, nF), speye(4));        % face-area 2-form mass
end

% Author: Diellor Basha, 2026
