function L = dirac_intrinsic_sq(V, F)
% OPERATORS.DIRAC_INTRINSIC_SQ  Intrinsic (immersion/edge-based) quaternionic Dirac squared.
%
%   L = rheome.operators.dirac_intrinsic_sq(V, F)
%
% L = D_int' * MF * D_int  [4nV x 4nV], where MF is the face-area 2-form mass. Ported
% VERBATIM from Brainstorm's tess_operators local_dirac_intrinsic_sq (itself from
% gptoolbox dirac_operator.m; Crane et al., "Spin Transformations of Discrete
% Surfaces"): each 4x4 block is quaternion left-multiplication by the opposite EDGE
% VECTOR (the immersion f) over twice the area. Unlike cotanL(x)I4, its square COUPLES
% the quaternion components -- it carries the intrinsic spin-connection (tangent-frame
% transport); its scalar (w) part equals the cotan Laplacian. Pairs with the extrinsic
% Gauss-map Dirac (rheome.operators.dirac_extrinsic) to give the full rotating-frame operator.
%
% Quaternion order [w,x,y,z], vertex-interleaved (index 4(v-1)+c).
%
% See also: rheome.operators.dirac_intrinsic (D itself), rheome.operators.dirac_frame, rheome.operators.dirac_extrinsic, rheome.eigen.dirac_frame
%
% Author: Diellor Basha, 2026 (port of Crane / gptoolbox via Brainstorm tess_operators)

    [D, MF] = rheome.operators.dirac_intrinsic(V, F);
    L  = D' * MF * D;
    L  = (L + L') / 2;                                      % symmetrize (1e-16 noise)
end

% Author: Diellor Basha, 2026
