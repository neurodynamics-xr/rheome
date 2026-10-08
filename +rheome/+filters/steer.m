function Wq = steer(W, q0)
% FILTERS.STEER  Steer a Dirac quaternion field by RIGHT quaternion multiplication.
%
%   Wq = rheome.filters.steer(W, q0)   % W [4*nV x ...], q0 [4x1] unit quaternion -> steered field
%
% The 4-fold multiplicity of the Dirac eigenmodes IS the orientation freedom of the current:
% right-multiplying the full-quaternion field by a unit quaternion q0 rigidly rotates every
% dipole (a global re-aiming of the seed/source direction) while leaving the spatial pattern
% intact. Operates on the quaternion field directly (real w retained) so it composes exactly
% with rheome.filters.frame_analysis / rheome.filters.frame_synthesis on a Dirac basis; project to a
% 3-vector afterwards with rheome.filters.tovec. Port of bst_eigenwavelet('Steer').
%
% ⚠⚠ THIS IS A SINGLE GLOBAL ROTATION AND CANNOT PRODUCE CIRCULATION. Right multiplication by one unit
% quaternion is right-H-linearity -- the symmetry that makes the Dirac spectrum 4-fold degenerate -- so
% every vector in the field turns the same way. A vortex needs the direction to rotate WITH POSITION.
% ⭐ For that use rheome.filters.spin, which CONJUGATES by a quaternion FIELD, q(v) V conj(q(v)); with the axis
% set to the surface normal it turns the in-plane part by any per-vertex angle and leaves the normal
% component exactly unchanged.
%
% INPUTS:
%   W   [4*nV x ...] full-quaternion Dirac field   q0  [4x1] steering quaternion (auto-normalised)
%
% OUTPUT:
%   Wq  [4*nV x ...] steered quaternion field (same size as W)
%
% See also: rheome.filters.tovec, rheome.filters.toquat, rheome.filters.frame_analysis
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('Steer'))

    sz = size(W);
    if mod(sz(1), 4) ~= 0
        error('filters:steer:rows', 'row count (%d) must be a multiple of 4 (quaternion field).', sz(1));
    end
    q0 = q0(:);  q0 = q0 / max(norm(q0), eps);
    R  = i_qrightmat(q0);                       % 4x4 right-multiplication (Hamilton) matrix
    Wq = reshape(R * reshape(W, 4, []), sz);
end

% --- 4x4 matrix R such that R*q == q (x) p  (Hamilton product), q,p,result = [w x y z] ---
function R = i_qrightmat(p)
    pw = p(1); px = p(2); py = p(3); pz = p(4);
    R = [ pw, -px, -py, -pz; ...
          px,  pw,  pz, -py; ...
          py, -pz,  pw,  px; ...
          pz,  py, -px,  pw ];
end

% Author: Diellor Basha, 2026
