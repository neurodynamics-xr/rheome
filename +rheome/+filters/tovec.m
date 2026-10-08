function V = tovec(W)
% FILTERS.TOVEC  Full-quaternion Dirac field -> physical 3-vector (drop the real w part).
%
%   V = rheome.filters.tovec(W)   % W [4*nV x ...]  ->  V [3*nV x ...]
%
% The Dirac basis carries a full quaternion [w x y z] per vertex; the physical current
% dipole is the imaginary 3-vector (x,y,z). Use this to project a quaternion field (e.g. the
% output of rheome.filters.frame_analysis / rheome.filters.frame_synthesis on a Dirac basis) back to a
% 3-vector for display or geometry. Inverse of rheome.filters.toquat. Port of bst_eigenwavelet('ToVec').
%
% See also: rheome.filters.toquat, rheome.filters.steer, rheome.filters.frame_analysis
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('ToVec'))

    sz = size(W);
    if mod(sz(1), 4) ~= 0
        error('filters:tovec:rows', 'row count (%d) must be a multiple of 4 (quaternion field).', sz(1));
    end
    nV = sz(1) / 4;
    Wr = reshape(W, 4, nV, []);
    V  = reshape(Wr(2:4, :, :), [3*nV, sz(2:end)]);
end

% Author: Diellor Basha, 2026
