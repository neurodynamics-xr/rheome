function W = toquat(V)
% FILTERS.TOQUAT  Physical 3-vector -> pure quaternion Dirac field (real part w = 0).
%
%   W = rheome.filters.toquat(V)   % V [3*nV x ...]  ->  W [4*nV x ...]
%
% Embeds a physical current-dipole 3-vector (x,y,z) as a pure quaternion [0 x y z] per vertex,
% the layout the Dirac basis expects. Feed the result to rheome.filters.frame_analysis / rheome.filters.apply
% on a Dirac basis. Inverse of rheome.filters.tovec. Port of bst_eigenwavelet('ToQuat').
%
% See also: rheome.filters.tovec, rheome.filters.steer, rheome.filters.frame_analysis
%
% Author: Diellor Basha, 2026 (port of bst_eigenwavelet('ToQuat'))

    sz = size(V);
    if mod(sz(1), 3) ~= 0
        error('filters:toquat:rows', 'row count (%d) must be a multiple of 3 (3-vector field).', sz(1));
    end
    nV = sz(1) / 3;
    Vr = reshape(V, 3, nV, []);
    Wr = zeros([4, nV, size(Vr, 3)], 'like', V);
    Wr(2:4, :, :) = Vr;
    W  = reshape(Wr, [4*nV, sz(2:end)]);
end

% Author: Diellor Basha, 2026
