function J = reconstruct(C, dbasis)
% FORWARD.RECONSTRUCT  Map Dirac mode coefficients back to per-vertex 3-vectors.
%
%   J = rheome.forward.reconstruct(C, dbasis)
%
% Inverse spectral transform. Because Phi is B-orthonormal, a coefficient set maps
% back to a cortical field by  f = Phi * C, keeping the quaternion vector part (the
% imaginary x/y/z slots; the w slot is dropped). Faithful port of bst_dirac
% RECONSTRUCT.
%
% INPUTS:
%   C       [nModes x m] mode coefficients (m columns = frames / channels / sing. vecs)
%   dbasis  Dirac eigenbasis from rheome.eigen.dirac_frame (.Phi, .nVert, .nModes)
%
% OUTPUT:
%   J       [3nV x m] per-vertex 3-vector field(s), rows [x1,y1,z1, x2,y2,z2, ...]
%
% See also: rheome.forward.dirac, rheome.eigen.dirac_frame
%
% Author: Diellor Basha, 2026

    if size(C, 1) ~= dbasis.nModes
        error('forward:reconstruct:size', 'C has %d rows but the basis has %d modes.', size(C,1), dbasis.nModes);
    end
    nV = dbasis.nVert;
    R  = dbasis.Phi * C;          % [4nV x m] quaternion field
    m  = size(R, 2);
    J  = zeros(3*nV, m);
    J(1:3:end, :) = R(2:4:end, :);   % x  (drop w rows 1:4:end)
    J(2:3:end, :) = R(3:4:end, :);   % y
    J(3:3:end, :) = R(4:4:end, :);   % z
end

% Author: Diellor Basha, 2026
