function vort = curl(J, S, fg)
% DIFFERENTIAL.CURL  Surface curl (scalar vorticity) of an ambient source vector field.
%
%   vort = rheome.differential.curl(J, S)
%   vort = rheome.differential.curl(J, S, fg)   % reuse a precomputed rheome.operators.face_gradient(S)
%
% On a 2-D surface there is no vector curl: the curl of a field is the SCALAR vorticity
% (the out-of-plane component). For an ambient R^3 per-vertex field it is (grad x J).n per
% face, area-weighted to vertices -- faithful to Brainstorm's bst_curl ambient branch.
% Positive = counter-clockwise rotation seen from outside; a pure gradient (irrotational)
% field gives ~0. Uses the full 3-vector via rheome.operators.face_gradient.
%
% INPUTS:
%   J  [3nV x nT] ambient vectors, rows [x1,y1,z1, x2,y2,z2, ...]
%   S  surface struct (rheome.io.read.surface / rheome.utils.hemisphere): needs .Vertices, .Faces
%   fg (optional) a precomputed rheome.operators.face_gradient(S.Vertices,S.Faces); pass it to avoid
%      rebuilding the operator per call (e.g. per-frame loops / rheome.detect.operator). S may be [].
%
% OUTPUT:
%   vort  [nV x nT] per-vertex scalar vorticity (CCW +, CW -)
%
% See also: rheome.operators.face_gradient, rheome.differential.divergence, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(S.Vertices, S.Faces); end
    if size(J, 1) ~= 3*fg.nV
        error('differential:curl:size', 'J has %d rows but expects 3*nV = %d.', size(J,1), 3*fg.nV);
    end
    Jx = J(1:3:end, :);  Jy = J(2:3:end, :);  Jz = J(3:3:end, :);
    % per-face curl vector  (grad x J):  component d of grad of field Jc is  G_d * Jc
    cvx = fg.Gy * Jz - fg.Gz * Jy;      % (grad x J)_x = dJz/dy - dJy/dz
    cvy = fg.Gz * Jx - fg.Gx * Jz;      % (grad x J)_y = dJx/dz - dJz/dx
    cvz = fg.Gx * Jy - fg.Gy * Jx;      % (grad x J)_z = dJy/dx - dJx/dy
    Nf  = fg.FaceNormal;
    omF = cvx .* Nf(:,1) + cvy .* Nf(:,2) + cvz .* Nf(:,3);   % [nF x nT] vorticity = curl . n
    vort = fg.W * omF;                                        % [nV x nT] area-weighted to vertices
end

% Author: Diellor Basha, 2026
