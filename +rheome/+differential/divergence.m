function div = divergence(J, S, fg, tangential)
% DIFFERENTIAL.DIVERGENCE  Surface divergence of an ambient source vector field.
%
%   div = rheome.differential.divergence(J, S)
%   div = rheome.differential.divergence(J, S, fg)   % reuse a precomputed rheome.operators.face_gradient(S)
%   div = rheome.differential.divergence(J, S, fg, true)   % tangent-projected field (see TANGENTIAL)
%   div = rheome.differential.divergence(J, S, fg, "curv") % divergence of the normal part J_n n alone (see TANGENTIAL)
%
% Divergence of an UNCONSTRAINED (ambient R^3) per-vertex vector field on the cortical
% surface -- e.g. the reconstructed Dirac / MNE source current. Positive = a source (current
% flowing outward), negative = a sink. Faithful to Brainstorm's bst_divergence ambient branch:
% the flat per-triangle divergence  Gx*Jx + Gy*Jy + Gz*Jz (rheome.operators.face_gradient), area-
% weighted to vertices. A constant field gives 0 even on folds (the hat gradients sum to zero).
% NOTE it uses the FULL 3-vector: the normal component enters via curvature (mean-curvature
% coupling), it is NOT projected out.
%
% INPUTS:
%   J  [3nV x nT] ambient vectors, rows [x1,y1,z1, x2,y2,z2, ...] (one or more time samples)
%   S  surface struct (rheome.io.read.surface / rheome.utils.hemisphere): needs .Vertices, .Faces
%   fg (optional) a precomputed rheome.operators.face_gradient(S.Vertices,S.Faces); pass it to avoid
%      rebuilding the operator per call (e.g. per-frame loops / rheome.detect.operator). S may be [].
%   tangential (optional, default false) which field the divergence is taken of; true/false or one of three names (nxr
%      FaceOperators.div_curl(div=...), Editor flowFieldCpu tangential):
%        false | "3d"    the ambient divergence, no projection (the default, as the Editor's default)
%        true  | "tan"   the tangent flow J - J_n n, what is happening in the surface (div J_tan = div J - div(J_n n))
%        "curv"          the divergence of J_n n alone (J_n = J.n, n the area-weighted vertex normal): the source or sink
%                        the surface's bending makes of a normal current, continuum 2 H J_n (rheome.differential.curvaturesource).
%                        Exactly (linearity)  div3d = div_tan + div_curv.
%
% OUTPUT:
%   div  [nV x nT] per-vertex scalar divergence (sources +, sinks -)
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);  S = struct('Vertices', V, 'Faces', F);
%   J = reshape(repmat([0 0 1], size(V, 1), 1).', [], 1);   % a constant ambient field
%   assert(max(abs(rheome.differential.divergence(J, S))) < 1e-10)
%
% See also: rheome.operators.face_gradient, rheome.differential.curl, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(S.Vertices, S.Faces); end
    if size(J, 1) ~= 3*fg.nV
        error('differential:divergence:size', 'J has %d rows but expects 3*nV = %d.', size(J,1), 3*fg.nV);
    end
    Jx = J(1:3:end, :);  Jy = J(2:3:end, :);  Jz = J(3:3:end, :);
    if nargin < 4 || isempty(tangential), tangential = false; end
    if (ischar(tangential) || isstring(tangential)) && strcmpi(tangential, '3d'), tangential = false; end
    if ischar(tangential) || isstring(tangential) || tangential              % "tan" / "curv" / true need the normals
        vn = fg.W * fg.FaceNormal;  vn = vn ./ max(sqrt(sum(vn.^2, 2)), eps);   % area-weighted vertex normals
        [Jx, Jy, Jz] = dif_part(Jx, Jy, Jz, vn, tangential, 'divergence');
    end
    divF = fg.Gx * Jx + fg.Gy * Jy + fg.Gz * Jz;    % [nF x nT] per-face divergence
    div  = fg.W * divF;                             % [nV x nT] area-weighted to vertices
end

% Author: Diellor Basha, 2026
