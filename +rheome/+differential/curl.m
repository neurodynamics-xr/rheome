function vort = curl(J, S, fg, tangent)
% DIFFERENTIAL.CURL  Surface curl (scalar vorticity) of an ambient source vector field.
%
%   vort = rheome.differential.curl(J, S)
%   vort = rheome.differential.curl(J, S, fg)   % reuse a precomputed rheome.operators.face_gradient(S)
%   vort = rheome.differential.curl(J, S, fg, false)   % raw ambient field, no tangent projection (see below)
%   vort = rheome.differential.curl(J, S, fg, "curv")  % curl of the normal part J_n n alone (see below)
%
% On a 2-D surface there is no vector curl: the curl of a field is the SCALAR vorticity
% (the out-of-plane component). For an ambient R^3 per-vertex field it is (grad x J).n per
% face, area-weighted to vertices -- faithful to Brainstorm's bst_curl ambient branch.
% The field is first projected onto the tangent plane at each vertex (J - (J.n)n, n = the area-weighted
% vertex normal fg.W*FaceNormal, normalised), exactly as nxr-analysis (nxr._port.faces.FaceOperators.div_curl)
% and the Cortical Flow Editor (flowFieldCpu, tangential) do: one pipeline, one mathematics. The normal part of
% an ambient field is not rotation on the surface, and on a folded cortex it moves the curl by ~1e-4 relative.
% Positive = counter-clockwise rotation seen from outside; a pure gradient (irrotational)
% field gives ~0. Uses the full 3-vector via rheome.operators.face_gradient.
%
% INPUTS:
%   J  [3nV x nT] ambient vectors, rows [x1,y1,z1, x2,y2,z2, ...]
%   S  surface struct (rheome.io.read.surface / rheome.utils.hemisphere): needs .Vertices, .Faces
%   fg (optional) a precomputed rheome.operators.face_gradient(S.Vertices,S.Faces); pass it to avoid
%      rebuilding the operator per call (e.g. per-frame loops / rheome.detect.operator). S may be [].
%   tangent (optional, default true) which field the curl is taken of; true/false or one of three names (nxr
%      FaceOperators.div_curl(curl=...), Editor flowFieldCpu tangential):
%        true  | "tan"   the tangent flow J - J_n n (what is happening in the surface; the default)
%        false | "3d"    the raw ambient curl, no projection; what rheome.differential.diracparts reads off the Dirac
%                        (kept for that identity); not the Editor's / nxr's default
%        "curv"          the curl of J_n n alone (J_n = J.n, n the area-weighted vertex normal): zero in the continuum, so
%                        what is left is mesh error. Exactly (linearity)  curl3d = curl_tan + curl_curv.
%
% OUTPUT:
%   vort  [nV x nT] per-vertex scalar vorticity (CCW +, CW -)
%
% Example:
%   [V, F] = rheome.geom.icosphere(4);  S = struct('Vertices', V, 'Faces', F);
%   J = reshape(cross(repmat([0 0 1], size(V, 1), 1), V, 2).', [], 1);   % rotation about z
%   vort = rheome.differential.curl(J, S);                    % analytic: 2z on the unit sphere
%   assert(norm(vort - 2*V(:, 3)) / norm(2*V(:, 3)) < 0.01)
%
% See also: rheome.operators.face_gradient, rheome.differential.divergence, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(S.Vertices, S.Faces); end
    if size(J, 1) ~= 3*fg.nV
        error('differential:curl:size', 'J has %d rows but expects 3*nV = %d.', size(J,1), 3*fg.nV);
    end
    Nf  = fg.FaceNormal;
    vn  = fg.W * Nf;  vn = vn ./ max(sqrt(sum(vn.^2, 2)), eps);          % area-weighted vertex normals
    if nargin < 4 || isempty(tangent), tangent = true; end
    Jx = J(1:3:end, :);  Jy = J(2:3:end, :);  Jz = J(3:3:end, :);
    [Jx, Jy, Jz] = dif_part(Jx, Jy, Jz, vn, tangent, 'curl');
    % per-face curl vector  (grad x J):  component d of grad of field Jc is  G_d * Jc
    cvx = fg.Gy * Jz - fg.Gz * Jy;      % (grad x J)_x = dJz/dy - dJy/dz
    cvy = fg.Gz * Jx - fg.Gx * Jz;      % (grad x J)_y = dJx/dz - dJz/dx
    cvz = fg.Gx * Jy - fg.Gy * Jx;      % (grad x J)_z = dJy/dx - dJx/dy
    omF = cvx .* Nf(:,1) + cvy .* Nf(:,2) + cvz .* Nf(:,3);   % [nF x nT] vorticity = curl . n
    vort = fg.W * omF;                                        % [nV x nT] area-weighted to vertices
end

% Author: Diellor Basha, 2026
