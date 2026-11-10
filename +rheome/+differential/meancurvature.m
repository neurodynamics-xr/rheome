function H = meancurvature(S, fg)
% DIFFERENTIAL.MEANCURVATURE  Signed mean curvature at the vertices, from the barycentric-area cotangent Laplacian (nxr's formula).
%
%   H = rheome.differential.meancurvature(S)
%   H = rheome.differential.meancurvature(S, fg)   % reuse rheome.operators.face_gradient(S.Vertices, S.Faces)
%       S   surface: .Vertices [nV x 3], .Faces [nF x 3]
%       H   [nV x 1]  1/length of S.Vertices (1/m for metre positions)
%
% H = -1/2 (DX . n): DX_i = 1/A_i sum_j 1/2 (cot a_ij + cot b_ij)(X_j - X_i) is the cotangent Laplacian of the positions,
% A_i a third of the star area, n the unit area-weighted vertex normal (the normals divergence / curl project with). DX = -2 H n, so
% H = +1/R on a sphere of radius R with outward normals (convex outward positive; the sign follows the face winding). A zero-area
% face (cross product <= 1e-12) contributes no cotangent. The cotan Laplacian is not pointwise convergent at the valence-5
% vertices of an icosphere (~14 % off at every level), but converges in the mean.
% Same formula as rheome.scale.measure_ambient's inline H except for the lumped mass: that one divides by the mass the Laplacian
% routine returns; this one uses nxr's A_i exactly, so the numbers match nxr.
%
% Counterpart: nxr-analysis FaceOperators.mean_curvature (src/nxr/_port/faces.py, commit f7e564f; Diellor, 2026-11-02, one pipeline,
% one mathematics). nxr casts positions to float32 first; this does not (about 1e-7 relative).
%
% Example:
%   [V, F] = rheome.geom.icosphere(4);  H = rheome.differential.meancurvature(struct('Vertices', 0.06*V, 'Faces', F));
%   assert(abs(mean(H) * 0.06 - 1) < 0.01)      % +1/R
%
% See also: rheome.differential.curvaturesource, rheome.differential.divergence, rheome.differential.divcurl
%
% Author: Diellor Basha, 2026

    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V, 1);
    if nargin < 2 || isempty(fg), fg = rheome.operators.face_gradient(V, F); end
    n = fg.W * fg.FaceNormal;  n = n ./ max(sqrt(sum(n.^2, 2)), eps);
    ii = zeros(0, 1);  jj = ii;  ww = ii;  A = zeros(nV, 1);
    for k = 1:3
        i = F(:, k);  j = F(:, mod(k, 3) + 1);  l = F(:, mod(k + 1, 3) + 1);   % angle at l, opposite the edge i-j
        u = V(i, :) - V(l, :);  w = V(j, :) - V(l, :);
        cr = sqrt(sum(cross(u, w, 2).^2, 2));  ok = cr > 1e-12;
        cot = zeros(size(cr));  cot(ok) = sum(u(ok, :) .* w(ok, :), 2) ./ cr(ok);
        ii = [ii; i; j]; jj = [jj; j; i]; ww = [ww; 0.5*cot; 0.5*cot]; %#ok<AGROW>
        A = A + accumarray(F(:, k), 0.5 * cr / 3, [nV 1]);
    end
    L = sparse(ii, jj, ww, nV, nV);  L = L - spdiags(full(sum(L, 2)), 0, nV, nV);
    A(A <= 0) = 1;
    H = -0.5 * sum(((L * V) ./ A) .* n, 2);
end

% Author: Diellor Basha, 2026
