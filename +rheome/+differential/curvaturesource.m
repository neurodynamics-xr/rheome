function c = curvaturesource(J, S, fg)
% DIFFERENTIAL.CURVATURESOURCE  2 H J_n: the continuum value of the divergence the surface's bending gives a normal current.
%
%   c = rheome.differential.curvaturesource(J, S)
%   c = rheome.differential.curvaturesource(J, S, fg)
%       J  [3nV x nT] ambient vectors (rows x1,y1,z1,...);  S  surface (.Vertices .Faces);  fg  optional face_gradient
%       c  [nV x nT]   2 H J_n, H = rheome.differential.meancurvature, J_n = J . n (n the area-weighted unit vertex normal)
%
% The mesh value is rheome.differential.divergence(J, S, fg, "curv") = div(J_n n); it converges to this (on an ellipsoid 70/55/45 mm,
% ico 3-5, to 6 / 3.5 / 2 %). Exactly div3d = div_tan + div_curv, and div_curv ~ 2 H J_n.
% Counterpart: nxr-analysis FaceOperators.div_curv_2h (src/nxr/_port/faces.py, commit f7e564f).
%
% See also: rheome.differential.meancurvature, rheome.differential.divergence, rheome.differential.divcurl
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(S.Vertices, S.Faces); end
    if size(J, 1) ~= 3*fg.nV
        error('differential:curvaturesource:size', 'J has %d rows but expects 3*nV = %d.', size(J,1), 3*fg.nV);
    end
    n = fg.W * fg.FaceNormal;  n = n ./ max(sqrt(sum(n.^2, 2)), eps);
    Jn = n(:,1) .* J(1:3:end, :) + n(:,2) .* J(2:3:end, :) + n(:,3) .* J(3:3:end, :);
    c = 2 * rheome.differential.meancurvature(S, fg) .* Jn;
end

% Author: Diellor Basha, 2026
