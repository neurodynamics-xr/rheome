function M = divcurl(J, S, fg)
% DIFFERENTIAL.DIVCURL  All six divergence / curl maps of a current at once: tan, curv and 3d.
%
%   M = rheome.differential.divcurl(J, S)
%   M = rheome.differential.divcurl(J, S, fg)      % reuse rheome.operators.face_gradient(S.Vertices, S.Faces)
%       J  [3nV x nT] ambient vectors (rows x1,y1,z1,...);  S  surface (.Vertices .Faces)
%       M  struct of [nV x nT] maps: div_tan, div_curv, div3d, curl_tan, curl_curv, curl3d  (1/length of S.Vertices: per metre for metres)
%
%   div_tan / curl_tan    the tangent flow's own divergence / curl (what happens in the surface);
%   div_curv              div(J_n n): the source / sink the surface's bending makes of a normal current (continuum 2 H J_n,
%                         rheome.differential.curvaturesource); curl_curv = curl(J_n n), zero in the continuum, so mesh error;
%   div3d / curl3d        the raw ambient field, no projection.
% Exactly (linearity, to rounding):  div3d = div_tan + div_curv,  curl3d = curl_tan + curl_curv.
% Defaults of rheome.differential.curl / divergence are unchanged (curl "tan", divergence "3d").
% Counterpart: nxr-analysis FaceOperators.div_curl_all (src/nxr/_port/faces.py, commit f7e564f; Diellor, 2026-11-02: "keep all three").
%
% See also: rheome.differential.divergence, rheome.differential.curl, rheome.differential.meancurvature, rheome.differential.curvaturesource
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(S.Vertices, S.Faces); end
    M.div_tan   = rheome.differential.divergence(J, S, fg, "tan");
    M.div_curv  = rheome.differential.divergence(J, S, fg, "curv");
    M.div3d     = rheome.differential.divergence(J, S, fg, "3d");
    M.curl_tan  = rheome.differential.curl(J, S, fg, "tan");
    M.curl_curv = rheome.differential.curl(J, S, fg, "curv");
    M.curl3d    = rheome.differential.curl(J, S, fg, "3d");
end

% Author: Diellor Basha, 2026
