function wd = weak_differential(V, F, fg)
% OPERATORS.WEAK_DIFFERENTIAL  Weak-form curl and divergence as [nV x 3nV] operators.
%
%   wd = rheome.operators.weak_differential(V, F)
%   wd = rheome.operators.weak_differential(V, F, fg)   % reuse a precomputed rheome.operators.face_gradient
%
% Returns operators whose ROWS pair a per-vertex test function with the curl (or divergence) of an
% ambient field, WITHOUT ever differentiating that field. On a closed surface, integrating by parts:
%
%     int psi (curl J) dA = -int (N x grad psi) . J dA
%     int psi (div  J) dA = -int (grad psi)     . J dA
%
% so the derivative moves onto psi -- smooth, known, and with an exact hat gradient. For a basis Phi,
%
%     coefficients = Phi' * wd.Curl * J        (NO mass matrix: the pairing already integrates)
%
% ⚠ THIS IS NOT A DIFFERENT DISCRETISATION. Verified in curl_weakform.m: wd.Curl agrees with the
% strong-form face curl integrated per face to 7.2e-16. Integration by parts holds EXACTLY at the
% discrete level -- P1 fields, exact per-face quadrature, closed surface, no boundary term. What it
% actually removes is fg.W, the area-weighted face->vertex average that rheome.differential.curl applies
% before projection. Measured against the analytic curl(zhat x p) = 2z/R on an ico5 sphere:
%
%     Psi' M (fg.W * curl_face)   rel err 4.53e-4     <- rheome.differential.curl, the lumped route
%     Psi' M (curl integrated per face)    1.76e-4
%     Psi' wd.Curl (weak)                  1.76e-4     <- identical, 2.6x better
%
% ⚠ AND IT IS A TRADE, NOT A FREE WIN. Both routes are the same linear map, so they pass noise
% identically -- measured at 1-30% noise the weak form is a flat 5% WORSE, because the lumping is a
% smoothing that blurs signal and smooths noise in equal measure. Use the weak form when accuracy on
% smooth structure matters; keep the lumped route when you want a per-vertex MAP, which is what it
% is actually for.
%
% ⚠ SIGN. The minus from integrating by parts is FOLDED IN, so callers do not have to remember it.
% Getting it wrong is invisible to every magnitude, energy and spectrum check -- it shows up only
% against a reference that carries a sign, as a relative error of exactly 2.000.
%
% ⚠ WHICH DIVERGENCE. Like rheome.differential.divergence, .Div pairs against the AMBIENT field: a constant
% ambient J gives zero, because the hat gradients are in-plane and sum to zero per face. The SURFACE
% divergence of the tangential part differs from it by the mean-curvature coupling 2H(J.N) (sign per
% the OUTWARD normal). Feed it a tangential field and you get the surface answer; feed it an ambient
% one and you get the ambient answer. See facegrad_identities.m.
%
% INPUTS:
%   V  [nV x 3] vertices      F  [nF x 3] faces      fg (optional) rheome.operators.face_gradient(V,F)
%
% OUTPUT (struct wd):
%   .Curl [nV x 3nV] sparse   .Div [nV x 3nV] sparse   -- ambient J is interleaved [x1 y1 z1 x2 ...]
%   .fg   the face-gradient bundle used
%
% See also: rheome.differential.curl, rheome.differential.divergence, rheome.operators.face_gradient, curl_weakform
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(fg), fg = rheome.operators.face_gradient(V, F); end
    nV = fg.nV;  nF = fg.nF;

    N  = fg.FaceNormal;
    dA = spdiags(fg.FaceArea, 0, nF, nF);

    % Rotated gradient per face: (N x grad psi), component by component.
    RGx = spdiags(N(:,2),0,nF,nF)*fg.Gz - spdiags(N(:,3),0,nF,nF)*fg.Gy;
    RGy = spdiags(N(:,3),0,nF,nF)*fg.Gx - spdiags(N(:,1),0,nF,nF)*fg.Gz;
    RGz = spdiags(N(:,1),0,nF,nF)*fg.Gy - spdiags(N(:,2),0,nF,nF)*fg.Gx;

    % Face-average of a P1 field. J is linear inside the triangle and (N x grad psi) is CONSTANT
    % there, so A_f * (N x grad psi)_f . mean(J over f) is the EXACT integral, not a quadrature rule.
    P = sparse(repmat((1:nF)',3,1), double(F(:)), 1/3, nF, nV);

    wd.Curl = i_interleave(-RGx'*dA*P, -RGy'*dA*P, -RGz'*dA*P, nV);
    wd.Div  = i_interleave(-fg.Gx'*dA*P, -fg.Gy'*dA*P, -fg.Gz'*dA*P, nV);
    wd.fg   = fg;
end

% ----- helper: pack three [nV x nV] blocks into [nV x 3nV] for interleaved xyz -----
function C = i_interleave(Qx, Qy, Qz, nV)
    [ix, jx, vx] = find(Qx);
    [iy, jy, vy] = find(Qy);
    [iz, jz, vz] = find(Qz);
    C = sparse([ix; iy; iz], [3*(jx-1)+1; 3*(jy-1)+2; 3*(jz-1)+3], [vx; vy; vz], nV, 3*nV);
end

% Author: Diellor Basha, 2026
