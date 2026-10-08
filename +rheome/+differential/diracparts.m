function P = diracparts(J, S, D)
% DIFFERENTIAL.DIRACPARTS  Curl, divergence and normal gradient of a field, read off the intrinsic Dirac.
%
%   P = rheome.differential.diracparts(J, S)
%   P = rheome.differential.diracparts(J, S, D)      % reuse rheome.operators.dirac_intrinsic(S.Vertices, S.Faces)
%       J   [3nV x nT] ambient vectors (rows x1,y1,z1,...), or [4nV x nT] quaternions [w x y z]
%       S   surface: .Vertices .Faces
%
% One sparse product, D*J, and the three parts of the result per face (rheome.operators.dirac_intrinsic):
%
%   P.curlF   [nF x nT]    -Re (D J)                  = curl J . n_f    (vorticity, CCW +)
%   P.divF    [nF x nT]    -Im (D J) . n_f            = div J           (source +)
%   P.gradNF  [3nF x nT]   Im (D J) tangential        = grad (J . n_f)  (rows x1,y1,z1 per face)
%   P.curl P.div           [nV x nT] area-weighted to vertices, exactly as rheome.differential.curl /
%                          rheome.differential.divergence do, so they are interchangeable with them
%   P.energy  struct: face-area integrals of curlF^2, divF^2, |gradNF|^2 (each [1 x nT]) and
%             .total = ||D J||^2 in the face mass, which is their sum
%
% ⭐ WHY THROUGH THE DIRAC. The three parts are one operator, so a field's Dirac energy
% ||D J||^2 = J' * dirac_intrinsic_sq * J splits EXACTLY into rotational, divergent and normal-bending
% shares -- and for a basis of dirac_intrinsic_sq (tau = 0) that split holds mode by mode and so band
% by band: a wavelet band's energy is lambda-weighted and the three shares add up to it.
% ⚠ ONLY AT tau = 0. The production relative Dirac (tau = 0.5) adds the extrinsic Gauss-map block, so
% its eigenmodes are not singular vectors of D_int; D_int of such a mode is not "that mode's curl".
% ⚠⚠ gradNF IS EXTRINSIC, NOT "THE NORMAL COMPONENT'S" ALONE: with n_f frozen per face it also carries
% the shape operator acting on J, so a purely TANGENT field has gradNF energy = the integral of |S J|^2
% (unit sphere: the integral of |J|^2, to 0.02% at ico6). Curl and div are the clean parts; the third
% is normal variation plus the field's bending with the surface.
% ⚠⚠ ON REAL CORTEX, READ CURL AND DIV AT A SCALE, NEVER PER VERTEX. Neighbouring vertex normals differ
% by ~30 deg on a reference subject's mesh, and pointwise curl and div leak into each other at mesh scale: a
% pure source reads curl/div energy 0.13-0.30; low-passed to > 60 mm it reads 0.002-0.017
% (dirac_parts_omega.m). This is the same number rheome.differential.curl gives -- not a Dirac property.
% rheome.differential.helmholtz separates at every scale (91% / 0.7%).
% ⚠ With a quaternion input, a real part w enters gradNF as a rotated gradient of w (not a vector part).
%
% See also: rheome.operators.dirac_intrinsic, rheome.differential.curl, rheome.differential.divergence, rheome.differential.helmholtz
%
% Author: Diellor Basha, 2026

    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V, 1);  nF = size(F, 1);
    if nargin < 3 || isempty(D), D = rheome.operators.dirac_intrinsic(V, F); end
    nT = size(J, 2);
    if size(J, 1) == 3*nV
        Q = zeros(4*nV, nT);  Q(setdiff(1:4*nV, 1:4:4*nV), :) = J;
    elseif size(J, 1) == 4*nV
        Q = J;
    else
        error('differential:diracparts:size', 'J has %d rows; expects 3*nV = %d or 4*nV = %d.', size(J,1), 3*nV, 4*nV);
    end
    Y = D * Q;                                                % [4nF x nT]
    w = Y(1:4:end, :);  ix = Y(2:4:end, :);  iy = Y(3:4:end, :);  iz = Y(4:4:end, :);
    n = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
    dblA = sqrt(sum(n.^2, 2));  n = n ./ dblA;  A = dblA / 2;
    dn = ix .* n(:,1) + iy .* n(:,2) + iz .* n(:,3);
    P.curlF = -w;
    P.divF = -dn;
    gx = ix - dn .* n(:,1);  gy = iy - dn .* n(:,2);  gz = iz - dn .* n(:,3);
    P.gradNF = zeros(3*nF, nT);  P.gradNF(1:3:end, :) = gx;  P.gradNF(2:3:end, :) = gy;  P.gradNF(3:3:end, :) = gz;
    Wi = sparse(F(:), repmat((1:nF)', 3, 1), repmat(A, 3, 1), nV, nF);
    Wv = Wi ./ max(sum(Wi, 2), eps);                          % rheome.operators.face_gradient's .W
    P.curl = Wv * P.curlF;  P.div = Wv * P.divF;
    P.energy.curl = A' * P.curlF.^2;
    P.energy.div = A' * P.divF.^2;
    P.energy.gradN = A' * (gx.^2 + gy.^2 + gz.^2);
    P.energy.total = P.energy.curl + P.energy.div + P.energy.gradN;
end

% Author: Diellor Basha, 2026
