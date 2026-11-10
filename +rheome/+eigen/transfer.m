function [T, R] = transfer(subj, tmpl, opts)
% EIGEN.TRANSFER  A subject's Hodge-Dirac mode coefficients carried onto a template's: the transfer matrix T [K_template x K_subject].
%
%   [T, R] = rheome.eigen.transfer(subj, tmpl)
%   [T, R] = rheome.eigen.transfer(subj, tmpl, Sampling="nearest")
%
% subj, tmpl  one struct per side:
%   .pkg     its basis on its glued surface (rheome.eigen.hodge_package): P [V x K], mu, kind, mass
%   .keep    [V x 1] each glued vertex's input vertex, 0 on the bridge (rheome.geom.glue);  .X [V x 3] the glued positions
%   .Pin, .Fin  the input surface (positions as the glue took them, faces 1-based);  .Sphere [nIn x 3] its registered sphere
%   .Left    the axis the left hemisphere lies on: "+y" (Brainstorm SCS, a subject), "-x" (MNI: fsaverage)
% T  coefficients c on the subject (its field f = P c) to the template's: T c. Block-diagonal by kind (exact -> exact,
%    coexact -> coexact): T = A_tmpl R P_subj, with A the template's analysis (exact 2 P' M; coexact G^-1 P' M, G = P' M P on
%    the coexact modes, which are not orthonormal). A field projection, not mode matching: mode signs and rotations inside
%    near-degenerate eigenspaces cancel, and T is exact for what the two bands share.
% R  [V_tmpl x V_subj] sparse, the registration: each template kept vertex pulled from the subject's cortex through the spheres
%    (rheome.geom.sphereregister, a hemisphere at a time; the subject's bridge is not cortex), each template bridge vertex the
%    blend of its 3 nearest cut-out wall vertices, those pulled the same way.
%
% Port of nxr-convert eigen/register.py (transfer_matrix, registration_matrix, fold, bridge_rows); its check: rheome.eigen.transfercheck.
%
% See also: rheome.inverse.hodge_register, rheome.eigen.transfercheck, rheome.geom.sphereregister
%
% Author: Diellor Basha, 2026 (port of nxr-convert transfer_matrix)
arguments, subj (1,1) struct, tmpl (1,1) struct, opts.Sampling (1,1) string = "barycentric", end
R = registration(subj, tmpl, opts.Sampling);
T = fold(tmpl.pkg, R * double(subj.pkg.P)) .* (tmpl.pkg.kind(:) == subj.pkg.kind(:)');
end

function R = registration(s, t, sampling)
ms = nnz(s.keep > 0); mt = nnz(t.keep > 0);
if max(abs(t.Pin(t.keep(1:mt), :) - t.X(1:mt, :)), [], 'all') > 1e-6, error('rheome:transfer', 'the template basis is not on the template surface given (positions differ)'); end
if max(abs(s.Pin(s.keep(1:ms), :) - s.X(1:ms, :)), [], 'all') > 1e-6, error('rheome:transfer', 'the subject basis is not on the subject surface given (positions differ)'); end
hs = rheome.geom.hemispheres(s.Pin, s.Fin, s.Left); ht = rheome.geom.hemispheres(t.Pin, t.Fin, t.Left);
cortex = false(size(s.Pin, 1), 1); cortex(s.keep(1:ms)) = true;
cols = zeros(size(s.Pin, 1), 1); cols(s.keep(1:ms)) = 1:ms;
pull = @(pts) rheome.geom.sphereregister(s.Sphere, s.Fin, hs, cortex, cols, t.Sphere(pts, :), ht(pts), sampling);
[cut, B] = bridgerows(t);
R = [pull(t.keep(1:mt)); B * pull(cut)];
R(size(t.X, 1), size(s.X, 1)) = 0;                               % the subject's bridge columns: never sampled
R = R(1:size(t.X, 1), 1:size(s.X, 1));
end

function [cut, B] = bridgerows(side)
% the input vertices the glue cut, and each bridge vertex's inverse-distance blend of its 3 nearest of them (the editor's carry)
m = nnz(side.keep > 0); nb = size(side.X, 1) - m;
cut = setdiff((1:size(side.Pin, 1))', side.keep(1:m));
if nb == 0, B = sparse(0, numel(cut)); return, end
[nn, d] = knnsearch(side.Pin(cut, :), side.X(m+1:end, :), 'K', 3);
w = 1 ./ max(d, 1e-9); w = w ./ sum(w, 2);
B = sparse(repmat((1:nb)', 1, 3), nn, w, nb, numel(cut));
end

function E = fold(pkg, X)
% the template's analysis of fields X [V x m]: exact 2 P' M x; coexact G^-1 P' M x
P = double(pkg.P); co = pkg.kind(:) == 2;
E = P' * (pkg.mass(:) .* X);
E(~co, :) = 2 * E(~co, :);
if any(co)
    G = P(:, co)' * (pkg.mass(:) .* P(:, co)); G = (G + G') / 2;
    Rc = chol(G); E(co, :) = Rc \ (Rc' \ E(co, :));
end
end
