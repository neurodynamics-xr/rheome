function [J, F] = hodge_synthesis(w, pkg, X, G)
% INVERSE.HODGE_SYNTHESIS  The cortical current that Hodge-Dirac synthesis weights make: the map at a time from w = B*x(t).
%
%   [J, F] = rheome.inverse.hodge_synthesis(w, pkg, X, G)
%
% w     [2K x m] synthesis weights (rows 1..K a, K+1..2K b; e.g. B * rheome.spectral.column(...)), real
% pkg   the basis (rheome.eigen.hodge_package);  X, G  the glued surface it is on
% J     [V x 3 x m] the current at each glued vertex: its tangent part the flow -grad A - n x grad S, its normal part s2
% F     struct of the forms, [V x m] each: .A = sum_exact b P (the exact flow's potential), .S = sum_coexact b P (the stream
%       function), .normal = s2 = sum_coexact a P (the normal current)
%
% For an analytic column (complex w) call it on real(w) and imag(w): amplitude and phase at a vertex follow from the two.
% Port of nxr-convert eigen/fused.py synthesis (the editor's HodgeBrushFilter sums and currentFromSums); the inverse of
% rheome.inverse.hodge_weights: reading J back returns w.
%
% See also: rheome.inverse.hodge_weights, rheome.spectral.column
%
% Author: Diellor Basha, 2026 (port of nxr-convert synthesis)
arguments, w (:,:) double, pkg (1,1) struct, X (:,3) double, G (:,3) double, end
K = numel(pkg.mu); P = double(pkg.P);
if size(w, 1) ~= 2 * K, error('rheome:hodge_synthesis', '%d rows of weights for %d modes', size(w, 1), K); end
co = pkg.kind == 2; ex = pkg.kind == 0;
F.normal = P(:, co) * w(co, :);
F.A = P(:, ex) * w(K + find(ex), :);
F.S = P(:, co) * w(K + find(co), :);
[gA, n] = i_vgrad(double(single(X)), G, F.A); gS = i_vgrad(double(single(X)), G, F.S);
J = -gA - cross(repmat(n, 1, 1, size(w, 2)), gS, 2);
J = J + (reshape(F.normal, [], 1, size(w, 2)) - sum(J .* n, 2)) .* n;
end

function [g, n] = i_vgrad(P, F, x)
% grad x at the vertices [V x 3 x m] (area-weighted face gradients) and the area-weighted vertex normals, as hodge_weights
nV = size(P, 1); nF = size(F, 1);
a = P(F(:,1),:); b = P(F(:,2),:); c = P(F(:,3),:);
nr = cross(b - a, c - a, 2); twoA = vecnorm(nr, 2, 2); live = twoA > 1e-12;
N = zeros(nF, 3); N(live, :) = nr(live, :) ./ twoA(live);
s = twoA; s(~live) = 1;
gr = {cross(N, c - b, 2) ./ s, cross(N, a - c, 2) ./ s, cross(N, b - a, 2) ./ s};
for i = 1:3, gr{i}(~live, :) = 0; end
inc = sparse(F(:), repmat((1:nF)', 3, 1), repmat(twoA .* live / 2, 3, 1), nV, nF);
ws = full(sum(inc, 2)); iw = zeros(nV, 1); iw(ws > 0) = 1 ./ ws(ws > 0);
toV = spdiags(iw, 0, nV, nV) * inc;
vn = inc * N; ln = vecnorm(vn, 2, 2); ln(ln == 0) = Inf; n = vn ./ ln;
g = zeros(nV, 3, size(x, 2));
for d = 1:3
    Gd = sparse(repmat((1:nF)', 3, 1), F(:), [gr{1}(:, d); gr{2}(:, d); gr{3}(:, d)], nF, nV);
    g(:, d, :) = reshape(toV * (Gd * x), nV, 1, []);
end
end
