function B = hodge_weights(W, pkg, X, G, keep, Pin)
% INVERSE.HODGE_WEIGHTS  An imaging kernel fused into the Hodge-Dirac modes: sensors straight to each mode's synthesis weights.
%
%   B = rheome.inverse.hodge_weights(W, pkg, X, G, keep, Pin)
%
% W     [3*nIn x C] a free-orientation kernel on the input cortex (x, y, z of each vertex interleaved), e.g. rheome.io.read.nxrkernel
% pkg   the basis on the glued surface (rheome.eigen.hodge_package): mu, kind, P [V x K], mass
% X, G  the glued surface (rheome.geom.glue): positions [V x 3], faces [nG x 3]
% keep  [V x 1] each glued vertex's input vertex, 0 on the bridge (kept first);  Pin [nIn x 3] the input positions
%       (keep = [] and Pin = []: W is already on the glued surface, [3V x C], and is read as it is)
% B     [2K x C] rows 1..K the normal-current weights a (coexact modes; 0 on an exact one: an MEG current has no charge),
%       rows K+1..2K the flow weights b — the editor's fused matrix: w = B*x(t) for a sensor column x(t)
%
% THE CARRY. A kept vertex takes its own rows; a bridge vertex the inverse-distance blend of its 3 nearest cut-out wall vertices
% (in the input positions). THE FORMS. Each channel's current J on the glued surface gives its normal current n.J and its tangent
% part's divergence and curl (the ambient face operator: grad phi_i = N x opposite edge / 2A per face, area-weighted to the
% vertices). THE READINGS under the vertex mass M: e = P_co' M (n.J), o_co = P_co' M curl, o_ex = P_ex' M div. THE WEIGHTS:
% a = G^-1 e with G = P_co' M P_co (the coexact potentials are not orthonormal), b = A^-1 [o_ex; o_co] with A the flow system
% (column k: the divergence and curl of mode k's flow, -grad P_k exact or -n x grad P_k coexact, read the same way), so a flow
% synthesised from b reads back to b whatever the discretisation. Linear in W.
%
% Port of nxr-convert eigen/fused.py (fused_weights, carry_kernel; FaceOperators = the editor's faceGradient.ts), itself the
% editor's hodgeSynthesisWeights. Positions are rounded to single first, as the editor receives them.
%
% See also: rheome.eigen.hodge_package, rheome.io.read.nxrkernel, rheome.spectral.column
%
% Author: Diellor Basha, 2026 (port of nxr-convert fused_weights)
arguments, W (:,:) {mustBeNumeric}, pkg (1,1) struct, X (:,3) double, G (:,3) double, keep (:,1) double = [], Pin (:,3) double = zeros(0, 3), end
C = size(W, 2); nIn = size(Pin, 1);
if isempty(keep), keep = (1:size(X, 1))'; Pin = X; nIn = size(X, 1); end
if size(W, 1) ~= 3 * nIn, error('rheome:hodge_weights', 'a free-orientation kernel has 3 rows a vertex: %d rows for %d vertices', size(W, 1), nIn); end
Wv = permute(reshape(double(W), 3, nIn, C), [2 1 3]);            % [nIn x 3 x C]
m = nnz(keep > 0); kept = keep(1:m);
J = zeros(size(X, 1), 3, C);
J(1:m, :, :) = Wv(kept, :, :);
if size(X, 1) > m                                                % THE CARRY onto the bridge
    cut = setdiff((1:nIn)', kept);
    Pin = double(single(Pin));
    [nn, d] = knnsearch(Pin(cut, :), X(m+1:end, :), 'K', 3);
    w = 1 ./ max(d, 1e-9); w = w ./ sum(w, 2);
    for j = 1:3, J(m+1:end, :, :) = J(m+1:end, :, :) + w(:, j) .* Wv(cut(nn(:, j)), :, :); end
end
ops = i_ops(double(single(X)), G);
P = double(pkg.P); M = pkg.mass(:); K = numel(pkg.mu);
co = find(pkg.kind == 2); ex = find(pkg.kind == 0);
normal = squeeze(sum(J .* ops.n, 2));                            % [V x C]
[dv, cu] = i_divcurl(ops, J);
e = P(:, co)' * (M .* normal);  oco = P(:, co)' * (M .* cu);  oex = P(:, ex)' * (M .* dv);
Gr = P(:, co)' * (M .* P(:, co));
a = zeros(0, C); if ~isempty(co), Gr = (Gr + Gr') / 2; R = chol(Gr); a = R \ (R' \ e); end
order = [ex; co];
A = zeros(numel(order), numel(order));
for k0 = 1:256:numel(order)                                      % THE FLOW SYSTEM, a block of modes at a time
    ks = order(k0:min(k0+255, end));
    g = i_vgrad(ops, P(:, ks));                                  % [V x 3 x w]
    fl = -g;
    c2 = pkg.kind(ks) == 2;
    fl(:, :, c2) = -cross(repmat(ops.n, 1, 1, nnz(c2)), g(:, :, c2), 2);
    [d2, c3] = i_divcurl(ops, fl);
    A(1:numel(ex), k0:k0+numel(ks)-1) = P(:, ex)' * (M .* d2);
    A(numel(ex)+1:end, k0:k0+numel(ks)-1) = P(:, co)' * (M .* c3);
end
b = A \ [oex; oco];
B = zeros(2 * K, C);
B(co, :) = a;
B(K + order, :) = b;
end

function ops = i_ops(P, F)
% the face gradient of hat functions as sparse matrices, the area-weighted face -> vertex mean and the vertex normals (fused.py FaceOperators)
nV = size(P, 1); nF = size(F, 1);
a = P(F(:,1),:); b = P(F(:,2),:); c = P(F(:,3),:);
nrm = cross(b - a, c - a, 2); twoA = vecnorm(nrm, 2, 2); live = twoA > 1e-12;
N = zeros(nF, 3); N(live, :) = nrm(live, :) ./ twoA(live);
area = zeros(nF, 1); area(live) = twoA(live) / 2;
s = twoA; s(~live) = 1;
g = {cross(N, c - b, 2) ./ s, cross(N, a - c, 2) ./ s, cross(N, b - a, 2) ./ s};
for i = 1:3, g{i}(~live, :) = 0; end
rows = repmat((1:nF)', 1, 3); cols = F;
ops.grad = cell(1, 3); ops.curl = cell(1, 3);
for d = 1:3
    gd = [g{1}(:, d) g{2}(:, d) g{3}(:, d)];
    ops.grad{d} = sparse(rows(:), cols(:), gd(:), nF, nV);
    cg = [i_crossd(N, g{1}, d) i_crossd(N, g{2}, d) i_crossd(N, g{3}, d)];   % (N x grad phi) . J = N . (grad phi x J)
    ops.curl{d} = sparse(rows(:), cols(:), cg(:), nF, nV);
end
inc = sparse(cols(:), rows(:), repmat(area, 3, 1), nV, nF);
ws = full(sum(inc, 2)); iw = zeros(nV, 1); iw(ws > 0) = 1 ./ ws(ws > 0);
ops.toV = spdiags(iw, 0, nV, nV) * inc;
vn = inc * N; ln = vecnorm(vn, 2, 2); ln(ln == 0) = Inf;
ops.n = vn ./ ln;
end

function v = i_crossd(A, B, d)
% component d of A x B, row by row
x = cross(A, B, 2); v = x(:, d);
end

function [dv, cu] = i_divcurl(ops, J)
% the tangent current's divergence and curl at the vertices, [V x m] each
J = J - sum(J .* ops.n, 2) .* ops.n;
dv = 0; cu = 0;
for k = 1:3
    Jk = squeeze(J(:, k, :)); if isvector(Jk), Jk = Jk(:); end
    dv = dv + ops.grad{k} * Jk; cu = cu + ops.curl{k} * Jk;
end
dv = ops.toV * dv; cu = ops.toV * cu;
end

function g = i_vgrad(ops, x)
% grad x at the vertices [V x 3 x m]: the area-weighted mean of the face gradients
g = zeros(size(x, 1), 3, size(x, 2));
for k = 1:3, g(:, k, :) = reshape(ops.toV * (ops.grad{k} * x), [], 1, size(x, 2)); end
end
