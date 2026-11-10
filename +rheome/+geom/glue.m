function [X, G, keep, bridge] = glue(P, F, wall)
% GEOM.GLUE  The two hemispheres made one closed surface by anatomy: each medial wall cut out, a bridge grown between the rims.
%
%   [X, G, keep, bridge] = rheome.geom.glue(P, F, wall)
%       P [n x 3], F [nF x 3] (1-based) a surface of two closed pieces;  wall [n x 1] logical, the medial wall (the atlas' unassigned vertices)
%       X [m + R*N x 3]  the glued positions: the kept vertices (ascending input index), then the bridge's R rings of N vertices
%       G [nG x 3]       the glued faces: the two cut hemispheres' faces (input order), then the bridge's
%       keep [m + R*N x 1]  the input vertex of each glued vertex, 0 on the bridge
%       bridge           struct: vertices (R*N), rings (R), area, longestEdge
%
% Each hemisphere's wall is cut as one disk (its largest component, its islands filled, pinches and extra rims cut until one rim
% is left). The two rims are oriented by geometry (along the axis between their centres), resampled to N = round(mean perimeter /
% mean edge) points, aligned by the least-squares shift, and joined by R = round(median gap / mean edge) - 1 rings; each band is
% stitched by the triangulation of least total edge length. Checked as the editor checks: every half-edge once, every edge on two
% faces, one piece, chi = chi1 + chi2 - 2.
%
% Port of the Cortical Flow editor's glue.ts (glueMesh), held equal to it through nxr-convert eigen/glue.py: the same faces, the same bridge
% points (tests/tGlueHodgeParity). Every sum runs in the editor's order, so the positions match bit for bit.
%
% See also: rheome.geom.meshkey, rheome.operators.intrinsic_delaunay, rheome.eigen.hodge_dirac
%
% Author: Diellor Basha, 2026 (port of the Cortical Flow editor's glue.ts)
arguments, P (:,3) double, F (:,3) double, wall (:,1), end
n = size(P,1);
if numel(wall) ~= n, error('rheome:glue', 'glue: the medial-wall mask has %d entries for %d vertices', numel(wall), n); end
wall = logical(wall);
[lab, nb] = i_labels(n, F, true(n,1));
if nb ~= 2, error('rheome:glue', 'glue: the surface has %d pieces, not two', nb); end
pieces = {lab == 1, lab == 2};
chi = zeros(1,2); disks = cell(2,2);
for p = 1:2
    fs = F(pieces{p}(F(:,1)),:);
    [vs, ~, inv] = unique(fs);
    chi(p) = i_euler(reshape(inv, [], 3), numel(vs));
    [disks{p,1}, disks{p,2}] = i_cut_disk(F, pieces{p}, pieces{p} & wall, p);
end
ca = i_centre(P, disks{1,2}); cb = i_centre(P, disks{2,2});
ax = cb - ca;
A = i_orient(P, disks{1,2}, ax); B = i_orient(P, disks{2,2}, ax);
h = i_mean_edge(P, F);
pa = P(A,:); pb = P(B,:);
[cuma, pera] = i_arc(pa); [cumb, perb] = i_arc(pb);
N = max(3, floor((pera + perb) / 2 / h + 0.5));                 % Math.round
SA = zeros(N,3); SB = zeros(N,3);
for i = 0:N-1
    SA(i+1,:) = i_sample(pa, cuma, pera, i / N); SB(i+1,:) = i_sample(pb, cumb, perb, i / N);
end
shift = 0; best = Inf;
for j = 0:N-1
    c = 0;
    for i = 0:N-1
        d = i_dist(SA(i+1,:), SB(mod(i+j,N)+1,:)); c = c + d*d;
    end
    if c < best, best = c; shift = j; end
end
gaps = zeros(N,1);
for i = 0:N-1, gaps(i+1) = i_dist(SA(i+1,:), SB(mod(i+shift,N)+1,:)); end
gaps = sort(gaps);
R = max(1, floor(gaps(floor(N/2)+1) / h + 0.5) - 1);
used = unique([disks{1,1}(:); disks{2,1}(:)]);
new = zeros(n,1); new(used) = 1:numel(used);
m = numel(used); total = m + R*N;
X = zeros(total, 3); X(1:m,:) = P(used,:);
for r = 1:R
    t = r / (R + 1);
    for i = 0:N-1
        a = SA(i+1,:); b = SB(mod(i+shift,N)+1,:);
        X(m + (r-1)*N + i + 1, :) = a + t * (b - a);
    end
end
ring = @(r) m + (r-1)*N + (1:N);
br = i_stitch(new(A), ring(1), X);
for r = 1:R-1, br = [br; i_stitch(ring(r), ring(r+1), X)]; end %#ok<AGROW>
br = [br; i_stitch(ring(R), new(B), X)];
Bf = reshape(br, 3, [])';
la = disks{1,2}; x0 = new(la(1)); x1 = new(la(2));
if any(Bf == x0 & Bf(:,[2 3 1]) == x1, 'all'), Bf = Bf(:,[1 3 2]); end   % the bridge runs against the cut hemispheres' winding
G = [new([disks{1,1}; disks{2,1}]); Bf];
keep = [used; zeros(R*N,1)];
half = reshape(G(:,[1 2 2 3 3 1])', 2, [])';
hk = half(:,1) * total + half(:,2);
if numel(unique(hk)) ~= numel(hk), error('rheome:glue', 'glue: a half-edge is used twice'); end
if ~all(ismember(half(:,2) * total + half(:,1), hk)), error('rheome:glue', 'glue: an edge has one face'); end
[~, nc] = i_labels(total, G, true(total,1));
if nc ~= 1, error('rheome:glue', 'glue: the result is not one piece'); end
want = chi(1) + chi(2) - 2; got = i_euler(G, total);
if got ~= want, error('rheome:glue', 'glue: chi = %d, not %d', got, want); end
e1 = X(Bf(:,2),:) - X(Bf(:,1),:); e2 = X(Bf(:,3),:) - X(Bf(:,1),:);
el = [vecnorm(X(Bf(:,1),:) - X(Bf(:,2),:), 2, 2); vecnorm(X(Bf(:,2),:) - X(Bf(:,3),:), 2, 2); vecnorm(X(Bf(:,3),:) - X(Bf(:,1),:), 2, 2)];
bridge = struct('vertices', R*N, 'rings', R, 'area', 0.5 * sum(vecnorm(cross(e1, e2, 2), 2, 2)), 'longestEdge', max(el));
end

function [lab, nb] = i_labels(n, F, on)
% components through the faces' edges among ON, numbered by their lowest vertex (the editor's union-find, read in vertex order); 0 off ON
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = E(on(E(:,1)) & on(E(:,2)), :);
bins = conncomp(graph(sparse(E(:,1), E(:,2), 1, n, n) + sparse(E(:,2), E(:,1), 1, n, n)))';
lab = zeros(n,1); idx = find(on);
[~, ~, lab(idx)] = unique(bins(idx), 'stable');                  % renumbered in order of first appearance by vertex
nb = max([lab; 0]);
end

function keep = i_largest(n, F, on)
% the members of ON in its largest component (a tie: the one whose lowest vertex comes first)
lab = i_labels(n, F, on);
if ~any(lab), keep = false(n,1); return, end
[~, best] = max(accumarray(lab(lab > 0), 1));                    % first maximum: the earliest component
keep = lab == best;
end

function [kind, out] = i_rims(F, n)
% the faces' half-edges used once, chained into loops in the faces' winding (loops started in order of insertion); 'pinch' for a vertex on two rims
hs = reshape(F(:,[1 2 2 3 3 1])', 2, [])';                       % face by face, corner by corner
k = hs(:,1) * (n+1) + hs(:,2);
bd = hs(~ismember(hs(:,2) * (n+1) + hs(:,1), k), :);
nxt = zeros(n,1); seen = false(n,1);
for i = 1:size(bd,1)
    a = bd(i,1);
    if seen(a), kind = 'pinch'; out = a; return, end
    seen(a) = true; nxt(a) = bd(i,2);
end
kind = 'loops'; out = {};
for i = 1:size(bd,1)
    a0 = bd(i,1);
    if ~seen(a0), continue, end
    loop = a0; seen(a0) = false; a = nxt(a0);
    while a ~= a0
        loop(end+1,1) = a; seen(a) = false; a = nxt(a); %#ok<AGROW>
    end
    out{end+1} = loop; %#ok<AGROW>
end
end

function [faces, rim] = i_cut_disk(F, piece, start, which)
n = numel(piece);
cut = i_largest(n, F, start);
if ~any(cut), error('rheome:glue', 'glue: hemisphere %d has no medial wall in the atlas - the hole is cut by anatomy only', which); end
for it = 1:200
    keep = i_largest(n, F, piece & ~cut);
    cut = cut | (piece & ~keep);                                 % fill its islands: the patch is a disk
    inside = piece(F(:,1)) & ~cut(F(:,1)) & ~cut(F(:,2)) & ~cut(F(:,3));
    faces = F(inside,:);
    [kind, r] = i_rims(faces, n);
    if strcmp(kind, 'pinch'), cut(r) = true; continue, end
    if numel(r) == 1, rim = r{1}; return, end
    [~, o] = sort(cellfun(@numel, r));                           % stable, as Python's sorted: every loop but the (last) longest is cut
    for q = o(1:end-1), cut(r{q}) = true; end
end
error('rheome:glue', 'glue: the medial wall did not become a disk');
end

function d = i_dist(a, b)
x = a(1) - b(1); y = a(2) - b(2); z = a(3) - b(3);
d = sqrt(x*x + y*y + z*z);
end

function c = i_centre(P, L)
c = [0 0 0]; nl = numel(L);
for v = L(:)', c = c + P(v,:) / nl; end
end

function o = i_orient(P, L, ax)
c = i_centre(P, L); s = 0; nl = numel(L);
for i = 1:nl
    p = P(L(i),:) - c; q = P(L(mod(i,nl)+1),:) - c;
    s = s + (p(2)*q(3) - p(3)*q(2)) * ax(1) + (p(3)*q(1) - p(1)*q(3)) * ax(2) + (p(1)*q(2) - p(2)*q(1)) * ax(3);
end
o = L(:); if s < 0, o = flipud(o); end
[~, k] = min(o);                                                 % from its lowest-numbered vertex: the start is not the winding's either
o = o([k:end 1:k-1]);
end

function [cum, per] = i_arc(pts)
nl = size(pts,1); cum = zeros(nl+1,1);
for i = 1:nl, cum(i+1) = cum(i) + i_dist(pts(i,:), pts(mod(i,nl)+1,:)); end
per = cum(nl+1);
end

function x = i_sample(pts, cum, per, s)
L = s * per; k = 1; nl = size(pts,1);
while k < nl && cum(k+1) <= L, k = k + 1; end
a = pts(k,:); b = pts(mod(k,nl)+1,:);
seg = cum(k+1) - cum(k);
if seg > 0, t = (L - cum(k)) / seg; else, t = 0; end
x = a + t * (b - a);
end

function T = i_stitch(A, B0, X)
% the band of least total edge length between two loops running the same way (the editor's stitch); faces as a column, 3 per face
A = A(:); B0 = B0(:);
j0 = 1;
for j = 2:numel(B0)
    if i_dist(X(B0(j),:), X(A(1),:)) < i_dist(X(B0(j0),:), X(A(1),:)), j0 = j; end
end
B = B0([j0:end 1:j0-1]);
m = numel(A); nb = numel(B); W = nb + 1;
a = @(i) A(mod(i,m)+1); b = @(j) B(mod(j,nb)+1);
cost = inf((m+1)*W, 1); frm = zeros((m+1)*W, 1); cost(1) = 0;
for i = 0:m
    for j = 0:nb
        c = cost(i*W + j + 1);
        if i < m
            x = c + i_dist(X(a(i+1),:), X(b(j),:)); k = (i+1)*W + j + 1;
            if x < cost(k) || (x == cost(k) && frm(k) ~= 1), cost(k) = x; frm(k) = 1; end
        end
        if j < nb
            x = c + i_dist(X(a(i),:), X(b(j+1),:)); k = i*W + j + 2;
            if x < cost(k), cost(k) = x; frm(k) = 2; end
        end
    end
end
T = zeros(3*(m+nb), 1); q = 0; i = m; j = nb;
while i > 0 || j > 0
    if frm(i*W + j + 1) == 1, T(q+(1:3)) = [a(i-1) a(i) b(j)]; i = i - 1;
    else, T(q+(1:3)) = [a(i) b(j) b(j-1)]; j = j - 1;
    end
    q = q + 3;
end
end

function h = i_mean_edge(P, F)
s = 0;
for f = 1:size(F,1)
    for k = 1:3, s = s + i_dist(P(F(f,k),:), P(F(f,mod(k,3)+1),:)); end
end
h = s / numel(F);
end

function chi = i_euler(F, n)
e = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');
chi = n - size(e,1) + size(F,1);
end
