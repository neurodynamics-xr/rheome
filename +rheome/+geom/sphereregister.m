function R = sphereregister(srcSphere, srcFaces, srcHemi, srcCortex, srcColumns, dstPoints, dstHemi, sampling)
% GEOM.SPHEREREGISTER  Points on a registered sphere sampled from a source surface's cortex: the sparse pull R, a hemisphere at a time.
%
%   R = rheome.geom.sphereregister(srcSphere, srcFaces, srcHemi, srcCortex, srcColumns, dstPoints, dstHemi)
%   R = rheome.geom.sphereregister(..., "nearest")
%
% srcSphere [n x 3]  the source's registered sphere (FreeSurfer sphere.reg), srcFaces [nF x 3] 1-based
% srcHemi   [n x 1] 'L'/'R' (rheome.geom.hemispheres);  srcCortex [n x 1] logical, the vertices that may be sampled
% srcColumns [n x 1] the column of R each source vertex feeds (1-based; e.g. its glued index), 0 where unused
% dstPoints [m x 3]  points on the destination's registered sphere;  dstHemi [m x 1] 'L'/'R'
% R         [m x max(srcColumns)] sparse: each point the blend of the source triangle containing it (barycentric, Workbench's
%           BARYCENTRIC), or else of the cortex triangle it is nearest inside (where the source's medial wall was cut);
%           "nearest": that triangle's vertex of largest weight, weight 1 (FreeSurfer mri_surf2surf's nearest neighbour)
%
% A point is tested against the 12 triangles whose centroids are nearest on the sphere; the one where its smallest
% barycentric weight is largest wins, weights clipped at 0 and renormalised. Port of nxr-convert eigen/register.py
% sphere_registration and _barycentric.
%
% See also: rheome.geom.hemispheres, rheome.eigen.transfer
%
% Author: Diellor Basha, 2026 (port of nxr-convert sphere_registration)
arguments
    srcSphere (:,3) double, srcFaces (:,3) double, srcHemi (:,1) char, srcCortex (:,1) logical, srcColumns (:,1) double
    dstPoints (:,3) double, dstHemi (:,1) char, sampling (1,1) string {mustBeMember(sampling, ["barycentric" "nearest"])} = "barycentric"
end
Xs = srcSphere ./ vecnorm(srcSphere, 2, 2); Yd = dstPoints ./ vecnorm(dstPoints, 2, 2);
I = []; J = []; W = [];
for h = 'LR'
    pts = find(dstHemi == h);
    if isempty(pts), continue, end
    [T, w] = i_barycentric(Xs, srcFaces, Yd(pts, :), srcCortex & srcHemi == h);
    if sampling == "nearest"
        [~, k] = max(w, [], 2); T = T(sub2ind(size(T), (1:size(T, 1))', k)); w = ones(numel(pts), 1);
    end
    I = [I; repmat(pts, size(T, 2), 1)]; J = [J; srcColumns(T(:))]; W = [W; w(:)]; %#ok<AGROW>
end
R = sparse(I, J, W, size(dstPoints, 1), max(srcColumns));
end

function [T, w] = i_barycentric(X, F, Y, ok)
F = F(all(ok(F), 2), :);
C = (X(F(:,1),:) + X(F(:,2),:) + X(F(:,3),:)) / 3; C = C ./ vecnorm(C, 2, 2);
near = knnsearch(C, Y, 'K', 12);
bw = zeros(size(Y, 1), 3); bf = ones(size(Y, 1), 1); bm = -inf(size(Y, 1), 1);
for j = 1:12
    Tj = F(near(:, j), :);
    a = X(Tj(:,1),:); b = X(Tj(:,2),:); c = X(Tj(:,3),:);
    n = cross(b - a, c - a, 2);
    P = Y .* (sum(n .* a, 2) ./ sum(n .* Y, 2));                  % Y projected onto the triangle's plane, from the centre
    w = [sum(cross(c - b, P - b, 2) .* n, 2), sum(cross(a - c, P - c, 2) .* n, 2), sum(cross(b - a, P - a, 2) .* n, 2)] ./ sum(n .* n, 2);
    m = min(w, [], 2); better = m > bm;
    bm(better) = m(better); bw(better, :) = w(better, :); bf(better) = near(better, j);
end
bw = max(bw, 0); w = bw ./ sum(bw, 2);
T = F(bf, :);
end
