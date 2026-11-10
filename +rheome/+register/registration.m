function R = registration(srcSphere, srcFaces, srcHemi, dstPoints, dstHemi, opts)
% REGISTER.REGISTRATION  The sphere registration: points on one registered sphere sampled from another cortex, the sparse matrix R.
%
%   R = rheome.register.registration(srcSphere, srcFaces, srcHemi, dstPoints, dstHemi)
%   R = rheome.register.registration(..., SrcCortex=mask, SrcColumns=cols, Sampling="nearest")
%   R = rheome.register.registration(..., Sampling="shepard", Neighbours=8, Power=2)
%
% srcSphere [n x 3]  the source's registered sphere (FreeSurfer sphere.reg), srcFaces [nF x 3] 1-based
% srcHemi   [n x 1]  'L'/'R' (rheome.register.hemispheres)
% dstPoints [m x 3]  points on the destination's registered sphere;  dstHemi [m x 1] 'L'/'R'
% SrcCortex  [n x 1] logical, the source vertices that may be sampled (default: all)
% SrcColumns [n x 1] the column of R each source vertex feeds (1-based; e.g. its glued index), 0 where unused (default: 1:n)
% Sampling   "barycentric" (default): the blend of the source triangle that holds the point (Workbench's BARYCENTRIC), or else
%            of the cortex triangle it is nearest inside (where the source's medial wall was cut);
%            "nearest": that triangle's vertex of largest weight, weight 1 (FreeSurfer mri_surf2surf's nearest neighbour);
%            "shepard": Brainstorm's ("Project on default anatomy"): rheome.register.shepard_weights of the Neighbours source
%            vertices nearest the point (a hemisphere's cortex vertices), Power
% R [m x max(SrcColumns)] sparse
%
% A point is tested against the 12 triangles whose centroids are nearest on the sphere; the one where its smallest barycentric
% weight is largest wins, weights clipped at 0 and renormalised. As nxr.register.sphere.registration.
%
% See also: rheome.register.hemispheres, rheome.register.shepard_weights, rheome.register.mode_registration
%
% Author: Diellor Basha, 2026
arguments
    srcSphere (:,3) double, srcFaces (:,3) double, srcHemi (:,1) char, dstPoints (:,3) double, dstHemi (:,1) char
    opts.SrcCortex (:,1) logical = true(size(srcSphere, 1), 1)
    opts.SrcColumns (:,1) double = (1:size(srcSphere, 1))'
    opts.Sampling (1,1) string {mustBeMember(opts.Sampling, ["barycentric" "nearest" "shepard"])} = "barycentric"
    opts.Neighbours (1,1) double {mustBeInteger, mustBePositive} = 8
    opts.Power (1,1) double = 2
end
cortex = opts.SrcCortex; cols = opts.SrcColumns;
Xs = srcSphere ./ vecnorm(srcSphere, 2, 2); Yd = dstPoints ./ vecnorm(dstPoints, 2, 2);
I = []; J = []; W = [];
for h = 'LR'
    pts = find(dstHemi == h);
    if isempty(pts), continue, end
    if opts.Sampling == "shepard"
        pool = find(cortex & srcHemi == h);
        [T, w] = rheome.register.shepard_weights(Xs(pool, :), Yd(pts, :), opts.Neighbours, opts.Power);
        T = pool(T);
    else
        [T, w] = i_barycentric(Xs, srcFaces, Yd(pts, :), cortex & srcHemi == h);
    end
    if opts.Sampling == "nearest"
        [~, k] = max(w, [], 2); T = T(sub2ind(size(T), (1:size(T, 1))', k)); w = ones(numel(pts), 1);
    end
    I = [I; repmat(pts, size(T, 2), 1)]; J = [J; cols(T(:))]; W = [W; w(:)]; %#ok<AGROW>
end
R = sparse(I, J, W, size(dstPoints, 1), max(cols));
end

function [T, w] = i_barycentric(X, F, Y, ok)
F = F(all(ok(F), 2), :);
C = (X(F(:,1),:) + X(F(:,2),:) + X(F(:,3),:)) / 3; C = C ./ vecnorm(C, 2, 2);
near = knnsearch(C, Y, 'K', min(12, size(F, 1)));
bw = zeros(size(Y, 1), 3); bf = ones(size(Y, 1), 1); bm = -inf(size(Y, 1), 1);
for j = 1:size(near, 2)
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
