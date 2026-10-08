function P = peaks(X, S, opts)
% DETECT.PEAKS  Local maxima of a scalar map on a mesh, per frame, with height, width and value.
%
%   P = rheome.detect.peaks(X, S)
%   P = rheome.detect.peaks(X, S, MinHeight=0.5, Relative=true, MaxPeaks=8)
%   P = rheome.detect.peaks(X, S, MinSeparation=0.01)      % non-maximum suppression within 10 mm
%   P = rheome.detect.peaks(X, S, Tiles=G)                 % also address each peak to a tile (rheome.geom.tiles)
%
% A peak is a vertex whose value beats every vertex in its one-ring. Ties are broken by vertex
% index, so a plateau yields exactly one peak rather than none or all of them. That is the whole
% definition: the map is assumed to be smooth at the scale that matters, and smoothing it is the
% caller's decision, not this function's.
%
% ⭐ PEAKS ARE FOUND AT VERTICES AND ONLY THEN ADDRESSED TO TILES. The maximum of a map within a
% tile is not a peak -- if it sits on the tile's boundary it is the flank of a neighbour's peak, and
% taking one per tile reports a "peak" in every tile in every frame. `Tiles` attaches the address of
% each true peak; it never selects one.
%
% ⭐ POSITION IS REFINED BELOW THE EDGE LENGTH. `.pos` is the value-weighted centroid of the peak and
% its one-ring, using only the part of each neighbour above the ring's lowest value. A peak moving
% at 0.05 m/s and 300 frames/s travels 0.17 mm per frame, a fiftieth of an edge, so a vertex-only
% position is a staircase and every speed read off it is quantised. `.vertex` keeps the lattice
% answer for anything that needs an index.
%
% ⚠⚠ ONE-RING MAXIMA SPLIT A FLAT TOP, so `MinSeparation` is not optional in practice. A basis-truncated
% atom has a plateau with ripples at the mesh's irregularity, and two vertices a few millimetres apart
% both beat their rings at near-equal height. Measured (plant_track_omega.m, a 25 mm mexhat moving at
% 0.25 m/s, no noise): 20% of frames held two or three peaks 1-5 mm apart, which the tracker turned
% into four tracks for one object. Suppressing any peak within MinSeparation of a higher one fixes it.
% Set it to the map's own smoothness scale (the atom width in a plant, r50 = 52 mm on the inverse).
% ⚠ The suppression radius is EUCLIDEAN, so it also suppresses a peak on the opposite bank of a sulcus
% that is millimetres away through the fold and far along the sheet. Below the instrument's resolution
% those two are not separable anyway; for a noiseless map with a small radius it does not arise.
%
% ⚠ WIDTH IS AN AREA, AND IT IS SHARED WHEN PEAKS ARE. `.widthMM` is the equivalent-disc diameter
% 2*sqrt(area/pi) of the connected half-height superlevel set holding the peak -- area from the
% mass matrix, because Euclidean extent across a folded sheet measures the folding. Two peaks whose
% half-height regions touch report the same region. Through an instrument this is a rendering width,
% not a source size (plant_scale_omega.m: ~100 mm recovered whatever was planted below 10 dB).
%
% INPUTS
%   X          [nV x nT] real scalar map (e.g. |curl| of an analytic current)
%   S          surface (.Vertices, .Faces)
%   MinHeight  peaks below this are dropped (0)
%   Relative   true: MinHeight is a fraction of each frame's maximum (false)
%   MaxPeaks   keep at most this many per frame, highest first (Inf)
%   MinSeparation  drop a peak within this Euclidean distance (m) of a higher one (0 = off)
%   Mass       [nV x nV] or [nV x 1] vertex areas for .widthMM; [] = barycentric from S
%   Tiles      a rheome.geom.tiles struct; adds .tile and .node
%
% OUTPUT (struct P)
%   .table      one row per peak: frame vertex pos(1x3) value height widthMM area [tile node]
%   .features   {nT x 1} struct(pos,type,strength) per frame -- rheome.detect.track's frameFeatures
%   .count      [1 x nT] peaks per frame
%   .meanEdge   mean mesh edge length (m)
%
% `height` is value minus the lowest value on the peak's one-ring: a local prominence, cheap and
% monotone in the true topographic prominence for isolated peaks. `value` is the map at the vertex.
%
% See also: rheome.detect.track, rheome.geom.tiles, rheome.geom.tree, rheome.detect.blobscale
%
% Author: Diellor Basha, 2026

    arguments
        X double
        S struct
        opts.MinHeight (1,1) double = 0
        opts.Relative  (1,1) logical = false
        opts.MaxPeaks  (1,1) double {mustBePositive} = Inf
        opts.MinSeparation (1,1) double {mustBeNonnegative} = 0
        opts.Mass = []
        opts.Tiles = []
    end
    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V, 1);
    if size(X, 1) ~= nV
        error('detect:peaks:size', 'X has %d rows; the surface has %d vertices.', size(X,1), nV);
    end
    nT = size(X, 2);

    E = unique(sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2), 'rows');   % undirected edges
    I = [E(:,1); E(:,2)];  J = [E(:,2); E(:,1)];                         % both directions
    nE = numel(I);
    inc = sparse(1:nE, I, 1, nE, nV);                                     % directed edge -> its tail
    Adj = sparse(I, J, true, nV, nV);
    meanEdge = mean(vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2));

    if isempty(opts.Mass),      a = i_vertexarea(V, F);
    elseif isvector(opts.Mass), a = double(opts.Mass(:));
    else,                       a = full(sum(opts.Mass, 2));
    end

    % ⭐ one sparse product decides every vertex of every frame: a vertex is NOT a peak if any
    % neighbour beats it, ties going to the lower index.
    beaten = X(J,:) > X(I,:) | (X(J,:) == X(I,:) & J < I);               % [nE x nT]
    isPk = (inc' * double(beaten)) == 0;                                 % [nV x nT]

    rows = cell(nT, 1);  feats = cell(nT, 1);  count = zeros(1, nT);
    g = graph(Adj);
    for t = 1:nT
        x = X(:, t);
        v = find(isPk(:, t));
        thr = opts.MinHeight;  if opts.Relative, thr = thr * max(x); end
        v = v(x(v) >= thr);
        [~, o] = sort(x(v), 'descend');  v = v(o);
        if opts.MinSeparation > 0 && numel(v) > 1                        % greedy, highest first
            D = pdist2(V(v,:), V(v,:));  keep = true(numel(v), 1);
            for k = 2:numel(v)
                keep(k) = ~any(D(k, 1:k-1) < opts.MinSeparation & keep(1:k-1)');
            end
            v = v(keep);
        end
        v = v(1:min(numel(v), opts.MaxPeaks));
        count(t) = numel(v);
        if isempty(v)
            feats{t} = struct('pos', zeros(0,3), 'type', {{}}, 'strength', zeros(0,1));
            continue
        end
        pos = zeros(numel(v), 3);  lo = zeros(numel(v), 1);  area = zeros(numel(v), 1);
        for k = 1:numel(v)
            nb = find(Adj(:, v(k)));
            lo(k) = min(x(nb));
            w = x([v(k); nb]) - lo(k);                                   % the ring above its floor
            pos(k,:) = (w' * V([v(k); nb], :)) / max(sum(w), realmin);
            % half-height superlevel set, grown from the peak along the mesh
            keep = x >= x(v(k)) / 2;
            c = conncomp(subgraph(g, find(keep)));
            idx = find(keep);  area(k) = sum(a(idx(c == c(idx == v(k)))));
        end
        val = x(v);
        R = table(repmat(t, numel(v), 1), v, pos, val, val - lo, 1e3 * 2 * sqrt(area / pi), area, ...
            'VariableNames', {'frame','vertex','pos','value','height','widthMM','area'});
        if ~isempty(opts.Tiles)
            R.tile = opts.Tiles.tileOf(v);
            R.node = opts.Tiles.node_id(R.tile);
        end
        rows{t} = R;
        feats{t} = struct('pos', pos, 'type', {repmat({'peak'}, numel(v), 1)}, 'strength', val);
    end
    P.table = vertcat(rows{:});
    if isempty(P.table)
        P.table = table(zeros(0,1), zeros(0,1), zeros(0,3), zeros(0,1), zeros(0,1), zeros(0,1), ...
            zeros(0,1), 'VariableNames', {'frame','vertex','pos','value','height','widthMM','area'});
    end
    P.features = feats;
    P.count = count;
    P.meanEdge = meanEdge;
end

% Barycentric vertex areas: a third of each incident triangle.
function a = i_vertexarea(V, F)
    fa = 0.5 * vecnorm(cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2), 2, 2);
    a = accumarray(F(:), repmat(fa / 3, 3, 1), [size(V,1) 1]);
end

% Author: Diellor Basha, 2026
