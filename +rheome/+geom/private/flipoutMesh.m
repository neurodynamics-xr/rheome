function M = flipoutMesh(V, F)
% FLIPOUTMESH  The halfedge mesh and the input geometry rheome.geom.flipOut works on (private).
%
%   M = flipoutMesh(V, F)      V [nV x 3], F [nF x 3] 1-based, oriented, manifold
%
% A port of geometry-central's mesh construction and of the geometry quantities FlipOut reads, in the same element
% order and the same floating-point order (see rheome.geom.flipOut for the source, commit and licence):
%   - ManifoldSurfaceMesh::ManifoldSurfaceMesh(polygons)        (manifold_surface_mesh.cpp): halfedges, edges, the
%     boundary loops, v.halfedge(), f.halfedge(). Edge e has halfedges 2e-1 (e.halfedge(), the one met first in the
%     face list) and 2e (its twin); everything is 1-based here, in the same order as geometry-central's 0-based ids.
%   - EmbeddedGeometryInterface::computeEdgeLengths / computeCornerAngles / computeFaceAreas
%     (embedded_geometry_interface.cpp), IntrinsicGeometryInterface::computeVertexAngleSums /
%     computeCornerScaledAngles / computeHalfedgeVectorsInFace / computeHalfedgeVectorsInVertex
%     (intrinsic_geometry_interface.cpp): the input geometry (a VertexPositionGeometry) that traceGeodesic reads.
%   - IntrinsicTriangulation::IntrinsicTriangulation (intrinsic_triangulation.cpp: the intrinsic copy's
%     vertexAngleSums, from edge lengths) and SignpostIntrinsicTriangulation::SignpostIntrinsicTriangulation
%     (signpost_intrinsic_triangulation.cpp: the initial signpost angles), so each path starts from a copy.
%
% Fields: nV nF nE nH nB, next tail face (per halfedge; face > nF is a boundary loop), vHe (per vertex), fHe (per
% face and loop), isInt (halfedge interior), vBnd (vertex on the boundary), P (positions), L (edge lengths),
% angleSum, hvF [nH x 2], hvV [nH x 2] (input geometry), angleSumI, signpost (the intrinsic triangulation's start).
%
% Author: Diellor Basha, 2026 (a port of geometry-central, MIT, Copyright (c) 2017-2019 Nicholas Sharp and the
% geometry-central contributors; licence: +rheome/+geom/LICENSE-geometry-central)

    V = double(V);  F = double(F);
    nF = size(F, 1);  nV = max(F(:));
    if size(F, 2) ~= 3, error('rheome:flipOut:faces', 'flipOut needs a triangle mesh.'); end
    if nV ~= size(V, 1) || numel(unique(F(:))) ~= nV
        error('rheome:flipOut:vertices', 'Every vertex must be in a face (geometry-central refuses unreferenced vertices).');
    end

    % == ManifoldSurfaceMesh(polygons): the face list walked in order; an edge is created where its first halfedge is met
    K = 3 * nF;
    tails = reshape(F.', [], 1);                                   % occurrence k = 3(f-1)+j: face f, corner j
    tips = reshape(F(:, [2 3 1]).', [], 1);
    fOf = repelem((1:nF).', 3);
    [~, firstK, eOf] = unique(sort([tails tips], 2), 'rows', 'stable');
    nE = numel(firstK);  nH = 2 * nE;
    isFirst = false(K, 1);  isFirst(firstK) = true;
    if any(accumarray(eOf, 1) > 2), error('rheome:flipOut:manifold', 'An edge is in more than two faces.'); end
    heOf = 2 * eOf - 1;  heOf(~isFirst) = heOf(~isFirst) + 1;
    if any(tails(~isFirst) ~= tips(firstK(eOf(~isFirst)))), error('rheome:flipOut:orientation', 'The faces are not consistently oriented.'); end

    tail = zeros(nH, 1);  tail(2*(1:nE)-1) = tails(firstK);  tail(2*(1:nE)) = tips(firstK);
    face = zeros(nH, 1);  face(heOf) = fOf;
    nxt = zeros(nH, 1);  kNext = (1:K).' + 1;  kNext(3:3:end) = kNext(3:3:end) - 3;
    nxt(heOf) = heOf(kNext);
    fHe = heOf(1:3:end);
    [~, lastK] = unique(tails, 'last');                            % v.halfedge(): the last one assigned
    vHe = heOf(lastK);

    % == Resolve boundary loops (walked as geometry-central does; v.halfedge() of a boundary vertex is its interior
    %    boundary halfedge)
    tw = @(h) h + 1 - 2 * (mod(h, 2) == 0);
    nB = 0;
    for iHe = 1:nH
        if face(iHe) ~= 0, continue, end
        loopInd = nF + nB + 1;  fHe(loopInd, 1) = iHe;  nB = nB + 1;
        curr = iHe;
        while true
            face(curr) = loopInd;
            currT = tw(curr);  vHe(tail(currT)) = currT;
            prev = curr;
            curr = tw(nxt(tw(curr)));
            while face(curr) ~= 0
                if curr == iHe, break, end
                curr = tw(nxt(curr));
            end
            nxt(curr) = prev;
            if curr == iHe, break, end
        end
    end
    isInt = face <= nF;
    vBnd = ~isInt(tw(vHe));

    M = struct('nV', nV, 'nF', nF, 'nE', nE, 'nH', nH, 'nB', nB, 'next', nxt, 'tail', tail, 'face', face, ...
        'vHe', vHe, 'fHe', fHe, 'isInt', isInt, 'vBnd', vBnd, 'P', V);

    % == The input geometry (VertexPositionGeometry)
    hI = find(isInt);                                              % the corners, in halfedge order
    eH = 2 * (1:nE).' - 1;
    d = V(tail(eH), :) - V(tail(nxt(eH)), :);
    L = sqrt(d(:,1).*d(:,1) + d(:,2).*d(:,2) + d(:,3).*d(:,3));   % EmbeddedGeometryInterface::computeEdgeLengths
    M.L = L;
    edgeOf = @(h) ceil(h / 2);                                     % the edge of halfedge h

    % EmbeddedGeometryInterface::computeCornerAngles: the angle at B = tail(c), between C = tip(c) and A = tail(prev(c))
    pB = V(tail(hI), :);  pC = V(tail(nxt(hI)), :);  pA = V(tail(nxt(nxt(hI))), :);
    u = unit3(pC - pB);  w = unit3(pA - pB);
    q = clampv(u(:,1).*w(:,1) + u(:,2).*w(:,2) + u(:,3).*w(:,3), -1, 1);
    corner = nan(nH, 1);  corner(hI) = acos(q);

    % IntrinsicGeometryInterface::computeVertexAngleSums: corners added in halfedge order
    M.angleSum = seqsum(tail(hI), corner(hI), nV);

    % EmbeddedGeometryInterface::computeFaceAreas: N = sum of cross(pA, pB) around the face from f.halfedge()
    h1 = fHe(1:nF);  h2 = nxt(h1);  h3 = nxt(h2);
    N = cross3(V(tail(h1),:), V(tail(h2),:));
    N = N + cross3(V(tail(h2),:), V(tail(h3),:));
    N = N + cross3(V(tail(h3),:), V(tail(h1),:));
    area = 0.5 * sqrt(N(:,1).*N(:,1) + N(:,2).*N(:,2) + N(:,3).*N(:,3));

    % IntrinsicGeometryInterface::computeHalfedgeVectorsInFace
    lAB = L(edgeOf(h1));  lBC = L(edgeOf(h2));  lCA = L(edgeOf(h3));
    hgt = 2 .* area ./ lAB;
    wd = sqrt(max(0, lCA.*lCA - hgt.*hgt));
    flip = lBC.*lBC > (lAB.*lAB + lCA.*lCA);  wd(flip) = -wd(flip);
    hvF = nan(nH, 2);
    hvF(h1, :) = [lAB, zeros(nF, 1)];
    hvF(h2, :) = [wd - lAB, hgt - 0];
    hvF(h3, :) = [-wd, -hgt];
    M.hvF = hvF;

    % IntrinsicGeometryInterface::computeCornerScaledAngles and computeHalfedgeVectorsInVertex (orbit ccw from v.halfedge())
    scaled = nan(nH, 1);
    vt = tail(hI);  sFull = 2.0 * pi ./ M.angleSum(vt);  sB = pi ./ M.angleSum(vt);
    bt = vBnd(vt);  sFull(bt) = sB(bt);
    scaled(hI) = sFull .* corner(hI);
    hvV = nan(nH, 2);
    cur = vHe;  coordSum = zeros(nV, 1);  live = true(nV, 1);
    while any(live)
        c = cur(live);
        hvV(c, :) = [cos(coordSum(live)) .* L(edgeOf(c)), sin(coordSum(live)) .* L(edgeOf(c))];
        coordSum(live) = coordSum(live) + scaled(c);
        stop = ~isInt(c);
        nc = c;  nc(~stop) = tw(nxt(nxt(c(~stop))));
        cur(live) = nc;
        idx = find(live);  live(idx(stop | nc == vHe(idx))) = false;
    end
    M.hvV = hvV;

    % == The intrinsic triangulation's start: IntrinsicGeometryInterface::computeCornerAngles (from lengths) summed in
    %    halfedge order (its vertexAngleSums), and the signpost angles of SignpostIntrinsicTriangulation's constructor
    cornerI = nan(nH, 1);  cornerI(hI) = lawOfCosines(L, hI, nxt, edgeOf);
    M.angleSumI = seqsum(tail(hI), cornerI(hI), nV);
    signpost = zeros(nH, 1);
    cur = vHe;  running = zeros(nV, 1);  live = true(nV, 1);
    while any(live)
        c = cur(live);
        signpost(c) = running(live);
        stop = ~isInt(c);
        r = running(live);  r(~stop) = r(~stop) + cornerI(c(~stop));  running(live) = r;
        nc = c;  nc(~stop) = tw(nxt(nxt(c(~stop))));
        cur(live) = nc;
        idx = find(live);  live(idx(stop | nc == vHe(idx))) = false;
    end
    M.signpost = signpost;
end

function a = lawOfCosines(L, h, nxt, edgeOf)
% EdgeLengthGeometry::cornerAngle (edge_length_geometry.ipp) at corner h
    lA = L(edgeOf(h));  lOpp = L(edgeOf(nxt(h)));  lB = L(edgeOf(nxt(nxt(h))));
    q = (lA .* lA + lB .* lB - lOpp .* lOpp) ./ (2. .* lA .* lB);
    a = acos(clampv(q, -1, 1));
end

function s = seqsum(g, x, n)
% sum x by group g, adding in the order given (C++'s += loop), 0 for an empty group
    [gs, o] = sort(g);  xs = x(o);                                 % sort is stable: the order within a group is kept
    first = [true; diff(gs) > 0];  start = find(first);
    grp = cumsum(first);  k = (1:numel(gs)).' - start(grp) + 1;
    T = zeros(n, max(k));  T(sub2ind(size(T), gs, k)) = xs;
    s = T(:, 1);
    for j = 2:size(T, 2), s = s + T(:, j); end
end

function u = unit3(v)
    r = 1 ./ sqrt(v(:,1).*v(:,1) + v(:,2).*v(:,2) + v(:,3).*v(:,3));
    u = v .* r;
end

function c = cross3(u, v)
    c = [u(:,2).*v(:,3) - u(:,3).*v(:,2), u(:,3).*v(:,1) - u(:,1).*v(:,3), u(:,1).*v(:,2) - u(:,2).*v(:,1)];
end

function x = clampv(x, lo, hi)
% geometry-central's clamp (NaN passes through)
    x(x > hi) = hi;  x(x < lo) = lo;
end
