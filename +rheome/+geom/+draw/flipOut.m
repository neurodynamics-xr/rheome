function [P, info] = flipOut(S, from, to, mesh)
% GEOM.DRAW.FLIPOUT  The FlipOut geodesic between two vertices: the edge-flip shortening of a Dijkstra path, on the faces.
%
%   [P, info] = rheome.geom.draw.flipOut(S, from, to)
%   [P, info] = rheome.geom.draw.flipOut(S, from, to, mesh)     % mesh = rheome.geom.draw.flipMesh(S), or info.mesh
%
% A locally shortest path (a geodesic) on the surface from vertex `from` to vertex `to`, found by FlipOut (Sharp &
% Crane 2020): start from the shortest path along the mesh's edges (Dijkstra), then repeatedly straighten the
% sharpest bend in it by flipping edges of an intrinsic triangulation, until no bend is less than pi on either side.
% The path is then traced back onto the input faces: each point is a vertex, a crossing of an edge, or (rarely) a
% point in a face. This is the line the viewer draws (nxr-compute's tracePathSurface, geometry-central in WASM),
% reproduced here to 1e-9 m: the same Dijkstra path, the same order of straightening, the same trace.
%
% A PORT of FlipOut from geometry-central. Pure MATLAB; no MEX.
%   Paper:  Nicholas Sharp and Keenan Crane, "You Can Find Geodesic Paths in Triangle Meshes by Just Flipping Edges",
%           ACM Transactions on Graphics 39(6), SIGGRAPH Asia 2020.
%   Code:   geometry-central (https://github.com/nmwsharp/geometry-central), commit
%           8579865362abd2185ef97d7ef0dae09338ab8ff1 (nxr-compute's deps/geometry-central). Ported:
%           src/surface/flip_geodesics.cpp (FlipEdgeNetwork, FlipEdgePath), signpost_intrinsic_triangulation.cpp,
%           intrinsic_triangulation.cpp (+ .ipp), intrinsic_geometry_interface.cpp, embedded_geometry_interface.cpp,
%           edge_length_geometry.ipp, mesh_graph_algorithms.cpp (shortestEdgePath), trace_geodesic.cpp,
%           surface_mesh.cpp (SurfaceMesh::flip), manifold_surface_mesh.cpp (construction from a face list),
%           barycentric_coordinate_helpers.ipp, surface_point.ipp, elementary_geometry.ipp; the call sequence of
%           nxr-compute's src/geodesic_path.cpp (tracePathSurface). Each ported function names its source below.
%           Helpers: private/flipoutMesh.m (mesh and geometry), private/flipoutTrace.m (traceGeodesic).
%   Licence of the ported code: MIT License, Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central
%           contributors. Its permission notice is in +rheome/+geom/+draw/LICENSE-geometry-central and applies to this file
%           and its two helpers:
%             Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
%             associated documentation files (the "Software"), to deal in the Software without restriction,
%             including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
%             and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
%             subject to the following conditions: The above copyright notice and this permission notice shall be
%             included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
%             WITHOUT WARRANTY OF ANY KIND (see the licence file for the full text).
%
% INPUTS:
%   S       surface (.Vertices [nV x 3], .Faces [nF x 3], 1-based, oriented, manifold; a boundary is allowed)
%   from, to   vertex indices (1-based)
%   mesh    (optional) rheome.geom.draw.flipMesh(S), or info.mesh from an earlier call on S (the halfedge mesh
%           and its geometry, built once)
%
% OUTPUTS:
%   P     [n x 3] the path's points, from `from` to `to`, on the input faces
%   info  .length    the path's length (mesh units): FlipEdgeNetwork::length(), the sum of its intrinsic edges
%         .vertices  [n x 3] 1-based: a vertex [v v v], an edge crossing [a b a], a face point [a b c]
%         .weights   [n x 3] so P(i,:) = weights(i,:) * Vertices(vertices(i,:),:): a vertex [1 0 0], an edge
%                    [1-t t 0], a face its barycentric coordinates
%         .faces     [n-1 x 1] the input face each segment lies in
%         .nFlips .nShortenIters (FlipEdgeNetwork's counters)  .mesh (pass it back for the next path)
%         .cache (deprecated, one release: the same as .mesh)
%
% See also: rheome.geom.draw.flipMesh, rheome.geom.draw.line, rheome.geom.edgegraph
%
% Author: Diellor Basha, 2026 (port; the algorithm and the code ported are geometry-central's, above)

    if nargin < 4 || isempty(mesh)
        M = flipoutMesh(S.Vertices, S.Faces);
    elseif ~isFlipMesh(mesh)
        error('rheome:flipOut:mesh', 'mesh must be rheome.geom.draw.flipMesh(S) (a %s was passed).', class(mesh));
    elseif mesh.nV ~= size(S.Vertices, 1) || mesh.nF ~= size(S.Faces, 1)
        error('rheome:flipOut:mesh', 'The mesh was built for another surface (%d vertices, %d faces; S has %d, %d).', ...
            mesh.nV, mesh.nF, size(S.Vertices, 1), size(S.Faces, 1));
    else
        M = mesh;
    end
    if ~(isscalar(from) && isscalar(to) && from >= 1 && to >= 1 && from <= M.nV && to <= M.nV && from == round(from) && to == round(to))
        error('rheome:flipOut:vertex', 'from and to are vertex indices in 1..%d.', M.nV);
    end
    from = double(from);  to = double(to);

    if from == to                                                  % tracePathSurface: the point itself
        P = M.P(from, :);
        info = struct('length', 0, 'vertices', [from from from], 'weights', [1 0 0], 'faces', zeros(0, 1), ...
            'nFlips', 0, 'nShortenIters', 0, 'mesh', M, 'cache', M);
        return
    end

    % ===== FlipEdgeNetwork::constructFromDijkstraPath (flip_geodesics.cpp)
    dijkstraPath = shortestEdgePath(M, from, to);
    if isempty(dijkstraPath)
        error('rheome:flipOut:unreachable', 'Vertex %d cannot reach vertex %d on this surface.', from, to);
    end

    % ===== The intrinsic triangulation (SignpostIntrinsicTriangulation over a copy of the mesh)
    nxt = M.next;  tl = M.tail;  fc = M.face;  vHe = M.vHe;  fHe = M.fHe;
    iL = M.L;  sAng = M.signpost;  sumI = M.angleSumI;  vBnd = M.vBnd;
    isOrig = true(M.nE, 1);
    TEST_EPS = 1e-6;                                               % IntrinsicTriangulation::triangleTestEPS
    EPS_ANGLE = 1e-5;                                              % FlipEdgeNetwork::EPS_ANGLE
    SHORTEST = 0;  LEFT = 1;  RIGHT = 2;                           % SegmentAngleType

    % ===== FlipEdgeNetwork state: segments (FlipEdgePath::pathHeInfo; ids 1, 2, ... in geometry-central's order,
    %       0 = INVALID_IND), the paths at each edge (pathsAtEdge: front = e.halfedge() side), the wedge queue
    cap = 256;
    segHe = zeros(cap, 1);  segPrev = zeros(cap, 1);  segNext = zeros(cap, 1);  segAlive = false(cap, 1);
    nSeg = 0;
    atEdge = cell(M.nE, 1);
    qA = zeros(cap, 1);  qT = zeros(cap, 1);  qS = zeros(cap, 1);  qN = 0;
    nFlips = 0;  nShortenIters = 0;

    % FlipEdgeNetwork::FlipEdgeNetwork -> FlipEdgePath::FlipEdgePath (open path: first tail ~= last tip)
    loopPrev = 0;
    for k = 1:numel(dijkstraPath)
        newHe = dijkstraPath(k);
        sid0 = newSegment(newHe, loopPrev, 0);
        pushOutsideSegment(newHe, sid0);
        if loopPrev ~= 0, segNext(loopPrev) = sid0; end
        addToWedgeAngleQueue(sid0);
        loopPrev = sid0;
    end

    iterativeShorten();

    % ===== FlipEdgeNetwork::getPathPolyline: each intrinsic halfedge traced along the input
    heList = getHalfedgeList();
    pts = zeros(0, 5);
    for k = 1:numel(heList)
        tr = traceIntrinsicHalfedgeAlongInput(heList(k));
        if ~isempty(pts) && onSameElement(pts(end, :), tr(1, :)), pts(end, :) = []; end
        pts = [pts; tr]; %#ok<AGROW>
    end

    % ===== tracePathSurface (nxr-compute): positions, vertex/weight placement, length
    [P, Vt, W] = toPositions(M, pts);
    lenSum = 0;
    for k = 1:numel(heList), lenSum = lenSum + iL(ceil(heList(k) / 2)); end   % FlipEdgeNetwork::length()
    info = struct('length', lenSum, 'vertices', Vt, 'weights', W, 'faces', segmentFaces(M, pts), ...
        'nFlips', nFlips, 'nShortenIters', nShortenIters, 'mesh', M, 'cache', M);

    % =================================================================================================================
    % Nested: the network (flip_geodesics.cpp) and the intrinsic triangulation (signpost / intrinsic triangulation)
    % =================================================================================================================

    function id = newSegment(he, prev, next)
    % FlipEdgeNetwork::getNextUniquePathSegmentInd + a pathHeInfo entry
        nSeg = nSeg + 1;  id = nSeg;
        if id > numel(segHe)
            segHe(2*id) = 0;  segPrev(2*id) = 0;  segNext(2*id) = 0;  segAlive(2*id) = false;
        end
        segHe(id) = he;  segPrev(id) = prev;  segNext(id) = next;  segAlive(id) = true;
    end

    function s = getOutsideSegment(he)
    % port of FlipEdgeNetwork::getOutsideSegment
        q = atEdge{ceil(he / 2)};
        if isempty(q), s = 0; elseif mod(he, 2) == 1, s = q(1); else, s = q(end); end
    end

    function popOutsideSegment(he)
    % port of FlipEdgeNetwork::popOutsideSegment
        e = ceil(he / 2);
        if mod(he, 2) == 1, atEdge{e}(1) = []; else, atEdge{e}(end) = []; end
    end

    function pushOutsideSegment(he, s)
    % port of FlipEdgeNetwork::pushOutsideSegment
        e = ceil(he / 2);
        if mod(he, 2) == 1, atEdge{e} = [s, atEdge{e}]; else, atEdge{e} = [atEdge{e}, s]; end
    end

    function [leftAngle, rightAngle] = measureSideAngles(hePrev, heNext)
    % port of FlipEdgeNetwork::measureSideAngles
        v = tl(heNext);
        s = sumI(v);
        angleIn = sAng(tw(hePrev));
        angleOut = sAng(heNext);
        isB = vBnd(v);
        if angleIn < angleOut
            rightAngle = angleOut - angleIn;
        elseif isB
            rightAngle = inf;
        else
            rightAngle = (s - angleIn) + angleOut;
        end
        if angleOut < angleIn
            leftAngle = angleIn - angleOut;
        elseif isB
            leftAngle = inf;
        else
            leftAngle = (s - angleOut) + angleIn;
        end
    end

    function a = minWedgeAngle(id)
    % port of FlipEdgeNetwork::minWedgeAngle(const FlipPathSegment&) (-> locallyShortestTestWithType's angle)
        if segPrev(id) == 0, a = pi; return, end
        [l, r] = measureSideAngles(segHe(segPrev(id)), segHe(id));
        a = min(l, r);                                             % std::fmin (no NaN here)
    end

    function addToWedgeAngleQueue(id)
    % port of FlipEdgeNetwork::addToWedgeAngleQueue (with locallyShortestTestWithBoth)
        if id == 0, return, end
        prevID = segPrev(id);
        if prevID == 0, return, end
        [l, r] = measureSideAngles(segHe(prevID), segHe(id));
        if l < r
            minA = l;  maxA = r;
            if l > (pi - EPS_ANGLE), minT = SHORTEST; else, minT = LEFT; end
            if r > (pi - EPS_ANGLE), maxT = SHORTEST; else, maxT = RIGHT; end
        else
            minA = r;  maxA = l;
            if r > (pi - EPS_ANGLE), minT = SHORTEST; else, minT = RIGHT; end
            if l > (pi - EPS_ANGLE), maxT = SHORTEST; else, maxT = LEFT; end
        end
        if minT == SHORTEST, return, end
        heapPush(minA, minT, id);
        if maxT == SHORTEST, return, end
        heapPush(maxA, maxT, id);
    end

    function addAllWedgesToAngleQueue()
    % port of FlipEdgeNetwork::addAllWedgesToAngleQueue (the queue is a heap: the insertion order does not matter)
        for sa = find(segAlive(1:nSeg)).'
            if segPrev(sa) ~= 0, addToWedgeAngleQueue(sa); end
        end
    end

    function tf = wedgeIsClear(id, type)
    % port of FlipEdgeNetwork::wedgeIsClear (straightenAroundMarkedVertices = true)
        prevID = segPrev(id);
        heNext = segHe(id);  hePrev = segHe(prevID);
        tf = false;
        if type == LEFT
            if getOutsideSegment(hePrev) ~= prevID, return, end
            if getOutsideSegment(heNext) ~= id, return, end
            heCurr = nxt(hePrev);
            while heCurr ~= heNext
                if ~isempty(atEdge{ceil(heCurr / 2)}), return, end
                heCurr = nxt(tw(heCurr));
            end
        else
            if getOutsideSegment(tw(hePrev)) ~= prevID, return, end
            if getOutsideSegment(tw(heNext)) ~= id, return, end
            heCurr = tw(nxt(nxt(tw(hePrev))));
            while heCurr ~= heNext
                if ~isempty(atEdge{ceil(heCurr / 2)}), return, end
                heCurr = tw(nxt(nxt(heCurr)));
            end
        end
        tf = true;
    end

    function locallyShortenAt(id, type)
    % port of FlipEdgeNetwork::locallyShortenAt
        prevID = segPrev(id);
        if prevID == 0, return, end
        heNext = segHe(id);  hePrev = segHe(prevID);
        if type == SHORTEST, return, end
        nShortenIters = nShortenIters + 1;
        if prevID == id
            error('rheome:flipOut:loop', 'A single-edge loop (FlipEdgeNetwork::processSingleEdgeLoop) is not ported: open paths only.');
        end
        initPathLength = iL(ceil(hePrev / 2)) + iL(ceil(heNext / 2));
        if type == LEFT
            sPrev = hePrev;  sNext = heNext;  reversed = false;
        else
            sPrev = tw(heNext);  sNext = tw(hePrev);  reversed = true;
        end
        sCurr = nxt(sPrev);  sPrevTwin = tw(sPrev);
        while sCurr ~= sNext                                       % flip until a shorter path exists
            if sCurr == sPrevTwin
                sCurr = nxt(tw(sCurr));
                continue
            end
            if flipEdgeIfPossible(ceil(sCurr / 2))
                nFlips = nFlips + 1;
                sCurr = tw(nxt(tw(sCurr)));                        % re-check the previous edge
            else
                sCurr = nxt(tw(sCurr));
            end
        end
        newPathLength = 0.;
        newPath = zeros(1, 0);
        sCurr = nxt(sPrev);
        while true
            newPath(end+1) = tw(nxt(sCurr)); %#ok<AGROW>
            newPathLength = newPathLength + iL(ceil(nxt(sCurr) / 2));
            if sCurr == sNext, break, end
            sCurr = nxt(tw(sCurr));
        end
        if newPathLength > initPathLength, return, end
        if reversed, newPath = tw(fliplr(newPath)); end
        replacePathSegment(id, type, newPath);
    end

    function replacePathSegment(nextID, type, newHalfedges)
    % port of FlipEdgePath::replacePathSegment
        heNext = segHe(nextID);  prevID = segPrev(nextID);  nextNextID = segNext(nextID);
        hePrev = segHe(prevID);  prevPrevID = segPrev(prevID);
        if type == LEFT
            popOutsideSegment(hePrev);  popOutsideSegment(heNext);
        else
            popOutsideSegment(tw(hePrev));  popOutsideSegment(tw(heNext));
        end
        segAlive(prevID) = false;  segAlive(nextID) = false;
        if prevPrevID == nextID
            error('rheome:flipOut:loop', 'A length-2 loop is not ported: open paths only.');
        end
        lp = prevPrevID;
        for h = newHalfedges
            nid = newSegment(h, lp, 0);
            if type == LEFT, pushOutsideSegment(tw(h), nid); else, pushOutsideSegment(h, nid); end
            if lp ~= 0, segNext(lp) = nid; end
            addToWedgeAngleQueue(nid);
            lp = nid;
        end
        if lp ~= 0, segNext(lp) = nextNextID; end
        if nextNextID ~= 0
            segPrev(nextNextID) = lp;
            addToWedgeAngleQueue(nextNextID);
        end
        if type == LEFT
            addToWedgeAngleQueue(getOutsideSegment(heNext));
            addToWedgeAngleQueue(getOutsideSegment(hePrev));
        else
            addToWedgeAngleQueue(getOutsideSegment(tw(heNext)));
            addToWedgeAngleQueue(getOutsideSegment(tw(hePrev)));
        end
    end

    function iterativeShorten()
    % port of FlipEdgeNetwork::iterativeShorten (no iteration or length limit)
        nIterations = 0;
        while qN > 0
            [minAngle, angleType, sid] = heapPop();
            if ~segAlive(sid), continue, end                       % stale: the segment is gone
            if minWedgeAngle(sid) ~= minAngle, continue, end       % stale: the angle has changed
            if ~wedgeIsClear(sid, angleType), continue, end
            locallyShortenAt(sid, angleType);
            nIterations = nIterations + 1;
            if mod(nIterations, 1000) == 0                         % purgeStaleQueueEntries
                qN = 0;
                addAllWedgesToAngleQueue();
            end
        end
    end

    function hes = getHalfedgeList()
    % port of FlipEdgePath::getHalfedgeList (open path)
        alive = find(segAlive(1:nSeg));
        last = alive(segNext(alive) == 0);  last = last(1);
        hes = zeros(1, 0);  w = last;
        while w ~= 0
            hes(end+1) = segHe(w); %#ok<AGROW>
            w = segPrev(w);
        end
        hes = fliplr(hes);
    end

    % ----- SignpostIntrinsicTriangulation (signpost_intrinsic_triangulation.cpp, intrinsic_triangulation.ipp)

    function flipped = flipEdgeIfPossible(e)
    % port of SignpostIntrinsicTriangulation::flipEdgeIfPossible
        flipped = false;
        he = 2 * e - 1;
        if ~(fc(tw(he)) <= M.nF), return, end                      % isFixed: a boundary edge (no marked edges)
        [p0, p1, p2, p3] = layoutDiamond(he);
        A1 = cross2(p1 - p0, p3 - p0);
        A2 = cross2(p3 - p2, p1 - p2);
        areaEPS = TEST_EPS * (A1 + A2);
        if A1 < areaEPS || A2 < areaEPS, return, end
        d = p1 - p3;  newLength = sqrt(d(1)*d(1) + d(2)*d(2));
        if ~isfinite(newLength), return, end
        if ~meshFlip(e), return, end
        iL(e) = newLength;
        updateAngleFromCWNeighor(he);
        updateAngleFromCWNeighor(tw(he));
        isOrig(e) = false;
        flipped = true;
    end

    function [p0, p1, p2, p3] = layoutDiamond(iHe)
    % port of IntrinsicTriangulation::layoutDiamond (intrinsic_triangulation.ipp)
        A1h = nxt(iHe);  A2h = nxt(A1h);  B0h = tw(iHe);  B1h = nxt(B0h);  B2h = nxt(B1h);
        l01 = iL(ceil(A1h / 2));  l12 = iL(ceil(A2h / 2));  l23 = iL(ceil(B1h / 2));
        l30 = iL(ceil(B2h / 2));  l02 = iL(ceil(iHe / 2));
        p3 = [0. 0.];  p0 = [l30 0.];
        p2 = layoutTriangleVertex(p3, p0, l02, l23);
        p1 = layoutTriangleVertex(p2, p0, l01, l12);
    end

    function ok = meshFlip(e)
    % port of SurfaceMesh::flip(e, preventSelfEdges = false) (surface_mesh.cpp), implicit twins
        ok = false;
        ha1 = 2 * e - 1;  ha2 = nxt(ha1);  ha3 = nxt(ha2);
        hb1 = 2 * e;  hb2 = nxt(hb1);  hb3 = nxt(hb2);
        if nxt(ha3) ~= ha1 || nxt(hb3) ~= hb1, return, end         % not a triangle
        if ha2 == hb1 || hb2 == ha1, return, end                   % incident on a degree-1 vertex
        va = tl(ha1);  vb = tl(hb1);  vc = tl(ha3);  vd = tl(hb3);
        fa = fc(ha1);  fb = fc(hb1);
        if vHe(va) == ha1, vHe(va) = hb2; end
        if vHe(vb) == hb1, vHe(vb) = ha2; end
        fHe(fa) = ha1;  fHe(fb) = hb1;
        nxt(ha1) = hb3;  nxt(hb3) = ha2;  nxt(ha2) = ha1;
        nxt(hb1) = ha3;  nxt(ha3) = hb2;  nxt(hb2) = hb1;
        tl(ha1) = vc;  tl(hb1) = vd;
        fc(ha3) = fb;  fc(hb3) = fa;
        ok = true;
    end

    function updateAngleFromCWNeighor(he)
    % port of SignpostIntrinsicTriangulation::updateAngleFromCWNeighor
        if ~(fc(he) <= M.nF), sAng(he) = sumI(tl(he)); return, end
        if ~(fc(tw(he)) <= M.nF), sAng(he) = 0.; return, end
        cwHe = nxt(tw(he));
        neighAngle = sAng(cwHe);
        cAngle = cornerAngle(cwHe);
        sAng(he) = standardizeAngle(tl(he), neighAngle + cAngle);
    end

    function a = cornerAngle(heA)
    % port of EdgeLengthGeometry::cornerAngle (edge_length_geometry.ipp), on the intrinsic lengths
        heOpp = nxt(heA);  heB = nxt(heOpp);
        lOpp = iL(ceil(heOpp / 2));  lA = iL(ceil(heA / 2));  lB = iL(ceil(heB / 2));
        q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
        if q > 1.0, q = 1.0; elseif q < -1.0, q = -1.0; end
        a = acos(q);
    end

    function a = standardizeAngle(v, a)
    % port of SignpostIntrinsicTriangulation::standardizeAngle: std::fmod(angle, vertexAngleSums[v]) (exact)
        if vBnd(v), return, end
        s = sumI(v);
        if a >= s
            if a < 2 * s, a = a - s; else, a = mod(a, s); end      % a - s is exact for s <= a < 2s (Sterbenz)
        end
    end

    function pts = traceIntrinsicHalfedgeAlongInput(he)
    % port of SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(he, trimEnd = true)
        vA = tl(he);  vB = tl(tw(he));                             % vertexLocations: every vertex is an input vertex
        if isOrig(ceil(he / 2))
            pts = [1 vA 0 0 0; 1 vB 0 0 0];
            return
        end
        % halfedgeVector(he): Vector2::fromAngle(signpostAngle * scaleFac) * edgeLength
        if vBnd(vA), scaling = sumI(vA) / pi; else, scaling = sumI(vA) / (2. * pi); end
        scaleFac = 1.0 / scaling;
        th = sAng(he) * scaleFac;  len = iL(ceil(he / 2));
        traceVec = [cos(th) * len, sin(th) * len];
        pts = flipoutTrace(M, vA, traceVec, M.nF * 10, vB);
    end

    % ----- the wedge queue: std::priority_queue<tuple<double, SegmentAngleType, FlipPathSegment>, greater<>>
    %       a min-heap on (angle, type, segment id), duplicates kept
    function heapPush(a, t, s)
        qN = qN + 1;
        if qN > numel(qA), qA(2*qN) = 0;  qT(2*qN) = 0;  qS(2*qN) = 0; end
        i = qN;
        while i > 1
            p = floor(i / 2);
            if lessKey(a, t, s, qA(p), qT(p), qS(p))
                qA(i) = qA(p);  qT(i) = qT(p);  qS(i) = qS(p);  i = p;
            else
                break
            end
        end
        qA(i) = a;  qT(i) = t;  qS(i) = s;
    end

    function [a, t, s] = heapPop()
        a = qA(1);  t = qT(1);  s = qS(1);
        la = qA(qN);  lt = qT(qN);  ls = qS(qN);
        qN = qN - 1;
        if qN == 0, return, end
        i = 1;
        while true
            c = 2 * i;
            if c > qN, break, end
            if c < qN && lessKey(qA(c+1), qT(c+1), qS(c+1), qA(c), qT(c), qS(c)), c = c + 1; end
            if lessKey(qA(c), qT(c), qS(c), la, lt, ls)
                qA(i) = qA(c);  qT(i) = qT(c);  qS(i) = qS(c);  i = c;
            else
                break
            end
        end
        qA(i) = la;  qT(i) = lt;  qS(i) = ls;
    end
end

% =====================================================================================================================
% Local functions
% =====================================================================================================================

function path = shortestEdgePath(M, startVert, endVert)
% port of shortestEdgePath (mesh_graph_algorithms.cpp): Dijkstra over the edges, a min-heap on (distance, halfedge)
    incoming = zeros(M.nV, 1);                                     % incomingHalfedge (0: not discovered)
    cap = 1024;  hd = zeros(cap, 1);  hh = zeros(cap, 1);  n = 0;
    enqueue(startVert, 0.);
    path = zeros(1, 0);
    while n > 0
        currDist = hd(1);  currIn = hh(1);
        pop();
        currVert = M.tail(twinh(currIn));
        if currVert == startVert || incoming(currVert) ~= 0, continue, end
        incoming(currVert) = currIn;
        if currVert == endVert
            w = currVert;
            while w ~= startVert
                h = incoming(w);
                path(end+1) = h; %#ok<AGROW>
                w = M.tail(h);
            end
            path = fliplr(path);
            return
        end
        enqueue(currVert, currDist);
    end

    function enqueue(v, dist)
        first = M.vHe(v);  he = first;
        while true                                                 % v.outgoingHalfedges(): he.twin().next()
            u = M.tail(twinh(he));
            if ~(u == startVert || incoming(u) ~= 0)
                push(dist + M.L(ceil(he / 2)), he);
            end
            he = M.next(twinh(he));
            if he == first, break, end
        end
    end

    function push(d, h)
        n = n + 1;
        if n > numel(hd), hd(2*n) = 0;  hh(2*n) = 0; end
        i = n;
        while i > 1
            p = floor(i / 2);
            if d < hd(p) || (d == hd(p) && h < hh(p))
                hd(i) = hd(p);  hh(i) = hh(p);  i = p;
            else
                break
            end
        end
        hd(i) = d;  hh(i) = h;
    end

    function pop()
        ld = hd(n);  lh = hh(n);  n = n - 1;
        if n == 0, return, end
        i = 1;
        while true
            c = 2 * i;
            if c > n, break, end
            if c < n && (hd(c+1) < hd(c) || (hd(c+1) == hd(c) && hh(c+1) < hh(c))), c = c + 1; end
            if hd(c) < ld || (hd(c) == ld && hh(c) < lh)
                hd(i) = hd(c);  hh(i) = hh(c);  i = c;
            else
                break
            end
        end
        hd(i) = ld;  hh(i) = lh;
    end
end

function tf = lessKey(a1, t1, s1, a2, t2, s2)
% std::tuple's < on (angle, type, segment)
    tf = a1 < a2 || (a1 == a2 && (t1 < t2 || (t1 == t2 && s1 < s2)));
end

function h = tw(h)
% Halfedge::twin() with implicit twins: 2e-1 <-> 2e
    h = h + 1 - 2 * (mod(h, 2) == 0);
end

function h = twinh(h)
    h = h + 1 - 2 * (mod(h, 2) == 0);
end

function c = cross2(u, v)
    c = u(1) * v(2) - u(2) * v(1);
end

function pC = layoutTriangleVertex(pA, pB, lBC, lCA)
% port of layoutTriangleVertex (elementary_geometry.ipp)
    d = pB - pA;  lAB = sqrt(d(1)*d(1) + d(2)*d(2));
    s = (lAB + lBC + lCA) / 2.0;
    tArea = sqrt(max(0., s * (s - lAB) * (s - lBC) * (s - lCA)));
    h = 2. * tArea / lAB;
    w = (lAB * lAB - lBC * lBC + lCA * lCA) / (2. * lAB);
    r = 1. / lAB;  vABn = [d(1) * r, d(2) * r];
    perp = [-vABn(2), vABn(1)];
    pC = (pA + [vABn(1) * w, vABn(2) * w]) + [perp(1) * h, perp(2) * h];
end

function tf = onSameElement(a, b)
% port of onSameElement (surface_point.cpp)
    tf = a(1) == b(1) && a(2) == b(2);
end

function [P, Vt, W] = toPositions(M, pts)
% nxr-compute's tracePathSurface: SurfacePoint::interpolate (surface_point.ipp) and the vertex/weight placement
    n = size(pts, 1);  P = zeros(n, 3);  Vt = zeros(n, 3);  W = zeros(n, 3);
    for i = 1:n
        p = pts(i, :);
        switch p(1)
            case 1
                v = p(2);  P(i, :) = M.P(v, :);  Vt(i, :) = [v v v];  W(i, 1) = 1.0;
            case 2
                he = 2 * p(2) - 1;  a = M.tail(he);  b = M.tail(he + 1);  t = p(3);
                P(i, :) = (1. - t) * M.P(a, :) + t * M.P(b, :);
                Vt(i, :) = [a b a];  W(i, :) = [1.0 - t, t, 0];
            case 3
                h = M.fHe(p(2));  a = M.tail(h);  b = M.tail(M.next(h));  c = M.tail(M.next(M.next(h)));
                P(i, :) = (p(3) * M.P(a, :) + p(4) * M.P(b, :)) + p(5) * M.P(c, :);
                Vt(i, :) = [a b c];  W(i, :) = p(3:5);
        end
    end
end

function f = segmentFaces(M, pts)
% the input face each segment lies in: a face holding both of its ends (the lowest-numbered when two do)
    n = size(pts, 1);  f = zeros(max(n - 1, 0), 1);
    cand = cell(n, 1);
    for i = 1:n
        p = pts(i, :);
        switch p(1)
            case 1
                first = M.vHe(p(2));  h = first;  c = [];
                while true
                    if M.face(h) <= M.nF, c(end+1) = M.face(h); end %#ok<AGROW>
                    h = M.next(twinh(h));
                    if h == first, break, end
                end
            case 2
                c = [M.face(2*p(2) - 1), M.face(2*p(2))];  c = c(c <= M.nF);
            case 3
                c = p(2);
        end
        cand{i} = c;
    end
    for i = 1:n-1
        s = intersect(cand{i}, cand{i+1});
        if isempty(s), f(i) = NaN; else, f(i) = s(1); end
    end
end
