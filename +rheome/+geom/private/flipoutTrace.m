function pts = flipoutTrace(M, v, traceVec, maxIters, target)
% FLIPOUTTRACE  Trace a vector from a vertex straight across the input mesh, trimmed to end at a vertex (private).
%
%   pts = flipoutTrace(M, v, traceVec, maxIters, target)
%
% The trace of an intrinsic edge along the input surface: a port of geometry-central's traceGeodesic from a vertex
% (trace_geodesic.cpp: traceGeodesic, traceGeodesic_fromVertex, traceGeodesic_fromEdge, traceInFaceTowardsEdge,
% traceInFaceFromEdge, traceInFaceBarycentric, traceGeodesic_iterative; barycentric_coordinate_helpers.ipp) with
% TraceOptions{includePath = true, maxIters}, followed by the trimming in
% SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput (signpost_intrinsic_triangulation.cpp) and
% trimTraceResult (trace_geodesic.cpp). M is from flipoutMesh; v and target are vertices; traceVec is the vector in
% v's tangent space (geometry-central's polar convention).
%
% pts [n x 5]: one SurfacePoint per row, [type elem w1 w2 w3]: type 1 a vertex (elem), 2 a point on edge elem at
% tEdge = w1 (from the tail of e.halfedge()), 3 a point in face elem with barycentric (w1 w2 w3) in the order of
% f.halfedge(), f.halfedge().next(), ... (geometry-central's faceCoords).
%
% Author: Diellor Basha, 2026 (a port of geometry-central, MIT, Copyright (c) 2017-2019 Nicholas Sharp and the
% geometry-central contributors; licence: +rheome/+geom/LICENSE-geometry-central)

    pts = trace(M, v, traceVec, maxIters);
    % SignpostIntrinsicTriangulation::traceIntrinsicHalfedgeAlongInput(he, trimEnd = true)
    [trimmed, ok] = trimTraceResult(M, pts, target);
    if ok, pts = [trimmed; 1 target 0 0 0]; end                    % else: the untrimmed trace, as geometry-central keeps it
end

% ---------------------------------------------------------------------------------------------------------------
function pts = trace(M, v, traceVec, maxIters)
% port of traceGeodesic(geom, SurfacePoint startP, Vector2 traceVec, options) for a vertex startP (trace_geodesic.cpp)
    pts = [1 v 0 0 0];
    n2 = traceVec(1)*traceVec(1) + traceVec(2)*traceVec(2);
    if n2 == 0, return, end
    len = sqrt(n2);
    r = traceVec * (1. / len);                                     % unit(traceVec)
    sub = fromVertex(M, v, r, len);
    pts = iterative(M, pts, sub, maxIters);
end

function pts = iterative(M, pts, sub, maxIters)
% port of traceGeodesic_iterative
    iter = 0;
    while ~sub.terminated
        if iter >= maxIters
            pts(end+1, :) = edgePointFromHalfedge(sub.crossHe, sub.tCross); %#ok<AGROW>
            return
        end
        pts(end+1, :) = edgePointFromHalfedge(sub.crossHe, sub.tCross); %#ok<AGROW>
        sub = inFaceFromEdge(M, sub.crossHe, sub.tCross, sub.dir, sub.len);
        iter = iter + 1;
    end
    pts(end+1, :) = sub.endPoint;
end

function sub = fromVertex(M, v, dir, len)
% port of traceGeodesic_fromVertex
    wedgeHe = 0;  rel = [0 0];
    minCross = inf;  minHe = 0;  minDir = [0 0];
    first = M.vHe(v);  curr = first;
    while true
        if ~M.isInt(curr), break, end
        nextHe = twin(M.next(M.next(curr)));
        s = normalize2(M.hvV(curr, :));
        e = normalize2(M.hvV(nextHe, :));
        crossStart = s(1)*dir(2) - s(2)*dir(1);
        crossEnd = e(1)*dir(2) - e(2)*dir(1);
        if crossStart > 0. && crossEnd <= 0.
            wedgeHe = curr;  rel = cdiv(dir, s);
            break
        end
        if abs(crossStart) < minCross, minCross = abs(crossStart);  minHe = curr;  minDir = [1, 1e-12]; end
        if abs(crossEnd) < minCross, minCross = abs(crossEnd);  minHe = nextHe;  minDir = [1, -1e-12]; end
        curr = nextHe;
        if curr == first, break, end
    end
    if wedgeHe == 0                                                % no interval passed: trace along the closest halfedge
        if isCanonical(minHe), d = minDir; else, d = -minDir; end  % convertVecToEdge
        sub = fromEdge(M, ceil(minHe / 2), convertTToEdge(minHe, 0.), d, len);
        return
    end
    if M.vBnd(v), sumA = pi; else, sumA = 2. * pi; end
    rel = cpow(rel, M.angleSum(v) / sumA);
    rel = normalize2(rel);
    startDir = normalize2(M.hvF(wedgeHe, :));
    dirInFace = cmul(rel, startDir);
    sub = inFaceTowardsEdge(M, M.next(wedgeHe), dirInFace, len);
end

function sub = fromEdge(M, e, tEdge, dir, len)
% port of traceGeodesic_fromEdge
    tEdge = clamp(tEdge, 0., 1.);
    he = 2 * e - 1;
    if dir(2) >= 0.
        traceHe = twin(he);  hDir = -dir;  tEdge = 1.0 - tEdge;
    else
        traceHe = he;
        if ~M.isInt(twin(traceHe))
            sub = endAt([2 e tEdge 0 0], dir, 0);
            return
        end
        hDir = dir;
    end
    sub = inFaceFromEdge(M, traceHe, tEdge, hDir, len);
end

function sub = inFaceTowardsEdge(M, towardsHe, dir, len)
% port of traceInFaceTowardsEdge
    f = M.face(towardsHe);
    rootHe = M.next(M.next(towardsHe));
    vc = vertexCoords(M, f);
    vecBaryCanonical = cart2bary(vc, dir) * len;
    b = permuteFromCanonical(M, vecBaryCanonical, rootHe);
    b(1) = min(b(1), 1e-12);  b(2) = max(b(2), 0.);  b(3) = max(b(3), 0.);   % fmin / fmax (no NaN here)
    dif = -(b(1) + b(2) + b(3));
    if dif > 0
        b(2) = b(2) + dif / 2;  b(3) = b(3) + dif / 2;
    else
        b(1) = b(1) + dif;
    end
    iHe = indexInTriangle(M, rootHe);
    startPoint = [0 0 0];  startPoint(iHe + 1) = 1.0;
    fixed = permuteToCanonical(M, b, rootHe);
    hittable = [false false false];  hittable(mod(iHe + 1, 3) + 1) = true;
    sub = inFaceBarycentric(M, f, startPoint, fixed, dir, len, hittable);
end

function sub = inFaceFromEdge(M, fromHe, tCrossFrom, dirH, lenH)
% port of traceInFaceFromEdge
    faceHe = twin(fromHe);
    f = M.face(faceHe);
    vc = vertexCoords(M, f);
    v = -dirH;
    v(2) = max(v(2), 1e-9);
    heDir = normalize2(M.hvF(faceHe, :));
    dirInFace = cmul(heDir, v);
    b = permuteFromCanonical(M, cart2bary(vc, dirInFace), faceHe);
    b(3) = max(b(3), 1e-12);
    dif = -(b(1) + b(2) + b(3));
    if dif > 0
        b(3) = b(3) + dif;
    else
        b(1) = b(1) + dif / 3.;  b(2) = b(2) + dif / 3.;  b(3) = b(3) + dif / 3.;
    end
    tCrossFrom = clamp(tCrossFrom, 0., 1.);
    iHe = indexInTriangle(M, faceHe);
    startPoint = [0 0 0];
    startPoint(iHe + 1) = tCrossFrom;
    startPoint(mod(iHe + 1, 3) + 1) = 1.0 - tCrossFrom;
    fixed = permuteToCanonical(M, b, faceHe);
    hittable = [true true true];  hittable(iHe + 1) = false;
    sub = inFaceBarycentric(M, f, startPoint, fixed * lenH, dirInFace, lenH, hittable);
end

function sub = inFaceBarycentric(M, f, startPoint, vecBary, dir, len, hittable)
% port of traceInFaceBarycentric
    vc = vertexCoords(M, f);
    endPoint = startPoint + vecBary;
    if endPoint(1) >= 0 && endPoint(2) >= 0 && endPoint(3) >= 0
        sub = endAt([3 f endPoint], dir, 0);
        return
    end
    tRay = inf;  crossHe = 0;  iOpp = -777;
    currHe = M.fHe(f);
    for i = 0:2
        currHe = M.next(currHe);                                   % always opposite the i'th vertex
        tRaw = -startPoint(i+1) / vecBary(i+1);
        if ~hittable(mod(i + 1, 3) + 1) || vecBary(i+1) >= 0, continue, end
        if tRaw < tRay
            tRay = tRaw;  crossHe = currHe;  iOpp = i;
        end
    end
    tRay = clamp(tRay, 0., 1. - 1e-9);
    if crossHe == 0
        sub = endAt([3 f startPoint], dir, 0);
        return
    end
    onEdge = startPoint + tRay * vecBary;
    a = onEdge(mod(iOpp + 2, 3) + 1);  b = onEdge(mod(iOpp + 1, 3) + 1);
    tCross = clamp(a / (b + a), 0., 1.);
    lenRemaining = (1.0 - tRay) * len;
    crossingEdgeVec = vc(mod(iOpp + 2, 3) + 1, :) - vc(mod(iOpp + 1, 3) + 1, :);
    remainingDir = cdiv(dir, normalize2(crossingEdgeVec));
    if ~M.isInt(twin(crossHe))                                     % a boundary: stop on it
        sub = endAt(edgePointFromHalfedge(crossHe, tCross), remainingDir, lenRemaining);
        return
    end
    sub = struct('terminated', false, 'crossHe', crossHe, 'tCross', tCross, 'dir', remainingDir, ...
        'len', lenRemaining, 'endPoint', []);
end

function [pts, ok] = trimTraceResult(M, pts, target)
% port of trimTraceResult (trace_geodesic.cpp)
    while size(pts, 1) > 1
        b = pts(end, :);
        if b(1) == 2
            he = 2 * b(2) - 1;
            if M.tail(he) == target || M.tail(twin(he)) == target, pts(end, :) = []; continue, end
        end
        if b(1) == 3 || b(1) == 1, pts(end, :) = []; continue, end
        break
    end
    ok = false;
    if isempty(pts), return, end
    b = pts(end, :);
    switch b(1)
        case 1
            if b(2) == target, ok = true; return, end
            first = M.vHe(b(2));  h = first;                       % adjacentVertices
            while true
                if M.tail(twin(h)) == target, ok = true; return, end
                h = M.next(twin(h));
                if h == first, break, end
            end
        case 2
            he = 2 * b(2) - 1;
            ok = M.tail(he) == target || M.tail(twin(he)) == target || ...
                M.tail(M.next(M.next(he))) == target || M.tail(M.next(M.next(twin(he)))) == target;
        case 3
            h = M.fHe(b(2));
            ok = any([M.tail(h), M.tail(M.next(h)), M.tail(M.next(M.next(h)))] == target);
    end
end

% --- helpers (barycentric_coordinate_helpers.ipp, vector2.ipp, trace_geodesic.cpp) ------------------------------
function sub = endAt(p, dir, lenLeft)
    sub = struct('terminated', true, 'crossHe', 0, 'tCross', 0, 'dir', dir, 'len', lenLeft, 'endPoint', p);
end

function p = edgePointFromHalfedge(he, t)
% SurfacePoint(Halfedge he, double tHalfedge)
    p = [2 ceil(he / 2) convertTToEdge(he, t) 0 0];
end

function t = convertTToEdge(he, t)
    if ~isCanonical(he), t = 1.0 - t; end
end

function tf = isCanonical(he)
    tf = mod(he, 2) == 1;                                          % he == he.edge().halfedge()
end

function h = twin(h)
    h = h + 1 - 2 * (mod(h, 2) == 0);
end

function vc = vertexCoords(M, f)
% vertexCoordinatesInTriangle
    h = M.fHe(f);
    vc = [0 0; M.hvF(h, :); -M.hvF(M.next(M.next(h)), :)];
end

function r = cart2bary(vc, v)
% cartesianVectorToBarycentric (geometry-central solves with Eigen's colPivHouseholderQr)
    A = [vc(1,1) vc(2,1) vc(3,1); vc(1,2) vc(2,2) vc(3,2); 1 1 1];
    x = A \ [v(1); v(2); 0];
    r = x.';
    s = r(1) + r(2) + r(3);
    r = r - s / 3.;                                                % normalizeBarycentricDisplacement
end

function i = indexInTriangle(M, he)
% halfedgeIndexInTriangle (0-based)
    h = M.fHe(M.face(he));
    if he == h, i = 0; return, end
    h = M.next(h);
    if he == h, i = 1; return, end
    i = 2;
end

function r = permuteToCanonical(M, b, refHe)
    k = indexInTriangle(M, refHe);  r = [0 0 0];
    for i = 0:2, r(mod(i + k, 3) + 1) = b(i + 1); end
end

function r = permuteFromCanonical(M, b, refHe)
    k = indexInTriangle(M, refHe);  r = [0 0 0];
    for i = 0:2, r(i + 1) = b(mod(i + k, 3) + 1); end
end

function u = normalize2(v)
    r = 1. / sqrt(v(1)*v(1) + v(2)*v(2));
    u = [v(1)*r, v(2)*r];
end

function c = cmul(a, b)
% Vector2::operator*(Vector2): complex product
    c = [a(1)*b(1) - a(2)*b(2), a(1)*b(2) + a(2)*b(1)];
end

function c = cdiv(a, b)
% Vector2::operator/(Vector2): complex quotient (then / denom, as a product by 1/denom)
    denom = b(1)*b(1) + b(2)*b(2);
    r = 1. / denom;
    c = [(a(1)*b(1) + a(2)*b(2)) * r, (a(2)*b(1) - a(1)*b(2)) * r];
end

function c = cpow(a, p)
% Vector2::pow(double): std::pow(complex, double) = exp(p log(a)) (libc++)
    lr = log(hypot(a(1), a(2)));  th = atan2(a(2), a(1));
    re = p * lr;  im = p * th;
    m = exp(re);
    c = [m * cos(im), m * sin(im)];
end

function x = clamp(x, lo, hi)
    if x > hi, x = hi; elseif x < lo, x = lo; end
end
