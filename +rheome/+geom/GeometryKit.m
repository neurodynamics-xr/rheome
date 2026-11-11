classdef GeometryKit < handle
% GEOM.GEOMETRYKIT  The vector heat method's transport, and the gauge of the connection Laplacian, on a surface (nxr GeometryKit).
%
%   K = rheome.geom.GeometryKit(V, F)
%   w = K.lift(X, vertices)               % embedded tangent vectors [k x 3] at vertices -> complex numbers in each vertex's gauge
%   X = K.embedAtVertices(vertices, w)    % complex numbers in the vertices' gauges -> embedded tangent vectors [k x 3]
%   t = K.turn(v)                         % [N x 1] the unit rotation that carries a vector at any vertex of v's component to v
%   t = K.turnsTo(anchors)                % [N x 1] the rotation that carries a vector at each vertex to its anchor (0 = none)
%   w = K.carry(v, X, vertices)           % vectors at vertices carried to v's tangent plane: complex, ONE plane, add them
%   X = K.embed(v, w)                     % a complex number in v's gauge as a tangent vector at v
%
% Vectors at different places on a folded surface are not added in 3-D (a fold turns what is parallel along the surface): they are CARRIED
% to one place along the surface, where they lie in one plane, and only then added. The carry is the vector heat method's transport
% ((L + M/t) X = M delta e^{i psi}, Sharp, Soliman and Crane 2019), solved on the intrinsic Delaunay triangulation, per connected component
% (a hemisphere). A vector is measured by its angle in the connection Laplacian's gauge, face by face through the INPUT faces' corners
% (the Editor's intrinsicGauge.ts: toIntrinsic / embedAt); a vector with no angle (zero, or normal to the surface) carries 0.
%
% Ports nxr.geom.logmap.GeometryKit (lift, turn, carry, turns_to, embed_at_vertices, embed) and nxr.operators.intrinsic_gauge. Vertices are 1-based.
% Component-wise, the gauge is the one of rheome.operators.connection_laplacian (mesh='intrinsic'): e1 along the signposts' reference.
%
% See also: rheome.cores.Pooling, rheome.operators.connection_laplacian, rheome.dynamics.opticalflow_connection
%
% Author: Diellor Basha, 2026

    properties (SetAccess = private)
        N        % number of vertices
        Label    % [N x 1] component of each vertex (0: none)
        Local    % [N x 1] index of each vertex within its component
        Comp     % struct array: idx, Cor (corners), Nrm, Dec (decomposition of L + M/t), M
    end
    properties (Access = private)
        Cache    % containers.Map: anchor -> turn (the last 16)
        Order    % anchors cached, oldest first
    end

    methods
        function K = GeometryKit(V, F, opt)
            arguments, V (:,3) double, F (:,3) double, opt.mollify (1,1) double = 1e-6, opt.flipEps (1,1) double = 1e-6, end
            n = size(V, 1);  K.N = n;
            g = graph([F(:,1); F(:,2); F(:,3)], [F(:,2); F(:,3); F(:,1)], [], n);
            lab = conncomp(g)';  K.Label = lab;  K.Local = zeros(n, 1);
            K.Cache = containers.Map('KeyType', 'double', 'ValueType', 'any');  K.Order = [];
            comps = struct('idx', {}, 'Cor', {}, 'Nrm', {}, 'Dec', {}, 'M', {});
            for c = 1:max(lab)
                idx = find(lab == c);
                if numel(idx) < 4, K.Label(idx) = 0;  comps(c).idx = [];  continue, end
                loc = zeros(n, 1);  loc(idx) = 1:numel(idx);  K.Local(idx) = 1:numel(idx);
                Fc = loc(F(lab(F(:,1)) == c, :));  P = V(idx, :);
                comps(c) = rheome.geom.GeometryKit.component(idx, P, Fc, opt);
            end
            K.Comp = comps;
        end

        function psi = lift(K, X, vertices)
            vertices = vertices(:);  X = reshape(X, [], 3);  psi = nan(numel(vertices), 1);
            lab = K.Label(vertices);
            for c = unique(lab(lab > 0))'
                m = lab == c;  C = K.Comp(c);
                a = rheome.geom.GeometryKit.toIntrinsic(C.Cor, C.Nrm, X(m, :), K.Local(vertices(m)));
                w = vecnorm(X(m, :), 2, 2) .* exp(1i * a);  w(~isfinite(a)) = NaN;
                psi(m) = w;
            end
        end

        function X = embedAtVertices(K, vertices, w)
            vertices = vertices(:);  w = w(:);  X = zeros(numel(vertices), 3);
            lab = K.Label(vertices);
            for c = unique(lab(lab > 0))'
                m = lab == c;  C = K.Comp(c);
                X(m, :) = rheome.geom.GeometryKit.embedAt(C.Cor, C.Nrm, w(m), K.Local(vertices(m)));
            end
        end

        function X = embed(K, v, w)
            X = K.embedAtVertices(v, w);
        end

        function t = turn(K, v)
            if isKey(K.Cache, v), t = K.Cache(v); return, end
            c = K.Label(v);  if c == 0, error('GeometryKit:component', 'vertex %d is on a component too small to analyse', v); end
            C = K.Comp(c);  i = K.Local(v);
            b = zeros(numel(C.idx), 1);  b(i) = C.M(i);
            X = C.Dec \ b;  mag = abs(X);
            r = conj(X) ./ max(mag, realmin);  r(~(mag > 1e-12)) = 0;
            t = nan(K.N, 1);  t(C.idx) = r;
            K.Cache(v) = t;  K.Order(end+1) = v;
            if numel(K.Order) > 16, remove(K.Cache, K.Order(1));  K.Order(1) = []; end
        end

        function w = carry(K, v, X, vertices)
            if nargin < 4, vertices = find(K.Label == K.Label(v)); end
            vertices = vertices(:);
            if any(K.Label(vertices) ~= K.Label(v)), error('GeometryKit:carry', 'the vertices are not all on the anchor''s component'); end
            t = K.turn(v);  w = K.lift(X, vertices) .* t(vertices);
            w(~isfinite(w)) = 0;
        end

        function t = turnsTo(K, anchors)
            % one transport solve a distinct anchor, the anchors of a component solved together in blocks
            anchors = anchors(:);  t = zeros(K.N, 1);
            for c = 1:numel(K.Comp)
                C = K.Comp(c);  if isempty(C.idx), continue, end
                a = unique(anchors(anchors > 0 & K.Label(max(anchors, 1)) == c));
                for s = 1:64:numel(a)
                    blk = a(s : min(s + 63, numel(a)));
                    B = zeros(numel(C.idx), numel(blk));
                    B(sub2ind(size(B), K.Local(blk), (1:numel(blk))')) = C.M(K.Local(blk));
                    X = C.Dec \ B;  mag = abs(X);
                    T = conj(X) ./ max(mag, realmin);  T(~(mag > 1e-12)) = 0;
                    for j = 1:numel(blk)
                        m = find(anchors == blk(j) & K.Label == c);
                        t(m) = T(K.Local(m), j);
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function C = component(idx, P, F, opt)
            n = size(P, 1);
            S = rheome.operators.intrinsic_delaunay(P, F, signposts=true, mollify=opt.mollify, flipEps=opt.flipEps);
            Cl = rheome.operators.connection_laplacian(P, F, mollify=opt.mollify, flipEps=opt.flipEps);
            Kn = rheome.operators.intrinsic_delaunay(P, F, signposts=true, mollify=opt.mollify, flipEps=Inf);     % the input halfedges, unflipped
            M = S.mass;  t = mean(S.len)^2;
            C.idx = idx;  C.M = M;  C.Nrm = Cl.normal;
            C.Dec = decomposition(Cl.A + spdiags(M / t, 0, n, n), 'lu');
            C.Cor = rheome.geom.GeometryKit.corners(P, F, Kn);
        end

        function cr = corners(P, F, Kn)
            % the corners of the INPUT faces, before the Delaunay flips (intrinsicGauge.ts cornersOf)
            nF = size(F, 1);  nH = 3*nF;  n = size(P, 1);
            nxt = @(h) h - mod(h-1,3) + mod(mod(h-1,3)+1,3);
            hs = (1:nH)';  h1 = nxt(hs);  h2 = nxt(h1);
            vt = Kn.F;  ln = Kn.len(Kn.edgeOf);
            a = ln(hs);  b = ln(h2);  o = ln(h1);
            cang = acos(max(-1, min(1, (a.^2 + b.^2 - o.^2) ./ (2*a.*b))));
            cone = 2*pi*ones(n, 1);  cone(Kn.boundary) = Kn.angleSum(Kn.boundary);
            sc = ones(n, 1);  ok = Kn.angleSum > 0;  sc(ok) = cone(ok) ./ Kn.angleSum(ok);
            cr.psi0 = cone(vt) .* Kn.angle ./ max(Kn.angleSum(vt), realmin);
            cr.width = cang .* sc(vt);  cr.scale = sc;
            T = reshape(F(:, [2 3 1])', [], 1);  Q = reshape(F', [], 1);
            e = P(T, :) - P(Q, :);  e = e ./ max(vecnorm(e, 2, 2), realmin);
            nn = cross(P(F(:,2),:) - P(F(:,1),:), P(F(:,3),:) - P(F(:,1),:), 2);  l = vecnorm(nn, 2, 2);  nn = nn ./ max(l, realmin);
            cr.edge = e;  cr.perp = cross(kron(nn, ones(3,1)), e, 2);
            [~, ord] = sortrows([vt, cr.psi0]);                      % by vertex, then psi0 (stable)
            cr.order = ord;
            cnt = accumarray(vt, 1, [n 1]);  cr.ptr = [0; cumsum(cnt)];
            top = zeros(n, 1);  has = cnt > 0;
            pe = cr.ptr(2:end);  last = ord(pe(has));
            top(has) = cr.psi0(last) + cr.width(last);
            cr.span = top;
            cr.stride = 4 * max(max(top), 2*pi) + 1;
            cr.keys = (repelem((0:n-1)', cnt)) * cr.stride + cr.psi0(ord);
        end

        function psi = toIntrinsic(cr, nrm, X, v)
            % the intrinsic angle of embedded vectors X at (component-local) vertices v: the first corner whose projected wedge holds it; in a gap the nearest wedge edge
            k = numel(v);  psi = nan(k, 1);
            cnt = cr.ptr(v + 1) - cr.ptr(v);  live = cnt > 0;
            N = nrm(v, :);  d = X - sum(X .* N, 2) .* N;  dl = vecnorm(d, 2, 2);
            live = live & dl > 1e-6 * vecnorm(X, 2, 2);
            idx = find(live);  if isempty(idx), return, end
            reps = cnt(idx);  row = repelem((1:numel(idx))', reps);
            start = cumsum([0; reps(1:end-1)]);  within = (0:sum(reps)-1)' - repelem(start, reps);
            pv = cr.ptr(v(idx));  h = cr.order(pv(row) + within + 1);
            t1 = d(idx, :) ./ dl(idx);  t2 = cross(N(idx, :), t1, 2);  t1 = t1(row, :);  t2 = t2(row, :);
            s = cr.scale(v(idx));  s = s(row);  wth = cr.width(h) ./ s;
            e = cr.edge(h, :);  p = cr.perp(h, :);
            e1 = sum(e .* t1, 2);  e2 = sum(e .* t2, 2);  p1 = sum(p .* t1, 2);  p2 = sum(p .* t2, 2);
            det = e1 .* p2 - e2 .* p1;  okd = abs(det) > 1e-12;  sd = det;  sd(~okd) = 1;
            th = atan2(-e2 ./ sd, p2 ./ sd);
            inside = okd & th >= -1e-12 & th <= wth + 1e-12;
            psiIn = cr.psi0(h) + min(max(th, 0), wth) .* s;
            co = cos(wth);  si = sin(wth);  q1 = co .* e1 + si .* p1;  q2 = co .* e2 + si .* p2;
            m0 = abs(atan2(e2, e1));  m1 = abs(atan2(q2, q1));
            miss = reshape([m0 m1]', [], 1);  cand = reshape([cr.psi0(h) cr.psi0(h) + cr.width(h)]', [], 1);
            rows2 = repelem(row, 2);
            bestMiss = accumarray(rows2, miss, [numel(idx) 1], @min);
            sel = find(miss <= bestMiss(rows2));
            firstBest = accumarray(rows2(sel), sel, [numel(idx) 1], @min);
            psiGap = cand(firstBest);
            ins = find(inside);
            firstIn = accumarray(row(ins), ins, [numel(idx) 1], @min);
            has = firstIn > 0;
            r = psiGap;  r(has) = psiIn(firstIn(has));
            psi(idx) = r;
        end

        function out = embedAt(cr, nrm, w, v)
            % complex numbers w in the gauges of (component-local) vertices v as tangent vectors there, at length |w|
            k = numel(v);  zl = abs(w);  x = angle(w);
            wrap = cr.span(v) > 2*pi - 1e-9;
            x(wrap) = x(wrap) - 2*pi*floor(x(wrap) / (2*pi));
            q = (v - 1) * cr.stride + x;
            % searchsorted(keys, q, side='right') - 1: the number of keys <= q, via one stable sort (keys before queries on ties)
            nk = numel(cr.keys);
            [~, o] = sort([cr.keys; q]);
            isq = o > nk;  before = cumsum(isq);
            pos = zeros(k, 1);  pos(o(isq) - nk) = (find(isq) - before(isq));       % keys ahead of each query
            kk = pos;                                                 % 1-based index of the last key <= q (0: none)
            a = cr.ptr(v) + 1;  b = cr.ptr(v + 1);
            kk = min(max(kk, a), max(b, a));
            h = cr.order(min(max(kk, 1), nk));
            hasc = b >= a;
            th = (x - cr.psi0(h)) ./ cr.scale(v);
            d = cos(th) .* cr.edge(h, :) + sin(th) .* cr.perp(h, :);
            N = nrm(v, :);  pj = d - sum(d .* N, 2) .* N;  pl = vecnorm(pj, 2, 2);
            use = pl > 1e-6;  d(use, :) = pj(use, :) ./ pl(use);
            out = zl .* d;
            out(zl == 0 | ~hasc, :) = 0;
        end
    end
end
