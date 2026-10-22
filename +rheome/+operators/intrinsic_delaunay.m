function S = intrinsic_delaunay(V, F, opt)
% OPERATORS.INTRINSIC_DELAUNAY  Intrinsic Delaunay triangulation of a surface and its cotangent Laplacian.
%
%   S = rheome.operators.intrinsic_delaunay(V, F)
%       S.F      [3nF x 1]  vertex of each halfedge (face f, corner k -> index 3(f-1)+k) after the flips
%       S.twin   [3nF x 1]  twin halfedge, 0 on a boundary
%       S.edgeOf [3nF x 1]  edge of each halfedge;  S.edgeHe [nE x 1] one halfedge of each edge
%       S.len    [nE x 1]   intrinsic edge lengths;  S.mass [nV x 1] mixed-Voronoi lumped mass of the INPUT mesh
%       S.L      [nV x nV]  cotangent Laplacian of the intrinsic triangulation (PSD, weights >= 0)
%       S.flips, S.nonDelaunay
%
% Same vertices and metric as the input; edges are flipped intrinsically (flipped length laid out
% from the two triangles) until every cotangent weight is >= -eps. Port of nxr-cortical-flow
% intrinsicDelaunay.ts, itself geometry-central's flipToDelaunay (Sharp, Soliman & Crane 2019),
% Bobenko & Springborn 2007. The mass is the mixed Voronoi area of the input triangulation
% (Meyer, Desbrun, Schroeder & Barr 2003), recorded before the first flip. Edge and halfedge
% numbering is the cortical-flow one (so matrices are comparable index by index). Closed or with
% boundary (a boundary edge is never flipped). Signposts / tufted cover are not ported.
%
% Author: Diellor Basha, 2026
arguments, V (:,3) double, F (:,3) double, opt.eps (1,1) double = 1e-6, end
nV = size(V,1); nF = size(F,1); H = 3*nF; epsw = opt.eps;
Fh = reshape(F', [], 1);                                     % halfedge h = 3(f-1)+k tails
nxt = @(h) h - mod(h-1,3) + mod(mod(h-1,3)+1,3);
hs = (1:H)'; Fn = Fh(nxt(hs));
k1 = Fh*(nV+1) + Fn; k2 = Fn*(nV+1) + Fh;
assert(numel(unique(k1)) == H, 'intrinsic_delaunay: an edge appears twice in one orientation (not oriented/manifold)');
[tf, loc] = ismember(k2, k1); twin = zeros(H,1); twin(tf) = loc(tf);
primary = twin == 0 | hs < twin;
eid = cumsum(primary); edgeOf = zeros(H,1);
edgeOf(primary) = eid(primary); t = twin > 0 & ~primary; edgeOf(t) = edgeOf(twin(t));
edgeHe = find(primary); nE = numel(edgeHe);
len = vecnorm(V(Fn(edgeHe),:) - V(Fh(edgeHe),:), 2, 2);
mass = i_mass(len, edgeOf, nF, Fh, nV);
queue = (1:nE)'; inq = true(nE,1); head = 1; flips = 0;
while head <= numel(queue)
    e = queue(head); head = head + 1; inq(e) = false;
    if ~(wgt(e) < -epsw), continue, end
    if ~flip(e), continue, end
    flips = flips + 1;
    h = edgeHe(e); tw = twin(h);
    for x = [nxt(h) nxt(nxt(h)) nxt(tw) nxt(nxt(tw))]
        n = edgeOf(x); if ~inq(n), inq(n) = true; queue(end+1,1) = n; end %#ok<AGROW>
    end
end
ew = arrayfun(@wgt, (1:nE)');
Fn = Fh(nxt(hs)); i = Fh(edgeHe); j = Fn(edgeHe); keep = i ~= j;
w = ew(keep); i = i(keep); j = j(keep);
S.L = sparse([i;j;i;j], [j;i;i;j], [-w;-w;w;w], nV, nV);
S.F = Fh; S.twin = twin; S.edgeOf = edgeOf; S.edgeHe = edgeHe; S.len = len; S.mass = mass;
S.nV = nV; S.nF = nF; S.nE = nE; S.flips = flips; S.nonDelaunay = nnz(ew < -epsw);

    function w = wgt(e)
        hh = edgeHe(e); tt = twin(hh); w = hcot(hh);
        if tt > 0, w = w + hcot(tt); end
    end
    function c = hcot(hh)
        a1 = nxt(hh); a2 = nxt(a1);
        lij = len(edgeOf(hh)); ljk = len(edgeOf(a1)); lki = len(edgeOf(a2));
        c = (-lij^2 + ljk^2 + lki^2) / (4*tri(lij,ljk,lki)) / 2;
    end
    function ok = flip(e)
        ok = false; hh = edgeHe(e); tt = twin(hh);
        if tt == 0, return, end
        a1 = nxt(hh); a2 = nxt(a1); b1 = nxt(tt); b2 = nxt(b1);
        if a1 == tt || b1 == hh, return, end
        va = Fh(hh); vb = Fh(a1); vc = Fh(a2); vd = Fh(b2);
        L = @(x) len(edgeOf(x));
        nl = flen(len(e), L(a1), L(a2), L(b1), L(b2));
        if ~isfinite(nl) || nl <= 0, return, end
        f1 = hh - mod(hh-1,3); f2 = tt - mod(tt-1,3);
        old = [b2 a1 a2 b1];
        outer = [edgeOf(old) twin(old)];
        if any(outer(:,2) == hh | outer(:,2) == tt), return, end
        Fh(f1) = vd; Fh(f1+1) = vb; Fh(f1+2) = vc;
        Fh(f2) = vc; Fh(f2+1) = va; Fh(f2+2) = vd;
        slots = [f1 f1+1 f2 f2+1];
        for q = 1:4
            s = slots(q); tw0 = outer(q,2); mv = find(old == tw0, 1);
            if ~isempty(mv), tw2 = slots(mv); else, tw2 = tw0; end
            edgeOf(s) = outer(q,1); twin(s) = tw2; edgeHe(outer(q,1)) = s;
            if isempty(mv) && tw0 > 0, twin(tw0) = s; end
        end
        edgeOf(f1+2) = e; edgeOf(f2+2) = e; twin(f1+2) = f2+2; twin(f2+2) = f1+2; edgeHe(e) = f1+2;
        len(e) = nl; ok = true;
    end
end

function a = tri(x, y, z)
s = (x+y+z)/2; a = sqrt(max(0, s*(s-x)*(s-y)*(s-z)));
end

function l = flen(lab, lbc, lca, lad, ldb)
ang = @(x,y,o) acos(max(-1, min(1, (x^2 + y^2 - o^2)/(2*x*y))));
th = ang(lab,lca,lbc) + ang(lab,lad,ldb);
l = sqrt(max(0, lca^2 + lad^2 - 2*lca*lad*cos(th)));
end

function M = i_mass(len, edgeOf, nF, Fh, nV)
% mixed Voronoi area: non-obtuse -> Voronoi corners; obtuse -> A/2 at the obtuse corner, A/4 elsewhere
h = reshape(1:3*nF, 3, [])';  l = len(edgeOf(h));            % l(:,k) = edge corner k -> k+1
s = sum(l,2)/2; A = sqrt(max(0, s.*(s-l(:,1)).*(s-l(:,2)).*(s-l(:,3))));
c = zeros(nF,3);
for k = 1:3
    a = l(:,k); b = l(:,mod(k+1,3)+1); o = l(:,mod(k,3)+1);
    c(:,k) = (a.^2 + b.^2 - o.^2) ./ (2*a.*b);
end
ct = c ./ sqrt(max(1e-300, 1 - c.^2));
P = zeros(nF,3);
for k = 1:3
    e2 = mod(k+1,3)+1;
    P(:,k) = (l(:,k).^2 .* ct(:,mod(k+1,3)+1) + l(:,e2).^2 .* ct(:,mod(e2+1,3)+1)) / 8;
end
[~, ob] = max(c < 0, [], 2); isob = any(c < 0, 2);
for k = 1:3
    ob_k = isob & ob == k;  oth = isob & ob ~= k;
    P(ob_k,k) = A(ob_k)/2; P(oth,k) = A(oth)/4;
end
M = accumarray(Fh(h(:)), P(:), [nV 1]);
end
