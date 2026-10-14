function [distance, solver] = geodesic(S, sourceVertices, solver)
% GEOM.GEODESIC  Single-/multi-source geodesic distance by Crane's heat method.
%
%   [distance, solver] = rheome.geom.geodesic(S, sourceVertices)          % builds the solver
%   [distance, solver] = rheome.geom.geodesic(S, sourceVertices, solver)  % reuse cached factors
%
% Distance ON THE SURFACE (not the ambient embedding) from a set of source vertices, via the
% heat method (Crane, Weischedel & Wardetzky 2013): diffuse a heat spike briefly, march down
% the normalised heat gradient, and integrate. Follows geometry-central's
% HeatMethodDistanceSolver (robust): everything is built INTRINSICALLY from edge lengths, first
% mollified (Sharp & Crane 2020) so that zero-area triangles (fsaverage has two) still factor, then
% flipped to the intrinsic Delaunay triangulation (nonnegative cotan weights); the heat step uses the
% LUMPED mass, t = h^2 (h = mean edge length).
%
% ⭐ WHY THE REWRITE (2026-10-14). The earlier version saturated at ~46 mm on fsaverage6 pial (true
% max ~240 mm) and was 20-30 % off within 40 mm, and its Cholesky failed on fsaverage. Three causes:
%   1. it zeroed the direction wherever |grad u| < eps -- an ABSOLUTE threshold, so every face where
%      the heat had decayed below ~1e-16 (beyond ~50 edges) got no direction and the Poisson
%      solution went flat there. Heat far from the source is tiny but still exact in double;
%      only its DIRECTION is used, so it is normalised wherever it is nonzero;
%   2. it used the Galerkin (consistent) mass and the extrinsic cotan weights, which are negative on
%      obtuse triangles: the heat operator loses the M-matrix sign pattern, the far heat cancels to
%      roundoff, flips sign and the direction field turns. Lumped mass on the intrinsic Delaunay
%      triangulation keeps it an M-matrix and the heat positive;
%   3. it pinned the Poisson nullspace with 1e-10*M, which is not scale invariant: in metres the
%      shift is ~1e-16 and the factorisation fails. The shift is now relative to the spectrum.
% Measured against geometry-central (potpourri3d) on fsaverage5/6/fsaverage: tests/test_geodesic.m
% and the request's evidence file.
%
% The returned solver caches the two Cholesky factorisations; pass it back for repeated source
% queries on the same surface (e.g. rheome.geom.karcher_mean), which is where the cost lives.
%
% INPUTS:
%   S              surface struct with .Vertices [nV x 3], .Faces [nF x 3]
%   sourceVertices [k x 1] vertex indices of the source set (distance 0 there)
%   solver         (optional) a solver returned by a previous call on the same S
%
% OUTPUT:
%   distance [nV x 1] geodesic distance in mesh units (0 at the sources, increasing away);
%            Inf on vertices no surface path reaches (another connected component, e.g. the
%            other hemisphere of a two-hemisphere cortex)
%   solver   cached factorisations for reuse
%
% See also: rheome.geom.edgegraph, rheome.geom.karcher_mean, rheome.geom.heatdistance
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(solver)
        solver = i_build(double(S.Vertices), double(S.Faces));
    end

    sourceIndicator = zeros(solver.nV, 1);
    sourceIndicator(sourceVertices) = 1;
    heat = solver.heatOp \ sourceIndicator;                % 1. brief heat diffusion

    gx = solver.Gx*heat;  gy = solver.Gy*heat;            % 2. per-face heat gradient (face frame)
    gradNorm = hypot(gx, gy);
    gradNorm(gradNorm == 0) = 1;                           % no heat at all here: no direction
    divergence = solver.Gx'*(solver.area.*(-gx./gradNorm)) ...   % 3. integrated divergence (= G'AX)
               + solver.Gy'*(solver.area.*(-gy./gradNorm));
    distance = solver.poissonOp \ divergence;              % 4. Poisson recovers distance
    distance = distance - mean(distance(sourceVertices));  % sources at 0 (as geometry-central)
    distance(heat == 0) = Inf;                             % unreached: another connected component
end

function solver = i_build(V, F)
% Intrinsic operators from (mollified) edge lengths on the intrinsic Delaunay triangulation, each face
% laid out in its own 2-D frame.
    nV = size(V, 1);  nF = size(F, 1);
    pairs = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];          % edge opposite corner 1, 2, 3
    [edges, ~, id] = unique(sort(pairs, 2), 'rows');
    FE = reshape(id, nF, 3);                                % FE(f,c): the edge opposite corner c
    len = vecnorm(V(edges(:,1),:) - V(edges(:,2),:), 2, 2);
    % intrinsic mollification (Sharp & Crane 2020): lengthen every edge by the same epsilon so each
    % triangle inequality holds with margin delta = 1e-6 h; makes degenerate triangles well defined
    L = len(FE);
    delta = 1e-6*mean(len);
    epsilon = max([0; delta - (L(:,2) + L(:,3) - L(:,1)); delta - (L(:,3) + L(:,1) - L(:,2)); ...
                   delta - (L(:,1) + L(:,2) - L(:,3))]);
    len = len + epsilon;
    h = mean(len);
    [F, FE, len] = i_delaunay(F, FE, len);

    % layout: p1 = (0,0), p2 = (l3,0), p3 = (x3,y3); hat-function gradient = rot90(opposite edge)/2A
    l1 = len(FE(:,1));  l2 = len(FE(:,2));  l3 = len(FE(:,3));
    x3 = (l2.^2 + l3.^2 - l1.^2)./(2*l3);
    y3 = sqrt(max(l2.^2 - x3.^2, 0));
    area = l3.*y3/2;
    e1 = [x3 - l3, y3];  e2 = [-x3, -y3];  e3 = [l3, zeros(nF,1)];   % p3-p2, p1-p3, p2-p1
    rows = repmat((1:nF)', 3, 1);  cols = F(:);
    solver.Gx = sparse(rows, cols, -[e1(:,2); e2(:,2); e3(:,2)]./repmat(2*area,3,1), nF, nV);
    solver.Gy = sparse(rows, cols,  [e1(:,1); e2(:,1); e3(:,1)]./repmat(2*area,3,1), nF, nV);
    solver.area = area;  solver.nV = nV;

    A = spdiags(area, 0, nF, nF);
    stiffness = solver.Gx'*A*solver.Gx + solver.Gy'*A*solver.Gy;     % = the cotan Laplacian
    stiffness = (stiffness + stiffness')/2;
    lumped = accumarray(cols, repmat(area/3, 3, 1), [nV 1]);
    solver.heatOp = decomposition(spdiags(lumped,0,nV,nV) + h^2*stiffness, 'chol');
    % pin the constant nullspace (one per connected component) with a shift 1e-10 below the top of
    % the spectrum: scale invariant, and its bias on the slowest mode is ~1e-5 on a hemisphere
    shift = 1e-10*mean(diag(stiffness))/mean(lumped);
    solver.poissonOp = decomposition(stiffness + shift*spdiags(lumped,0,nV,nV), 'chol');
end

function [F, FE, len] = i_delaunay(F, FE, len)
% Intrinsic Delaunay triangulation by edge flips (Fisher et al. 2007): flip every interior edge whose
% two opposite angles sum past pi. Same vertices, new connectivity and lengths; the cotan weights are
% then nonnegative, which keeps the heat positive far from the source (pial meshes are full of
% obtuse triangles). Edges on a boundary, non-manifold or between inconsistently oriented faces stay.
    nE = numel(len);  nF = size(F, 1);
    EF = zeros(nE, 2);  cnt = zeros(nE, 1);
    for f = 1:nF
        for c = 1:3
            e = FE(f,c);  cnt(e) = cnt(e) + 1;
            if cnt(e) <= 2, EF(e, cnt(e)) = f; end
        end
    end
    interior = cnt == 2;
    stack = find(interior);  top = numel(stack);  queued = interior;
    stack(end+1:end+4*nE) = 0;
    while top > 0
        e = stack(top);  top = top - 1;  queued(e) = false;
        f = EF(e,1);  g = EF(e,2);
        ca = find(FE(f,:) == e, 1);  cd = find(FE(g,:) == e, 1);
        if i_cot(len(FE(f,:)), ca) + i_cot(len(FE(g,:)), cd) >= -1e-12, continue, end
        a = F(f,ca);  b = F(f,mod(ca,3)+1);  c = F(f,mod(ca+1,3)+1);  d = F(g,cd);
        if F(g,mod(cd,3)+1) ~= c || F(g,mod(cd+1,3)+1) ~= b || a == d, continue, end
        eca = FE(f,mod(ca,3)+1);  eab = FE(f,mod(ca+1,3)+1);
        ebd = FE(g,mod(cd,3)+1);  edc = FE(g,mod(cd+1,3)+1);
        % lay the quad out on edge (b,c): b = 0, c = (lbc,0), a above, d below
        lbc = len(e);
        xa = (len(eab)^2 - len(eca)^2 + lbc^2)/(2*lbc);  ya =  sqrt(max(len(eab)^2 - xa^2, 0));
        xd = (len(ebd)^2 - len(edc)^2 + lbc^2)/(2*lbc);  yd = -sqrt(max(len(ebd)^2 - xd^2, 0));
        len(e) = hypot(xa - xd, ya - yd);
        F(f,:) = [a b d];  FE(f,:) = [ebd e eab];
        F(g,:) = [a d c];  FE(g,:) = [edc eca e];
        EF(ebd, EF(ebd,:) == g) = f;  EF(eca, EF(eca,:) == f) = g;
        for n = [eab eca ebd edc]
            if interior(n) && ~queued(n), top = top + 1;  stack(top) = n;  queued(n) = true; end
        end
    end
end

function c = i_cot(l, i)
% cotangent of the corner opposite edge i of a triangle with edge lengths l
    s = sum(l)/2;  A = sqrt(max(s*(s-l(1))*(s-l(2))*(s-l(3)), 0));
    c = (sum(l.^2) - 2*l(i)^2)/(4*max(A, realmin));
end

% Author: Diellor Basha, 2026
