function [velocityField, info] = opticalflow_connection(activityMap, surface, opts)
% DYNAMICS.OPTICALFLOW_CONNECTION  Horn-Schunck optical flow smoothed by the CONNECTION LAPLACIAN (the Cortical Flow Editor's method).
%
%   [velocityField, info] = rheome.dynamics.opticalflow_connection(activityMap, surface [,opts])
%
% A port of the Editor's toolbox/flow/opticalFlow.ts (with differential/intrinsicGauge.ts) and of nxr-neuro's
% opticalflow_connection.cpp: the same mathematics, the same numbers. Per frame t, with I the map normalised to unit peak,
% Idot = I(t+1) - I(t), grad I the vertex gradient and v a tangent vector at a vertex, written z = (a, b) in C's gauge:
%
%     min over v   sum_v M_v (Idot + v . grad I)^2  +  alpha  v' C v
%
% C is the connection Laplacian of the intrinsic Delaunay triangulation (rheome.operators.connection_laplacian_idt: cotan weights,
% each edge carrying a vector across by its parallel transport), M the mixed-Voronoi mass. v'Cv compares a vector with its
% neighbour's TRANSPORTED to it, so the smoothness does not depend on any frame. (The earlier penalty alpha (a'Ka + b'Kb) on the
% coefficients in per-vertex (1,1,1)xn frames did, and was removed 2026-11-11.)
%
% The system is in C's own gauge: vectors go in and out FACE BY FACE through the intrinsic corner angles (the signposts of the
% INPUT faces, before the Delaunay flips), never through a vertex's projected plane. The Editor's DATA TERM is reproduced: the
% gradient of the input faces (positions rounded to single and the face gradients stored in single, as the Editor holds them),
% each face's gradient placed at its corners' intrinsic angles and averaged by area; the data blocks are vertex-lumped
% (M_v g g'); a ridge of 1e-9 * alpha * mean(diag C) keeps a patch with a parallel field definite.
%
% INPUTS:
%   activityMap  [nV x nT] non-negative scalar activity per vertex per frame
%   surface      struct with .Vertices [nV x 3] and .Faces [nF x 3] (in metres: the velocity is in surface units per frame)
%   opts.alpha   smoothness weight, dimensionless (default 0.01, the Editor's). ⚠ vCv is not zero on a rigid rotation of a sphere (curvature term), so the
%                prior pulls the speed LOW as alpha grows (sphere bump: 95 % of truth at 0.01, 80 % at 0.1, 33 % at 1).
%   opts.peak    the normaliser (default max of the map)
%   opts.mollify, opts.flipEps  as rheome.operators.intrinsic_delaunay (1e-6)
% OUTPUT:
%   velocityField [3nV x (nT-1)] ambient tangent velocity, rows [x1;y1;z1;x2;...] per frame, in surface units PER FRAME
%   info          .alpha .peak .mollify .flipEps .flips
%
% See also: rheome.dynamics.opticalflow_scalar, rheome.operators.connection_laplacian_idt, rheome.operators.intrinsic_delaunay
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'alpha') || isempty(opts.alpha), opts.alpha = 0.01; end
    if ~isfield(opts, 'mollify') || isempty(opts.mollify), opts.mollify = 1e-6; end
    if ~isfield(opts, 'flipEps') || isempty(opts.flipEps), opts.flipEps = 1e-6; end
    alpha = opts.alpha;
    if ~(isscalar(alpha) && alpha > 0 && isfinite(alpha))
        error('dynamics:opticalflow_connection:alpha', 'alpha = %g: the smoothness weight must be positive and finite.', alpha);
    end
    V = double(surface.Vertices);  F = double(surface.Faces);  n = size(V, 1);  nF = size(F, 1);
    I = double(activityMap);  nT = size(I, 2);
    if size(I, 1) ~= n, error('dynamics:opticalflow_connection:size', 'The map has %d vertices, the surface %d.', size(I,1), n); end
    if nT < 2, error('dynamics:opticalflow_connection:frames', 'The flow is between frames: it needs at least 2, got %d.', nT); end
    if any(~isfinite(I(:))) || any(I(:) < 0), error('dynamics:opticalflow_connection:map', 'The map must be finite and non-negative (an envelope or a magnitude).'); end
    peak = max(I(:));  if isfield(opts, 'peak') && ~isempty(opts.peak), peak = opts.peak; end
    info = struct('alpha', alpha, 'peak', peak, 'mollify', opts.mollify, 'flipEps', opts.flipEps, 'flips', NaN);
    velocityField = zeros(3*n, nT-1);
    if peak == 0, return, end
    I = I / peak;

    % ---- the pencil: C, M of the flipped (Delaunay) mesh; the corners of the INPUT faces (no flips: flipEps = Inf) ----
    [C, ~, Sd] = rheome.operators.connection_laplacian_idt(V, F, 1, mollify=opts.mollify, flipEps=opts.flipEps);
    info.flips = Sd.flips;
    M = Sd.mass;
    K = rheome.operators.intrinsic_delaunay(V, F, signposts=true, mollify=opts.mollify, flipEps=Inf);   % input halfedges
    cr = i_corners(V, F, K);
    Nrm = i_normals(V, F);
    g = i_face_gradient(V, F);                                   % [nF x 3], as the Editor holds them (single)
    ridge = 1e-9 * alpha * (full(trace(C)) / (2*n));
    aC = alpha * C;
    d1 = 2*(1:n)' - 1;  d2 = 2*(1:n)';
    for t = 1:nT-1
        gf = i_face_vectors(F, g, I(:, t));                        % [nF x 3] gradient of frame t
        z0 = i_gradient(cr, gf, n);                                % [n x 2] in C's gauge
        dI = I(:, t+1) - I(:, t);
        g1 = z0(:,1);  g2 = z0(:,2);
        rhs = zeros(2*n, 1);  rhs(d1) = -M.*g1.*dI;  rhs(d2) = -M.*g2.*dI;
        if ~any(rhs), continue, end
        D = sparse([d1; d1; d2; d2], [d1; d2; d1; d2], [M.*g1.^2 + ridge; M.*g1.*g2; M.*g1.*g2; M.*g2.^2 + ridge], 2*n, 2*n);
        z = (aC + D) \ rhs;
        velocityField(:, t) = reshape(i_embed(cr, Nrm, z(d1), z(d2))', [], 1);
    end
end

% ---------------------------------------------------------------------------------------------------------------------------
function cr = i_corners(V, F, K)
% The corners of the input faces (intrinsicGauge.ts: cornersOf). Halfedge h = 3(f-1)+k leaves vertex F(f,k).
    nF = size(F, 1);  nH = 3*nF;  n = size(V, 1);
    vt = K.F;                                                    % tail vertex of each halfedge
    nxt = @(h) h - mod(h-1,3) + mod(mod(h-1,3)+1,3);
    hs = (1:nH)';  h1 = nxt(hs);  h2 = nxt(h1);
    ln = K.len(K.edgeOf);                                        % mollified lengths
    a = ln(hs);  b = ln(h2);  o = ln(h1);
    cang = acos(max(-1, min(1, (a.^2 + b.^2 - o.^2) ./ (2*a.*b))));
    cone = 2*pi*ones(n, 1);  cone(K.boundary) = K.angleSum(K.boundary);
    sc = ones(n, 1);  ok = K.angleSum > 0;  sc(ok) = cone(ok) ./ K.angleSum(ok);
    cr.psi0 = cone(vt) .* K.angle ./ max(K.angleSum(vt), realmin);
    cr.scale = sc;  cr.width = cang .* sc(vt);  cr.vt = vt;
    T = reshape(F(:, [2 3 1])', [], 1);  P = reshape(F', [], 1);
    e = V(T, :) - V(P, :);  e = e ./ max(vecnorm(e, 2, 2), realmin);
    Pf = V(F(:,2),:) - V(F(:,1),:);  Qf = V(F(:,3),:) - V(F(:,1),:);
    nn = cross(Pf, Qf, 2);  l = vecnorm(nn, 2, 2);  cr.area = l/2;  nn = nn ./ max(l, realmin);
    nh = kron(nn, ones(3,1));
    cr.edge = e;  cr.perp = cross(nh, e, 2);                    % N x edge
    cr.face = kron((1:nF)', ones(3,1));
    cr.nF = nF;  cr.n = n;
    % corners of a vertex in fan order (by psi0), as one global key
    span = max(cr.psi0 + cr.width) + 1;  cr.span = span;
    [~, ord] = sortrows([vt, cr.psi0]);
    cr.ord = ord;  cr.key = vt(ord)*span + cr.psi0(ord);
    last = zeros(n, 1);  lastIdx = accumarray(vt(ord), (1:nH)', [n 1], @max);
    last(:) = lastIdx;  cr.lastIdx = last;
    ext = zeros(n, 1);   % the fan's extent: psi0 + width of its last corner
    hl = ord(max(lastIdx, 1));  ext(:) = cr.psi0(hl) + cr.width(hl);
    cr.wrap = ext > 2*pi - 1e-9;
    cr.in = lastIdx > 0;
end

function g = i_face_gradient(V, F)
% assembleFaceGradient: single positions, double arithmetic, single storage. Returns the three vertex gradients per face [nF x 9].
    P = double(single(V));
    A = P(F(:,1),:);  B = P(F(:,2),:);  Cc = P(F(:,3),:);
    nn = cross(B - A, Cc - A, 2);  twoA = vecnorm(nn, 2, 2);  u = nn ./ twoA;
    g = [cross(u, Cc - B, 2), cross(u, A - Cc, 2), cross(u, B - A, 2)] ./ twoA;
    g(~(twoA > 0), :) = 0;
    g = double(single(g));
end

function gf = i_face_vectors(F, g, s)
    s1 = s(F(:,1));  s2 = s(F(:,2));  s3 = s(F(:,3));
    gf = s1.*g(:,1:3) + s2.*g(:,4:6) + s3.*g(:,7:9);
end

function z = i_gradient(cr, gf, n)
% intrinsicGradient: each face's gradient at its corners' intrinsic angles, averaged by area.
    gh = gf(cr.face, :);
    gl = vecnorm(gh, 2, 2);
    th = atan2(sum(gh .* cr.perp, 2), sum(gh .* cr.edge, 2));
    psi = cr.psi0 + th .* cr.scale(cr.vt);
    A = cr.area(cr.face);  A(~(A > 0)) = 0;
    w = A .* gl;
    num1 = accumarray(cr.vt, w .* cos(psi), [n 1]);  num2 = accumarray(cr.vt, w .* sin(psi), [n 1]);
    ws = accumarray(cr.vt, A, [n 1]);
    z = zeros(n, 2);  k = ws > 0;
    z(k,1) = num1(k) ./ ws(k);  z(k,2) = num2(k) ./ ws(k);
end

function out = i_embed(cr, Nrm, a, b)
% embedAt: z = (a, b) in C's gauge, unscaled into the corner's face, projected on the tangent plane at length |z|.
    n = cr.n;  zl = hypot(a, b);  x = atan2(b, a);
    wrap = cr.wrap;
    x(wrap) = x(wrap) - 2*pi*floor(x(wrap) / (2*pi));
    q = (1:n)'*cr.span + x;
    q = max(q, (1:n)'*cr.span);
    k = discretize(q, [cr.key; Inf]);  k(isnan(k)) = 1;
    h = cr.ord(k);
    th = (x - cr.psi0(h)) ./ cr.scale;
    d = cos(th) .* cr.edge(h,:) + sin(th) .* cr.perp(h,:);
    d = d - sum(d .* Nrm, 2) .* Nrm;  pl = vecnorm(d, 2, 2);
    raw = cos(th) .* cr.edge(h,:) + sin(th) .* cr.perp(h,:);
    use = pl > 1e-6;  d(use,:) = d(use,:) ./ pl(use);  d(~use,:) = raw(~use,:);
    out = zl .* d;
    out(zl == 0 | ~cr.in, :) = 0;
end

function N = i_normals(V, F)
% tangentBasis's normals (stored in single): the sum of the (unnormalised) face normals at each vertex, normalised.
    n = size(V, 1);
    fn = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
    N = zeros(n, 3);
    for k = 1:3, for c = 1:3, N(:,c) = N(:,c) + accumarray(F(:,k), fn(:,c), [n 1]); end, end
    N = double(single(N ./ max(vecnorm(N, 2, 2), realmin)));   % the Editor holds its normals in single
end
