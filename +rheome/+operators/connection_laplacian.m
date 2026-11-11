function C = connection_laplacian(V, F, opt)
% OPERATORS.CONNECTION_LAPLACIAN  Levi-Civita connection Laplacian on the tangent bundle.
%
%   C = rheome.operators.connection_laplacian(V, F)                                   % intrinsic Delaunay mesh (default)
%   C = rheome.operators.connection_laplacian(V, F, mesh='embedded', mass='galerkin')  % the pre-2026-11-10 operator
%   C = rheome.operators.connection_laplacian(V, F, mesh='intrinsic', mass='mixedVoronoi', mollify=1e-6, flipEps=1e-6)
%
% ONE connection Laplacian with a mesh option (PI, 2026-11-10: "one mesh and mass"). 'intrinsic' (default) is the
% construction of rheome.operators.connection_laplacian_idt (the Editor's connectionLaplacian.ts: signposts on the
% mollified intrinsic Delaunay triangulation, weights >= 0), returned as the complex n x n Hermitian A; 'embedded' is the
% surface's own triangles. On the intrinsic mesh .Rt lives on the intrinsic edges (flipped edges join vertices that are
% not neighbours in F), and .e1 is the unit tangent along the signposts' reference direction; the detectors and flow
% code that walk the faces of F with .Rt (detect.*, flow.*, operators.gauge / trivial_connection) ask for 'embedded'.
% The intrinsic form also works on a mesh with boundary; the embedded one needs a closed patch.
%
% Pure-MATLAB port of geometry-central's vertex connection Laplacian (the operator behind
% nxr-compute's 'laplacian','connection'). It lifts the scalar cotan Laplacian to the
% COMPLEX tangent bundle: a tangent vector at vertex v is a complex number in the vertex
% frame (e1,e2), and the Hermitian matrix A couples neighbours through the intrinsic
% Levi-Civita parallel transport. Its off-diagonal is  A_ij = -w_ij * rot_ij  (w = cotan
% weight, rot = unit-complex transport). The connection carries the surface's Gaussian
% curvature as holonomy, so winding numbers of a tangent field computed with it obey
% Poincare-Hopf exactly (sum of indices = Euler characteristic).
%
% INPUTS:
%   V [nV x 3] vertices,  F [nF x 3] triangles (a single connected, oriented, closed patch;
%   call it per hemisphere -- use rheome.utils.hemisphere)
%
% OUTPUT (struct C):
%   .A        [nV x nV] complex Hermitian connection Laplacian
%   .B        [nV x nV] real mass (rheome.operators.mass, type opt.mass; default mixedVoronoi)
%   .e1 .e2   [nV x 3]  per-vertex tangent frame (gauge); e2 = normal x e1
%   .normal   [nV x 3]  per-vertex unit normal (area-weighted)
%   ⚠ .e1/.e2 ARE NOT A USABLE GAUGE. They are built along each vertex's reference half-edge, an
%   arbitrary per-vertex choice with no relation between neighbours: median ambient angle between
%   neighbouring e1 is 79 deg and the transport residual is 1.35 rad (a reference subject left
%   hemisphere). The OPERATOR does not care, because .Rt corrects for it, and nothing here is
%   wrong. But any analysis that reads a per-COMPONENT value in this frame is reading the mesh's
%   half-edge ordering. ⭐ Call rheome.operators.gauge for a smooth or a parallel frame instead.
%
%   ⚠ THE CONNECTION IS NOT FLAT PER FACE. The charts use cornerScaledAngles -- angles rescaled so
%   each vertex sums to 2*pi -- which removes the angle defect from the VERTEX and pushes it into
%   the FACES: per-face holonomy has a non-integer part up to 0.95 rad. The bookkeeping is exact
%   and easy to misread: over the oriented face-edge incidence d1,
%       sum_f (d1*argR)_f = 0 identically,  sum_f frac = +chi turns,  sum_f round = -chi
%   so Gauss-Bonnet lives in the FRACTIONAL parts while the integer parts carry -chi, and the 4020
%   faces with a nonzero rounded defect are not 4020 singularities.
%
%   .Rt       [nV x nV] sparse unit-complex transport, Rt(i,j) = rot carrying a tangent
%                       vector from i to j (for the winding: see rheome.detect.criticalPoints)
%   .nV
%
% OPTIONS: mesh 'intrinsic' (default) | 'embedded';  mass 'mixedVoronoi' (default) | 'galerkin' | 'barycentric';
%   mollify, flipEps  the intrinsic mesh's relative mollification (1e-6) and flip tolerance (1e-6)
% Also in C: .mesh, .massType, .mollify, .flipEps, .mollifyEps, .flips (NaN for the embedded mesh).
%
% See also: rheome.operators.connection_laplacian_idt, rheome.operators.laplace_beltrami, rheome.operators.mass, rheome.utils.hemisphere, rheome.detect.criticalPoints
%
% Author: Diellor Basha, 2026 (formulas: geometry-central, Sharp/Crane)

    arguments
        V (:,3) double, F (:,3) double
        opt.mesh (1,1) string {mustBeMember(opt.mesh, ["intrinsic" "embedded"])} = "intrinsic"
        opt.mass (1,1) string {mustBeMember(opt.mass, ["mixedVoronoi" "galerkin" "barycentric"])} = "mixedVoronoi"
        opt.mollify (1,1) double = 1e-6
        opt.flipEps (1,1) double = 1e-6
    end
    if opt.mesh == "intrinsic", C = i_intrinsic(V, F, opt); return, end
    nV = size(V,1);  nF = size(F,1);

    % ---- half-edge structure (slot s: edge F(:,s) -> F(:,s+1)); he index = (s-1)*nF + f ----
    S1 = [1 2 3];  S2 = [2 3 1];
    tail = reshape(F(:, S1), [], 1);          % [3nF x 1]
    tip  = reshape(F(:, S2), [], 1);
    nHE  = 3*nF;
    slotNext = [2 3 1];                        % next slot within a face
    heNext = ( (slotNext-1)' * nF + (1:nF) )'; heNext = heNext(:);   % next(he)
    % twin: match (tail,tip) with (tip,tail)
    key    = tail + (tip-1)*nV;                % undirected pair encoded per direction
    keyRev = tip  + (tail-1)*nV;
    [~, locRev] = ismember(keyRev, key);
    twin = locRev;                             % twin(he)
    if any(twin==0), error('operators:connection_laplacian:openMesh', ...
            'Mesh has boundary/non-manifold edges; pass a closed hemisphere.'); end

    % ---- corner angle at the TAIL of each half-edge (angle at F(f,s) in face f) ----
    p_i = V(tail,:);  p_j = V(tip,:);  p_k = V(tip(heNext),:);   % third vertex = tip of next he
    eij = p_j - p_i;  eik = p_k - p_i;
    cang = acos( max(-1,min(1, sum(eij.*eik,2)./(vecnorm(eij,2,2).*vecnorm(eik,2,2)) )) );  % corner angle
    angSum = accumarray(tail, cang, [nV 1]);                    % sum of corner angles per vertex
    cScaled = cang .* (2*pi ./ angSum(tail));                   % cornerScaledAngles

    % ---- CCW orbit around each vertex -> cumulative angle (coordSum) per half-edge ----
    nextOut = twin(heNext(heNext));            % next outgoing he CCW: twin(next(next(he)))
    coordSum = nan(nHE,1);
    seen = false(nV,1);
    for s = 1:nHE                              % start a cycle at the first unseen vertex's he
        v = tail(s);
        if seen(v), continue; end
        seen(v) = true;
        he = s;  acc = 0;
        while true
            coordSum(he) = acc;
            acc = acc + cScaled(he);
            he = nextOut(he);
            if he == s, break; end
        end
    end

    % ---- transport per half-edge: rot = exp(i*(coordSum[twin] - coordSum[he] + pi)) ----
    rot = exp(1i * (coordSum(twin) - coordSum + pi));           % transport i->j along he

    % ---- signed cotan edge weight (half-weight per corner opposite the edge) ----
    cotOpp = 1 ./ tan(cang);                   % cot of angle at tail; opposite edge = (tip, tip_next)
    % the corner at 'tail' is opposite edge (tip -> tip_next); accumulate 1/2 cot there
    eo_i = tip;  eo_j = tip(heNext);
    Wc = sparse([eo_i; eo_j], [eo_j; eo_i], 0.5*[cotOpp; cotOpp], nV, nV);
    w_he = full(Wc(tail + (tip-1)*nV));        % cotan weight of each half-edge's edge

    % ---- assemble A: A(iTail,iTip) += -w * transport[twin(he)]  ;  A(iTail,iTail) += w ----
    A = sparse([tail; tail], [tail; tip], [w_he; -w_he.*rot(twin)], nV, nV);
    A = (A + A')/2;                            % Hermitian symmetrize (kills round-off)

    % transport matrix for the detector: Rt(i,j) = rot carrying i -> j
    Rt = sparse(tail, tip, rot, nV, nV);

    % ---- vertex frames: e1 along the reference half-edge (v.halfedge), e2 = n x e1 ----
    fn = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:), 2);   % 2A * face normal
    N  = i_vertex_normals(F, fn, nV);
    refHe = zeros(nV,1);                        % one outgoing he per vertex (the "v.halfedge")
    for s = 1:nHE, v = tail(s); if refHe(v)==0, refHe(v)=s; end, end
    d  = V(tip(refHe),:) - V(tail(refHe),:);    % reference edge direction
    e1 = d - sum(d.*N,2).*N;  e1 = e1 ./ max(vecnorm(e1,2,2), eps);
    e2 = cross(N, e1, 2);

    M = i_mass(V, F, opt);
    C = struct('A', A, 'B', (M+M')/2, 'e1', e1, 'e2', e2, 'normal', N, 'Rt', Rt, 'nV', nV, 'mesh', "embedded", 'massType', opt.mass, ...
               'mollify', NaN, 'flipEps', NaN, 'mollifyEps', NaN, 'flips', NaN);
end

function M = i_mass(V, F, opt)
    if opt.mass == "mixedVoronoi", M = rheome.operators.mass(V, F, 'mixedVoronoi', mollify=opt.mollify);
    else, M = rheome.operators.mass(V, F, opt.mass); end
end

function C = i_intrinsic(V, F, opt)
% The complex form of connection_laplacian_idt: vertex block (i,j) of the real 2n form [a -b; b a] is the complex a + ib.
    nV = size(V,1);
    [L, ~, S, T] = rheome.operators.connection_laplacian_idt(V, F, 1, mollify=opt.mollify, flipEps=opt.flipEps);
    A = L(1:2:end, 1:2:end) + 1i * L(2:2:end, 1:2:end);  A = (A + A') / 2;
    % transport i -> j along an intrinsic edge is the conjugate of the rotation in A(i,j) = -w e^{ir}
    Rt = sparse([T.i; T.j], [T.j; T.i], [exp(-1i*T.r); exp(1i*T.r)], nV, nV);
    fn = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:), 2);
    N  = i_vertex_normals(F, fn, nV);
    ref = S.ref;  e1 = zeros(nV,3);  ok = ref > 0;      % the signposts' reference direction, in the tangent plane
    d = V(max(ref,1),:) - V;  d(~ok,:) = 0;
    e1 = d - sum(d.*N,2).*N;  e1 = e1 ./ max(vecnorm(e1,2,2), eps);
    e2 = cross(N, e1, 2);
    M = i_mass(V, F, opt);
    C = struct('A', A, 'B', (M+M')/2, 'e1', e1, 'e2', e2, 'normal', N, 'Rt', Rt, 'nV', nV, 'mesh', "intrinsic", 'massType', opt.mass, ...
               'mollify', opt.mollify, 'flipEps', opt.flipEps, 'mollifyEps', S.mollifyEps, 'flips', S.flips);
end

% area-weighted unit vertex normals
function N = i_vertex_normals(F, fn, nV)
    N = zeros(nV,3);
    for k = 1:3, N = N + [accumarray(F(:,k), fn(:,1), [nV 1]), ...
                          accumarray(F(:,k), fn(:,2), [nV 1]), ...
                          accumarray(F(:,k), fn(:,3), [nV 1])]; end
    N = N ./ max(vecnorm(N,2,2), eps);
end

% Author: Diellor Basha, 2026
