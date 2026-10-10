function g = gauge(V, F, varargin)
% OPERATORS.GAUGE  A per-vertex tangent FRAME on a surface: a reference direction field.
%
%   g = rheome.operators.gauge(V, F)                       % smooth frame by vector diffusion
%   g = rheome.operators.gauge(V, F, Method="trivial")     % exactly parallel frame, two +1 poles
%   g = rheome.operators.gauge(V, F, Method="smoothest")   % Knoppel 2013: lowest connection-Laplacian eigenvector
%   g = rheome.operators.gauge(V, F, ..., Connection=C)    % reuse rheome.operators.connection_laplacian(V,F)
%
% WHY THIS EXISTS. rheome.operators.connection_laplacian returns .e1/.e2 and calls them a gauge, but
% they are built along each vertex's reference half-edge -- an arbitrary choice with no relation
% between neighbours. The OPERATOR does not care, because its transport .Rt corrects for it. Any
% analysis that reads a per-component value does care, because the components then describe the
% mesh's half-edge ordering. Measured on a reference subject left hemisphere (10242 vertices), poles at
% two faces 70 mm apart:
%
%   frame                        ambient angle between   transport residual   singular faces
%                                    neighbours           per edge vs LC      (per-face winding)
%   half-edge (connection .e1)         79.1 deg              1.351 rad              3920
%   vector diffusion                   21.8 deg              0.087 rad                20
%   smoothest eigenvector (Knoppel)    22.3 deg              0.086 rad                20
%   trivial connection, poles +1 +1    20.2 deg              0.079 rad                 2  (at the poles)
%
% ⭐ WHICH TO USE. 'trivial' when you want to CHOOSE where the frame breaks: it is exactly parallel
% (residual 0 against its own connection) everywhere but the prescribed faces, and on the cortex it
% is also the smoothest of the three by both measures above. On a sphere with poles at the poles it
% IS the geographic frame (0.36 deg median from east). 'diffusion' when you need no choice: one
% sparse solve, singularities wherever the seed puts them (t = 0.01, 1 and 100 agree to 0.1 deg;
% 'smoothest', its t -> inf limit, lands 0.5 deg away).
%
% ⚠ AN EARLIER VERSION CALLED THE TRIVIAL FRAME "AMBIENTLY ROUGHEST" (55.4 deg, 0.80 rad from LC)
% AND THAT WAS A BUG, NOT GEOMETRY. It solved d1*x = 2*pi*k - d1*argE, but d1*argE carries integer
% chart parts n_f besides the curvature, so the frame's index was k_f - n_f: it wound on 2560 faces.
% See rheome.operators.trivial_connection, which wraps the holonomy (Crane et al. 2010) and is what the
% prescribed indices now mean.
%
% ⚠⚠ NO SURFACE-WIDE SMOOTH FRAME EXISTS. Poincare-Hopf forces the indices of a tangent field to
% sum to the Euler characteristic (2 for a closed hemisphere mesh, which is what rheome.load.bases
% caches). Every method here carries singularities of total index chi and the frame is meaningless
% in their one-ring. .singular lists them, MEASURED from the frame. ⭐ 'trivial' is the only method
% that lets you CHOOSE where they go -- Singular=[f1 f2], SingularCharge=[+1 +1] (the default).
%
% INPUTS (name-value):
%   Method         "diffusion" (default) | "smoothest" | "trivial"
%   Time           diffusion time, default 1 (the result is insensitive to it; see above)
%   Seed           diffusion seed vertex, default round(nV/2)
%   Connection     a precomputed rheome.operators.connection_laplacian(V,F), to avoid rebuilding
%   Singular       ('trivial') face indices to place singularities at, default the two faces of
%                  largest separation
%   SingularCharge ('trivial') integer indices, default +1 at each of two faces
%                  ⚠ THEY MUST SUM TO chi (2 for a closed hemisphere). The face holonomies sum to
%                  2*pi*chi, so any other total is inconsistent; this function errors rather than
%                  let the least-squares solve put the remainder on the pinned face. A vortex/
%                  antivortex pair is [+1 -1] PLUS the budget elsewhere, e.g. [+1 -1 +1 +1].
%
% OUTPUTS (struct g):
%   .e1 .e2 [nV x 3] the tangent frame, orthonormal, e2 = normal x e1
%   .normal [nV x 3] the vertex normals (from the connection)
%   .phase  [nV x 1] the frame angle in the connection's own charts (e1 = cos*C.e1 + sin*C.e2)
%   .singular        face indices where the frame's winding is nonzero, and .charge -- measured
%                    from the frame for every method (for 'trivial' it equals the prescription)
%   .residual        median |wrapped transport error| against the gauge's OWN connection
%                    (0 for 'trivial', by construction)
%   .residualLC      the same against the unmodified Levi-Civita transport -- the number that is
%                    comparable across methods.
%   .neighbourAngle  median ambient angle between neighbouring frame vectors, degrees
%   .method
%
% REFERENCES. Crane, Desbrun & Schroeder (2010) Trivial connections on discrete surfaces, Comput.
% Graph. Forum 29(5) -- 'trivial'. Knoppel, Crane, Pinkall & Schroeder (2013) Globally optimal
% direction fields, ACM TOG 32(4), doi 10.1145/2461912.2462005 (Zotero X8KPFRAC) -- 'smoothest'.
% Both are in geometry-central and nxr-compute (smooth_field.cpp, gauge.cpp).
%
% See also: rheome.operators.trivial_connection, rheome.operators.connection_laplacian, rheome.detect.criticalPoints,
%           rheome.flow.directionfield
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Method', "diffusion");
    p.addParameter('Time', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('Seed', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Connection', []);
    p.addParameter('Singular', [], @isnumeric);
    p.addParameter('SingularCharge', [], @isnumeric);
    p.parse(varargin{:});
    o = p.Results;  method = lower(string(o.Method));

    F = double(F);  nV = size(V,1);  nF = size(F,1);
    C = o.Connection;
    if isempty(C), C = rheome.operators.connection_laplacian(V, F); end
    wr = @(z) mod(z + pi, 2*pi) - pi;

    switch method
        case "diffusion"
            s = o.Seed;  if isempty(s), s = round(nV/2); end
            u0 = zeros(nV,1);  u0(s) = 1;
            ph = angle((C.B + o.Time*C.A) \ (C.B*u0));
        case "smoothest"
            [ue,~] = eigs(C.A, C.B, 1, 'smallestabs', struct('tol',1e-8));
            ph = angle(ue);        % ⚠ eigs with ONE output returns the EIGENVALUE
        case "trivial"
            [ph, Ax] = i_trivial(V, F, C, o);
        otherwise
            error('operators:gauge:method', 'Method must be diffusion, smoothest or trivial.');
    end

    e1 = cos(ph).*C.e1 + sin(ph).*C.e2;
    e1 = e1 ./ max(vecnorm(e1,2,2), eps);
    g = struct('e1', e1, 'e2', cross(C.normal, e1, 2), 'normal', C.normal, ...
               'phase', ph, 'method', char(method));

    % diagnostics, and the singularities -- via PRINCIPAL VALUES, the way rheome.detect.criticalPoints
    % localises winding, because that is the sum Poincare-Hopf applies to
    [ii,jj] = find(triu(C.A ~= 0, 1));
    aR = @(r,c) angle(full(C.Rt(sub2ind([nV nV], r, c))));
    % ⚠ TWO DIFFERENT RESIDUALS. .residual asks whether the frame is parallel with respect to ITS
    % OWN connection -- zero by construction for 'trivial'. .residualLC asks the same question
    % against the unmodified Levi-Civita transport, the number comparable across methods.
    g.residualLC = median(abs(wr(ph(jj) - ph(ii) - aR(ii,jj))));
    if method == "trivial"
        g.residual = median(abs(wr(ph(jj) - ph(ii) - full(Ax(sub2ind([nV nV], ii, jj))))));
    else
        g.residual = g.residualLC;
    end
    g.neighbourAngle = median(acosd(max(-1, min(1, sum(e1(ii,:).*e1(jj,:), 2)))));
    a=F(:,1); b=F(:,2); c=F(:,3);
    q = round(( wr(ph(b)-ph(a)-aR(a,b)) + wr(ph(c)-ph(b)-aR(b,c)) ...
              + wr(ph(a)-ph(c)-aR(c,a)) ) / (2*pi));
    g.singular = find(q ~= 0);  g.charge = q(g.singular);
end

% ----- the trivial connection: rheome.operators.trivial_connection with the default poles -----
function [ph, Ax] = i_trivial(V, F, C, o)
    sing = o.Singular(:);  chg = o.SingularCharge(:);
    if isempty(sing)                                     % default: two well-separated faces
        fcen = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:))/3;
        [~, i1] = max(vecnorm(fcen - mean(fcen,1), 2, 2));
        [~, i2] = max(vecnorm(fcen - fcen(i1,:), 2, 2));
        sing = [i1; i2];
    end
    if isempty(chg), chg = ones(numel(sing),1); end
    if numel(chg) ~= numel(sing)
        error('operators:gauge:charge', 'SingularCharge must match Singular in length.');
    end
    nE = size(unique(sort([F(:,[2 3]); F(:,[3 1]); F(:,[1 2])], 2), 'rows'), 1);
    chi = size(V,1) - nE + size(F,1);
    if sum(chg) ~= chi
        error('operators:gauge:budget', ...
            ['SingularCharge must sum to chi = %d (Poincare-Hopf), not %d. The face holonomies ' ...
             'sum to 2*pi*chi, so another total is inconsistent -- the solve would not fail, it ' ...
             'would place the remainder at the pinned face.'], chi, sum(chg));
    end
    T = rheome.operators.trivial_connection(V, F, C, sing, chg);
    ph = T.phase;  Ax = T.Ax;
end

% Author: Diellor Basha, 2026
