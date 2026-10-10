function T = trivial_connection(V, F, C, faces, charges)
% OPERATORS.TRIVIAL_CONNECTION  Crane's trivial connection on the dual complex: a per-vertex frame
% that is exactly parallel except at prescribed faces, where its index is the prescribed integer.
%
%   T = rheome.operators.trivial_connection(V, F, C, faces, charges)
%       V [nV x 3], F [nF x 3], C = rheome.operators.connection_laplacian(V,F),
%       faces [n x 1] face indices, charges [n x 1] integers with sum(charges) = chi.
%
% THE EQUATION (Crane, Desbrun & Schroeder 2010, "Trivial connections on discrete surfaces"):
% find the least-norm edge 1-form x with, on every face f,
%
%       (d1 * x)_f  =  2*pi*k_f  -  K_f,         K_f = wrap( (d1 * argE)_f )  in (-pi, pi]
%
% where argE is the Levi-Civita transport of C on each edge and K_f is the face's holonomy, i.e.
% the discrete Gaussian curvature carried by that face. sum_f K_f = 2*pi*chi (Gauss-Bonnet), so
% the system is consistent exactly when sum(k) = chi -- the Poincare-Hopf budget.
% nxr-compute (src/direction_field.cpp, computeTrivialConnection) solves the same equation on the
% PRIMAL complex -- a face field with vertex singularities, rhs -K_v + 2*pi*sigma_v. This is its
% dual: the toolbox's frames live on vertices, so the singularities live on faces. Measured on the
% ico5 sphere with both at the poles: the two agree to 0.31 deg median.
%
% ⚠⚠ THE WRAP IS THE WHOLE POINT. d1*argE is NOT the curvature: C's charts give each edge a
% principal-value transport, so (d1*argE)_f = K_f + 2*pi*n_f with integer chart parts n_f
% (nonzero on 4020 cortex faces, sum -chi) and sum_f (d1*argE)_f = 0 identically. An earlier version
% solved d1*x = 2*pi*k - d1*argE, which prescribes index k_f - n_f: with k = [+1 -1] the frame
% wound on 2560 faces (55.4 deg between neighbours) and the "budget is zero, not chi" note was this
% artefact. With the wrap, k = [+1 +1] at two faces gives exactly two +1 singularities there.
%
% OUTPUT (struct T)
%   .phase [nV x 1]  frame angle in C's charts (integrated over a spanning tree; path-independent)
%   .Ax    [nV x nV] sparse, the adjusted transport angle on each edge (antisymmetric)
%   .x     [nE x 1]  the least-norm adjustment     .E [nE x 2]  .d1 [nF x nE]  .argE [nE x 1]
%   .K     [nF x 1]  face holonomy (curvature), radians     .n [nF x 1] chart integer parts
%   .k     [nF x 1]  the prescribed index per face           .chi
%
% See also: rheome.operators.gauge, rheome.flow.directionfield, rheome.operators.connection_laplacian
%
% Author: Diellor Basha, 2026

    F = double(F);  nV = size(V,1);  nF = size(F,1);
    faces = faces(:);  charges = charges(:);
    [E,~,~] = unique(sort([F(:,[2 3]); F(:,[3 1]); F(:,[1 2])], 2), 'rows');
    nE = size(E,1);  chi = nV - nE + nF;
    if numel(faces) ~= numel(charges)
        error('operators:trivial_connection:charge', 'charges must match faces in length.');
    end
    if sum(charges) ~= chi
        error('operators:trivial_connection:budget', ...
            ['The charges must sum to chi = %d (Poincare-Hopf), not %d: sum_f K_f = 2*pi*chi, so ' ...
             'any other total makes d1*x = 2*pi*k - K inconsistent.'], chi, sum(charges));
    end

    argE = angle(full(C.Rt(sub2ind([nV nV], E(:,1), E(:,2)))));
    dh = [F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)];
    [~, de] = ismember(sort(dh,2), E, 'rows');
    sgE = 2*double(dh(:,1) == E(de,1)) - 1;
    d1 = sparse(repmat((1:nF)',3,1), de, sgE, nF, nE);   % oriented face-edge incidence

    hol = d1*argE;
    n = round(hol/(2*pi));  K = hol - 2*pi*n;            % principal holonomy = curvature
    k = accumarray(faces, charges, [nF 1]);
    rhs = 2*pi*k - K;
    L = d1*d1';  L(1,1) = L(1,1) + 1;                    % pin the constant nullspace (sum rhs = 0)
    x = d1' * (L \ rhs);
    argE2 = argE + x;

    % integrate over a spanning tree: path-independent now, so the tree cannot matter
    Ax = sparse([E(:,1);E(:,2)], [E(:,2);E(:,1)], [argE2; -argE2], nV, nV);
    et = dfsearch(minspantree(graph(E(:,1), E(:,2)), 'Root', 1), 1, 'edgetonew');
    ph = zeros(nV,1);
    for r = 1:size(et,1), ph(et(r,2)) = ph(et(r,1)) + Ax(et(r,1), et(r,2)); end

    T = struct('phase', ph, 'Ax', Ax, 'x', x, 'E', E, 'd1', d1, 'argE', argE, ...
               'K', K, 'n', n, 'k', k, 'chi', chi);
end

% Author: Diellor Basha, 2026
