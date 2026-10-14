function [P, info] = geodesicpath(S, from, to, solver)
% GEOM.GEODESICPATH  A geodesic on the surface between two vertices, traced across the faces.
%
%   [P, info] = rheome.geom.geodesicpath(S, from, to)
%   [P, info] = rheome.geom.geodesicpath(S, from, to, solver)   % reuse a rheome.geom.geodesic solver
%
% A path ON the folded surface: steepest descent of the heat-method distance to `to` (rheome.geom.geodesic), from
% `from`, straight inside each triangle and crossing edges where the descent leaves it -- the Editor's traceGeodesic,
% on a distance that matches geometry-central. Its points lie on the mesh, so its length (the sum of its segments)
% is a length along the surface. Measured on fsaverage6 pial against the exact (edge-flip) geodesic: tGeodesicPath
% and the evidence of request 07dbab95. A path along mesh edges (rheome.geom.edgegraph) is 1-10 % longer.
%
% INPUTS:
%   S       surface (.Vertices, .Faces)        from, to   vertex indices
%   solver  (optional) from rheome.geom.geodesic on the same S
%
% OUTPUTS:
%   P     [n x 3] the path's points, from `from` to `to`: the two vertices and every edge crossing between them
%   info  .length (mesh units, along the surface)  .distance (the heat distance from `from` to `to`)
%         .faces [n-1 x 1] the face each segment lies in  .fallbacks (steps taken along an edge, where the descent
%         met a vertex or a flat; few on a fine mesh)  .solver
%
% See also: rheome.geom.geodesic, rheome.geom.edgegraph
%
% Author: Diellor Basha, 2026

    if nargin < 4, solver = []; end
    V = double(S.Vertices);  F = double(S.Faces);  nF = size(F, 1);
    [d, solver] = rheome.geom.geodesic(S, to, solver);
    if ~isfinite(d(from)), error('rheome:geodesicpath:unreachable', 'Vertex %d cannot reach vertex %d on this surface.', from, to); end
    % the faces across each edge and around each vertex (kept in the solver for the next path)
    if ~isfield(solver, 'topology')
        pairs = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
        [~, ~, eid] = unique(sort(pairs, 2), 'rows');
        [ef, o] = sort(eid);  fo = repmat((1:nF)', 3, 1);  fo = fo(o);
        first = [true; diff(ef) > 0];  E2F = [fo, fo];  last = [first(2:end); true];
        E2F = E2F(first, :);  E2F(:, 2) = fo(last);              % a boundary edge: the same face twice
        [vf, o] = sort(F(:));  fv = repmat((1:nF)', 3, 1);
        solver.topology = struct('FE', reshape(eid, nF, 3), 'E2F', E2F, 'V2F', {mat2cell(fv(o), accumarray(vf, 1, [size(V,1) 1]))});
    end
    FE = solver.topology.FE;  E2F = solver.topology.E2F;  V2F = solver.topology.V2F;
    e1 = V(F(:,2),:) - V(F(:,1),:);  e2 = V(F(:,3),:) - V(F(:,1),:);

    P = V(from, :);  faces = [];  fb = 0;
    f = 0;  b = [];  v = from;                                    % at a vertex: v; on an edge of face f: barycentric b
    for step = 1:20 * size(V, 1)
        if v == to, break; end
        if v > 0                                                  % at a vertex: the incident face the descent enters
            cand = V2F{v};  ok = false;
            for f = cand(:)'
                c = find(F(f,:) == v);  b = zeros(1, 3);  b(c) = 1;
                [lam, bn] = exitOf(f, b);
                if lam > 0 && bn(c) < 1 - 1e-9, ok = true; break; end
            end
            if ~ok                                                % no face descends from here: an edge step down
                nb = unique(F(cand, :));  [~, j] = min(d(nb));
                if d(nb(j)) >= d(v), error('rheome:geodesicpath:stuck', 'The descent stopped at vertex %d.', v); end
                fb = fb + 1;  P(end+1, :) = V(nb(j), :);  faces(end+1, 1) = cand(1);  v = nb(j);  continue %#ok<AGROW>
            end
        else
            [lam, bn] = exitOf(f, b);
        end
        if any(F(f,:) == to)                                      % the target's face: straight to it
            P(end+1, :) = V(to, :);  faces(end+1, 1) = f;  v = to;  break %#ok<AGROW>
        end
        if ~(lam > 0)                                             % flat or reversed: step to the face's lowest vertex
            [~, j] = min(d(F(f,:)));  fb = fb + 1;
            P(end+1, :) = V(F(f,j), :);  faces(end+1, 1) = f;  v = F(f,j);  continue %#ok<AGROW>
        end
        P(end+1, :) = bn * V(F(f,:), :);  faces(end+1, 1) = f; %#ok<AGROW>
        z = find(bn <= 1e-9);
        if numel(z) >= 2                                          % left through a vertex
            [~, j] = max(bn);  v = F(f, j);  continue
        end
        g = E2F(FE(f, z), :);  g = g(g ~= f);
        if isempty(g), error('rheome:geodesicpath:boundary', 'The descent left the surface across a boundary edge.'); end
        v = 0;  nb = zeros(1, 3);                                 % the same point, in the next face's barycentrics
        for c = 1:3, k = find(F(g(1),:) == F(f,c)); if ~isempty(k), nb(k) = bn(c); end, end
        f = g(1);  b = nb;
    end
    if v ~= to, error('rheome:geodesicpath:stuck', 'The descent did not reach vertex %d.', to); end
    info = struct('length', sum(vecnorm(diff(P), 2, 2)), 'distance', d(from), 'faces', faces, 'fallbacks', fb, 'solver', solver);

    function [lam, bn] = exitOf(f, b)
    % where the descent leaves face f from barycentric point b: lam > 0 the step, bn the exit's barycentrics
        A = [e1(f,:); e2(f,:)]';                                  % 3 x 2
        g3 = A * ((A' * A) \ [d(F(f,2)) - d(F(f,1)); d(F(f,3)) - d(F(f,1))]);   % the in-plane gradient of d
        st = (A' * A) \ (A' * (-g3));                             % the descent in (s, t) = (b2, b3)
        db = [-st(1) - st(2), st(1), st(2)];
        neg = db < -1e-14;
        if ~any(neg), lam = 0;  bn = b;  return, end
        lam = min(b(neg) ./ -db(neg));
        bn = max(b + lam * db, 0);  bn = bn / sum(bn);
        if any(b(~neg) < -1e-12), lam = 0; end
    end
end

% Author: Diellor Basha, 2026
