function [distance, solver] = geodesic(S, sourceVertices, solver)
% GEOM.GEODESIC  Single-/multi-source geodesic distance by Crane's heat method.
%
%   [distance, solver] = rheome.geom.geodesic(S, sourceVertices)          % builds the solver
%   [distance, solver] = rheome.geom.geodesic(S, sourceVertices, solver)  % reuse cached factors
%
% Distance ON THE SURFACE (not the ambient embedding) from a set of source vertices, via the
% heat method (Crane, Weischedel & Wardetzky 2013): diffuse a heat spike briefly, march down
% the normalised heat gradient, and integrate. Built entirely from the module's cotan stiffness
% + mass (rheome.operators.laplace_beltrami) and per-face gradient (rheome.operators.face_gradient). Ported
% from nxr-compute/src/geodesic.cpp (geometry-central HeatMethodDistanceSolver).
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
%   distance [nV x 1] geodesic distance (source ~ 0, increasing away), metres (mesh units)
%   solver   cached factorisations for reuse
%
% See also: rheome.geom.karcher_mean, rheome.operators.laplace_beltrami, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(solver)
        vertices = double(S.Vertices);  faces = double(S.Faces);
        [stiffness, mass] = rheome.operators.laplace_beltrami(vertices, faces);
        faceGradient = rheome.operators.face_gradient(vertices, faces);
        edges = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
        edgeVectors = vertices(edges(:,1),:) - vertices(edges(:,2),:);
        meanEdgeLength = mean(sqrt(sum(edgeVectors.^2, 2)));
        timeStep = meanEdgeLength^2;                       % Crane's default diffusion time
        solver.heatOp       = decomposition(mass + timeStep*stiffness, 'chol');
        solver.poissonOp    = decomposition(stiffness + 1e-10*mass, 'chol'); % pin constant nullspace
        solver.faceGradient = faceGradient;
        solver.nV = faceGradient.nV;
    end
    fg = solver.faceGradient;

    sourceIndicator = zeros(solver.nV, 1);
    sourceIndicator(sourceVertices) = 1;
    heat = solver.heatOp \ sourceIndicator;                % 1. brief heat diffusion

    gradX = fg.Gx*heat;  gradY = fg.Gy*heat;  gradZ = fg.Gz*heat;   % 2. per-face heat gradient
    gradNorm = sqrt(gradX.^2 + gradY.^2 + gradZ.^2);
    gradNorm(gradNorm < eps) = 1;
    directionX = -gradX./gradNorm;  directionY = -gradY./gradNorm;  directionZ = -gradZ./gradNorm;

    faceArea = fg.FaceArea;                                % 3. integrated divergence (matches L)
    divergenceToVertex = fg.Gx'*(faceArea.*directionX) ...
                       + fg.Gy'*(faceArea.*directionY) ...
                       + fg.Gz'*(faceArea.*directionZ);
    distance = solver.poissonOp \ divergenceToVertex;      % 4. Poisson recovers distance
    distance = distance - min(distance);                   % shift source set to ~0
end

% Author: Diellor Basha, 2026
