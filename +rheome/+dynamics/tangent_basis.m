function tb = tangent_basis(surface)
% DYNAMICS.TANGENT_BASIS  Per-vertex orthonormal tangent frame from the vertex normals.
%
%   tb = rheome.dynamics.tangent_basis(surface)   % surface needs .VertNormals
%
% Returns tb.e1, tb.e2 [nV x 3] -- an orthonormal basis of each vertex's tangent plane -- plus
% tb.normal (unit vertex normal). e2 = normal x e1. Same construction as Brainstorm's optical-flow
% basis_vertices (the "[z-y x-z y-x]" trick), with a fallback for near-[1 1 1] normals. Used to
% express a tangent velocity field v = a*e1 + b*e2 by its two scalar coordinates (a,b).
%
% See also: rheome.dynamics.opticalflow_scalar, rheome.dynamics.opticalflow_vector
%
% Author: Diellor Basha, 2026

    unitNormal = surface.VertNormals ./ max(vecnorm(surface.VertNormals, 2, 2), eps);
    frameE1 = [unitNormal(:,3) - unitNormal(:,2), ...
               unitNormal(:,1) - unitNormal(:,3), ...
               unitNormal(:,2) - unitNormal(:,1)];
    degenerate = abs(unitNormal * (ones(3,1)/sqrt(3))) > 0.97;      % normal ~ [1 1 1]
    frameE1(degenerate,:) = [unitNormal(degenerate,2), -unitNormal(degenerate,1), zeros(sum(degenerate),1)];
    frameE1 = frameE1 ./ max(vecnorm(frameE1, 2, 2), eps);
    frameE2 = cross(unitNormal, frameE1, 2);
    frameE2 = frameE2 ./ max(vecnorm(frameE2, 2, 2), eps);
    tb = struct('e1', frameE1, 'e2', frameE2, 'normal', unitNormal);
end

% Author: Diellor Basha, 2026
