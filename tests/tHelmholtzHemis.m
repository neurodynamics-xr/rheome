classdef tHelmholtzHemis < matlab.unittest.TestCase
% The Hodge solve on a two-hemisphere surface: each hemisphere is pinned on its own whether the surface
% gives its hemispheres as index sets, as a per-vertex index, or not at all (connected components), and the
% potentials equal those of each hemisphere solved alone.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function twoHemispheresSolveAsTheirParts(tc)
            [S1, J1] = i_sphere([0 0 0]);  [S2, J2] = i_sphere([3 0 0]);
            n = size(S1.Vertices, 1);
            S = struct('Vertices', [S1.Vertices; S2.Vertices], 'Faces', [S1.Faces; S2.Faces + n], ...
                       'VertNormals', [S1.VertNormals; S2.VertNormals]);
            J = [J1; J2];
            H1 = rheome.differential.helmholtz(J1, S1);  H2 = rheome.differential.helmholtz(J2, S2);
            want = [H1.Phi H1.Psi; H2.Phi H2.Psi];
            idx = [ones(n, 1); 2*ones(n, 1)];
            for hemi = {[], {(1:n)', (n+1:2*n)'}, idx}
                Sh = S;  if ~isempty(hemi{1}), Sh.Hemi = hemi{1}; end
                lastwarn('');
                H = rheome.differential.helmholtz(J, Sh);
                [~, id] = lastwarn();
                tc.verifyEmpty(id, 'the solve warned');
                tc.verifyEqual([H.Phi H.Psi], want, 'AbsTol', 1e-9 * max(abs(want(:))));
                phi = rheome.differential.poisson(H.Div, Sh);
                tc.verifyEqual(mean(phi(1:n)), 0, 'AbsTol', 1e-12);       % mean-zero per hemisphere
            end
        end

        function refusesAWrongLengthIndex(tc)
            [S, J] = i_sphere([0 0 0]);  S.Hemi = [1; 2];
            tc.verifyError(@() rheome.differential.helmholtz(J, S), 'differential:hemis:size');
        end
    end
end

% an icosphere at c with a tangent field (source + vortex about a vertex)
function [S, J] = i_sphere(c)
    [V, F] = rheome.geom.icosphere(3);
    n = V ./ vecnorm(V, 2, 2);
    S = struct('Vertices', V + c, 'Faces', F, 'VertNormals', n);
    g = exp(-sum((V - V(2,:)).^2, 2));  gr = (V - V(2,:)) .* g;  gr = gr - sum(gr .* n, 2) .* n;
    J = reshape((gr + cross(n, gr, 2))', [], 1);
end

% Author: Diellor Basha, 2026
