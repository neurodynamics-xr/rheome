function M = mass(V, F, type, opt)
% OPERATORS.MASS  Mass matrix of a triangle surface mesh.
%
%   M = rheome.operators.mass(V, F)                   % 'mixedVoronoi' (default)
%   M = rheome.operators.mass(V, F, 'galerkin')       % consistent (full) P1 FEM mass
%   M = rheome.operators.mass(V, F, 'barycentric')    % A/3 per incident triangle (diagonal)
%   M = rheome.operators.mass(V, F, 'mixedVoronoi', mollify=1e-6)
%
% The mass matrix M is the discrete inner product on scalar fields living on the
% mesh: <f,g> = f' M g. It is what makes the Laplace-Beltrami eigenmodes an
% M-orthonormal basis (phi_i' M phi_j = delta_ij).
%
%   mixedVoronoi (DEFAULT, "one mesh and mass", PI 2026-11-10)  diagonal; the mixed Voronoi area
%                (Meyer, Desbrun, Schroeder & Barr 2003) of the MOLLIFIED input mesh, recorded before any
%                flip: rheome.operators.intrinsic_delaunay's S.mass. The Cortical Flow Editor's and nxr's
%                "lumped mass" (the Hodge-Dirac star0).
%   galerkin     consistent P1 mass, (A/12) [2 1 1; 1 2 1; 1 1 2] per triangle (Brainstorm's pairing). Embedded.
%   barycentric  (A/3) I3 per triangle, the row sums of galerkin. Diagonal. Embedded. (Formerly named 'lumped';
%                nxr and the Editor call the mixed Voronoi mass "lumped", so that name is gone here.)
% Each assembles to a total mass equal to the surface area.
%
% INPUTS:
%   V    [nV x 3] vertex coordinates
%   F    [nF x 3] triangle vertex indices (1-based)
%   type 'mixedVoronoi' (default) | 'galerkin' | 'barycentric'
%   mollify  (name-value, mixedVoronoi only) relative mollification delta, default 1e-6 (intrinsic_delaunay's)
%
% OUTPUT:
%   M    [nV x nV] sparse symmetric positive-definite mass matrix
%
% See also: rheome.operators.laplace_beltrami, rheome.operators.intrinsic_delaunay
%
% Author: Diellor Basha, 2026

    arguments
        V (:,3) double, F (:,3) double
        type (1,1) string = "mixedVoronoi"
        opt.mollify (1,1) double = 1e-6
    end
    nV = size(V, 1);
    if lower(type) == "mixedvoronoi"
        S = rheome.operators.intrinsic_delaunay(V, F, mollify=opt.mollify);
        M = spdiags(S.mass, 0, nV, nV);
        return
    end

    % Per-face area from the cross product of two edge vectors.
    e1 = V(F(:,2), :) - V(F(:,1), :);
    e2 = V(F(:,3), :) - V(F(:,1), :);
    A  = 0.5 * sqrt(sum(cross(e1, e2, 2).^2, 2));   % [nF x 1]

    switch lower(type)
        case 'barycentric'
            % Each vertex collects A/3 from every incident triangle -> diagonal mass.
            w = repmat(A/3, 3, 1);
            idx = F(:);
            M = sparse(idx, idx, w, nV, nV);

        case 'galerkin'
            % Consistent P1 mass: diagonal entries A/6, off-diagonal (each edge) A/12.
            i = [F(:,1); F(:,2); F(:,3); ...                       % diagonal
                 F(:,1); F(:,2); F(:,1); F(:,3); F(:,2); F(:,3)];  % off-diagonal (both directions)
            j = [F(:,1); F(:,2); F(:,3); ...
                 F(:,2); F(:,1); F(:,3); F(:,1); F(:,3); F(:,2)];
            v = [A/6;  A/6;  A/6; ...
                 A/12; A/12; A/12; A/12; A/12; A/12];
            M = sparse(i, j, v, nV, nV);

        otherwise
            error('operators:mass:type', 'Unknown mass type ''%s'' (use ''mixedVoronoi'', ''galerkin'' or ''barycentric'').', type);
    end

    M = (M + M') / 2;   % enforce exact symmetry
end

% Author: Diellor Basha, 2026
