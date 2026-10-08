function M = mass(V, F, type)
% OPERATORS.MASS  Finite-element mass matrix of a triangle surface mesh.
%
%   M = rheome.operators.mass(V, F)              % 'galerkin' (default, matches Brainstorm)
%   M = rheome.operators.mass(V, F, 'galerkin')  % consistent (full) FEM mass
%   M = rheome.operators.mass(V, F, 'lumped')    % barycentric lumped (diagonal) mass
%
% The mass matrix M is the discrete inner product on scalar fields living on the
% mesh: <f,g> = f' M g. It is what makes the Laplace-Beltrami eigenmodes an
% M-orthonormal basis (phi_i' M phi_j = delta_ij). Brainstorm's LBO pairs the
% cotan stiffness with the GALERKIN (consistent) mass, so that is the default here.
%
% For a linear (P1) element on a triangle of area A the local mass matrices are:
%   galerkin : (A/12) * [2 1 1; 1 2 1; 1 1 2]     (consistent; couples neighbours)
%   lumped   : (A/3)  * I3                         (row-sum of galerkin; diagonal)
% Both assemble to a total mass equal to the surface area (sum(M(:)) == area).
%
% INPUTS:
%   V    [nV x 3] vertex coordinates
%   F    [nF x 3] triangle vertex indices (1-based)
%   type 'galerkin' (default) | 'lumped'
%
% OUTPUT:
%   M    [nV x nV] sparse symmetric positive-definite mass matrix
%
% See also: rheome.operators.laplace_beltrami
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(type), type = 'galerkin'; end
    nV = size(V, 1);

    % Per-face area from the cross product of two edge vectors.
    e1 = V(F(:,2), :) - V(F(:,1), :);
    e2 = V(F(:,3), :) - V(F(:,1), :);
    A  = 0.5 * sqrt(sum(cross(e1, e2, 2).^2, 2));   % [nF x 1]

    switch lower(type)
        case 'lumped'
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
            error('operators:mass:type', 'Unknown mass type ''%s'' (use ''galerkin'' or ''lumped'').', type);
    end

    M = (M + M') / 2;   % enforce exact symmetry
end

% Author: Diellor Basha, 2026
