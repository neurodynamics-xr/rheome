function [V, F] = icosphere(nSub)
% GEOM.ICOSPHERE  Geodesic sphere by recursive icosahedron subdivision (FreeSurfer-style).
%
%   [V, F] = rheome.geom.icosphere(nSub)
%
% Builds a near-uniform unit-sphere triangulation by subdividing an icosahedron nSub
% times (each triangle -> 4, midpoints projected to the sphere). This is exactly the
% construction FreeSurfer uses for its ico spheres:
%   nSub : 0    1    2    3    4     5(=ico5)  6(=ico7)
%   nV   : 12   42   162  642  2562  10242     40962
%   nF   : 20   80   320  1280 5120  20480     81920
%
% A pristine ico sphere (uniform triangles) makes the LBO eigenfunctions match the
% spherical harmonics to high accuracy -- the clean reference for the analytic-vs-
% numerical PDE comparison (unlike the distorted cortical registration sphere).
%
% INPUT:
%   nSub : number of subdivisions (5 for ico5). Default 5.
% OUTPUTS:
%   V    [nV x 3] unit-sphere vertices
%   F    [nF x 3] triangle indices (1-based)
%
% Example:
%   [V, F] = rheome.geom.icosphere(3);
%   assert(size(V, 1) == 642 && size(F, 1) == 1280 && max(abs(vecnorm(V, 2, 2) - 1)) < 1e-12)
%
% See also: demos.sphere_pde, rheome.operators.laplace_beltrami
%
% Author: Diellor Basha, 2026

    if nargin < 1 || isempty(nSub), nSub = 5; end

    % --- base icosahedron (12 vertices, 20 faces) ---
    t = (1 + sqrt(5)) / 2;
    V = [-1  t  0;  1  t  0; -1 -t  0;  1 -t  0;
          0 -1  t;  0  1  t;  0 -1 -t;  0  1 -t;
          t  0 -1;  t  0  1; -t  0 -1; -t  0  1];
    V = V ./ vecnorm(V, 2, 2);
    F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12;
         2 6 10; 6 12 5; 12 11 3; 11 8 7; 8 2 9;
         4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10;
         5 10 6; 3 5 12; 7 3 11; 9 7 8; 10 9 2];

    % --- subdivide ---
    for s = 1:nSub
        mid = containers.Map('KeyType','char','ValueType','double');
        nF  = size(F,1);
        Fn  = zeros(4*nF, 3);
        for f = 1:nF
            a = F(f,1); b = F(f,2); c = F(f,3);
            ab = i_mid(a,b);  bc = i_mid(b,c);  ca = i_mid(c,a);
            Fn(4*f-3,:) = [a  ab ca];
            Fn(4*f-2,:) = [b  bc ab];
            Fn(4*f-1,:) = [c  ca bc];
            Fn(4*f  ,:) = [ab bc ca];
        end
        F = Fn;
    end

    % nested helper: get-or-create the (normalized) midpoint of edge (i,j)
    function m = i_mid(i, j)
        key = sprintf('%d_%d', min(i,j), max(i,j));
        if isKey(mid, key)
            m = mid(key);
        else
            p = (V(i,:) + V(j,:)) / 2;
            p = p / norm(p);
            V(end+1,:) = p;         %#ok<AGROW>
            m = size(V,1);
            mid(key) = m;
        end
    end
end

% Author: Diellor Basha, 2026
