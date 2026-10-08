function S = fbd_sphere(nSub, K)
% FBD_SPHERE  An ico sphere with its LBO eigenbasis, cached.
%
%   S = fbd_sphere(nSub, K)      nSub 3 -> 642 verts, 4 -> 2562
%
% R = 100 mm, FreeSurfer's registration-sphere radius. The sphere is used because every
% quantity the demos plant has a closed form on it: a degree-l harmonic sits at
% lambda = l(l+1)/R^2, and geodesic distance is R*acos(v.s).
%
% Author: Diellor Basha, 2026

    persistent cache
    key = sprintf('n%dk%d', nSub, K);
    if isstruct(cache) && isfield(cache, key), S = cache.(key); return; end

    R = 0.100;
    [V, F] = rheome.geom.icosphere(nSub);
    V = R * (V ./ vecnorm(V, 2, 2));
    [L, M] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
    b = rheome.eigen.modes(L, M, K);

    S = struct('V', V, 'F', F, 'nV', b.nV, 'R', R, ...
               'Phi', b.Phi, 'Lambda', b.Lambda, 'Mass', b.Mass, ...
               'T', rheome.graphtransform.eigen(b.Phi, b.Mass, b.Lambda));
    cache.(key) = S;
end

% Author: Diellor Basha, 2026
