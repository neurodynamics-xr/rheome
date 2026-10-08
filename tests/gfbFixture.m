function fx = gfbFixture(nSub, K)
% GFBFIXTURE  Small real LBO basis for graphfilterbank tests.
%   fx = gfbFixture()        % ico4 (2562 verts), K = 100 -- 4.7 verts/half-wavelength
%   fx = gfbFixture(nSub, K)
%
% Uses this repo's packages ON PURPOSE: the fixture is test scaffolding, not the
% deliverable. @graphfilterbank itself must never reference them (see tGfbGenericity).
    if nargin < 1 || isempty(nSub), nSub = 4;   end
    if nargin < 2 || isempty(K),    K    = 100; end
    persistent cache
    key = sprintf('n%dk%d', nSub, K);
    if isstruct(cache) && isfield(cache, key), fx = cache.(key); return; end

    R = 0.100;                                   % m, FreeSurfer registration-sphere radius
    [V, F] = rheome.geom.icosphere(nSub);
    V = R * (V ./ vecnorm(V, 2, 2));
    [L, M] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
    b = rheome.eigen.modes(L, M, K);

    fx = struct('Phi', b.Phi, 'Lambda', b.Lambda, 'Mass', b.Mass, ...
                'nV', b.nV, 'V', V, 'F', F, 'L', L, 'M', M, 'R', R);
    cache.(key) = fx;
end

% Author: Diellor Basha, 2026
