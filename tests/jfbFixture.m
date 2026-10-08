function fx = jfbFixture()
% JFBFIXTURE  Tiny joint fixture: ico3, K = 40 modes, 64 samples at 64 Hz.
%
% Uses this repo's packages ON PURPOSE -- the fixture is test scaffolding, not the
% deliverable. @jointfilterbank itself must never reference them (see tGfbGenericity).
    persistent cache
    if ~isempty(cache), fx = cache; return; end

    R = 0.100;
    [V, F] = rheome.geom.icosphere(3);                 % 642 vertices
    V = R * (V ./ vecnorm(V, 2, 2));
    [L, M] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
    b = rheome.eigen.modes(L, M, 40);

    fx = struct();
    fx.Phi = b.Phi;  fx.Lambda = b.Lambda;  fx.Mass = b.Mass;  fx.nV = b.nV;
    fx.V = V;  fx.F = F;  fx.L = L;  fx.M = M;
    fx.N = 64;  fx.fs = 64;
    fx.T   = rheome.graphtransform.eigen(b.Phi, b.Mass, b.Lambda);
    fx.gfb = rheome.graphfilterbank(b.Lambda, 'Wavelet','itersine', 'NumFilters',3);
    cache = fx;
end

% Author: Diellor Basha, 2026
