function ladderWriteMeshes(out)
% The synthetic meshes of tSurfaceLadder's golden (ladderCaptureGolden.py labels them with nsp's ladder.py):
% a symmetric ico-3 sphere (exact ties everywhere), a jittered ico-4 sphere, and two disjoint ico-3
% spheres of unequal size (disconnected tiles: infinite and NaN sort keys).
rng(20261014);
[V3, F3] = rheome.geom.icosphere(3);  [V4, F4] = rheome.geom.icosphere(4);
V4 = V4 .* (1 + 0.05 * randn(size(V4, 1), 1)) * 0.07;
V2 = [V3 * 0.05; V3 * 0.03 + [0.2 0 0]];  F2 = [F3; F3 + size(V3, 1)];
meshes = {struct('V', V3 * 0.07, 'F', F3), struct('V', V4, 'F', F4), struct('V', V2, 'F', F2)};
save(out, 'meshes', '-v7');
end
