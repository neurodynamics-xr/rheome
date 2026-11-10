function s = solve(pkg, keep, positions)
% REGISTER.SOLVE  One Hodge-Dirac solve as the registration reads it: the basis with the glue it is on.
%
%   s = rheome.register.solve(pkg, keep, positions)        % pkg = rheome.eigen.hodge_package(...); [~, ~, keep] = rheome.geom.glue(...)
%
% pkg        the modes: P [V x K], mu, kind, mass (rheome.eigen.hodge_package)
% keep       [V x 1] each glued vertex's input vertex (1-based), 0 on the bridge, kept vertices first (rheome.geom.glue)
% positions  [V x 3] the glued positions
% s          struct with those fields; s.kept the number of kept vertices. As nxr.register.modes.Solve.
%
% See also: rheome.register.mode_registration
%
% Author: Diellor Basha, 2026
arguments, pkg (1,1) struct, keep (:,1) double, positions (:,3) double, end
if size(pkg.P, 1) ~= numel(keep) || size(positions, 1) ~= numel(keep)
    error('rheome:register:solve', 'the basis has %d vertices, keep %d, positions %d', size(pkg.P, 1), numel(keep), size(positions, 1));
end
s = struct('P', pkg.P, 'mu', pkg.mu, 'kind', pkg.kind, 'mass', pkg.mass, 'keep', keep, 'positions', positions, 'kept', nnz(keep > 0));
end
