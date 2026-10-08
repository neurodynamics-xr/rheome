function Y = localize(basis, vertices, gains)
% FILTERS.LOCALIZE  Localize a filter (bank) at cortical vertices -> atoms.
%
%   Y = rheome.filters.localize(basis, vertex, gains)
%   Y = rheome.filters.localize(basis, vertices, gains)
%
% Drops a unit delta at each requested vertex and filters it through the spectral
% gain(s). The result is the filter "localized" at that vertex -- a cortical atom.
% This is exactly the delta-seed step of Brainstorm's bst_eigenfilter('Atom') for the
% Laplace–Beltrami case: c0 = Phi' M delta,  atom = Phi * (g(lambda) .* c0).
%
% INPUTS:
%   basis     struct from rheome.eigen.modes
%   vertices  scalar vertex index, or [P x 1] indices to seed
%   gains     [K x 1] single filter, or [K x Nf] filterbank, or a function handle
%
% OUTPUT:
%   Y  single vertex : [nV x 1] (single filter) or [nV x Nf] (filterbank)
%      P vertices     : [nV x P] (single filter) or [nV x P x Nf] (filterbank)
%
% See also: rheome.filters.apply, rheome.filters.heat, rheome.filters.mexhat, rheome.show.surface
%
% Author: Diellor Basha, 2026

    vertices = vertices(:);
    if any(vertices < 1 | vertices > basis.nV | mod(vertices,1) ~= 0)
        error('filters:localize:vertex', 'vertices must be integer indices in 1..%d.', basis.nV);
    end
    P = numel(vertices);
    X = sparse(vertices, (1:P)', 1, basis.nV, P);   % delta at each seed
    Y = rheome.filters.apply(basis, full(X), gains);
end

% Author: Diellor Basha, 2026
