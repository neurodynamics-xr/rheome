function U = impulse(basis, vertex, G, domain, nKeep, dir)
% FILTERS.IMPULSE  Space-time impulse response of a joint filter at a vertex delta.
%
%   U = rheome.filters.impulse(basis, vertex, G, domain [,nKeep])          % scalar (LBO) basis
%   U = rheome.filters.impulse(dbasis, vertex, G, domain, nKeep, dir)      % Dirac (vector) basis
%
% Seeds a delta at a vertex (with a flat temporal spectrum) and propagates it through a
% joint filter, giving the filter's space-time impulse response -- the "dynamics" the
% filter describes. Domain-aware, exactly like bst_eigenfilter('Atom'):
%   'ts'  G = g(lambda,t)      -> coeffs = G .* c0
%   'js'  G = g(lambda,omega)  -> coeffs = real(ifft(G)) .* c0
% where c0 = Phi' * M * delta are the seed's mode coefficients.
%
% SCALAR basis (rheome.eigen.modes): the seed is a scalar delta at the vertex; returns the
% scalar field U [nV x nT].
% DIRAC basis (rheome.eigen.dirac_frame): the seed is an oriented DIPOLE 'dir' [3x1] at the vertex
% (embedded in the imaginary quaternion slots); returns the decoded VECTOR field
% [3nV x nT] (rows [x1,y1,z1, x2,...]), via rheome.forward.reconstruct.
%
% INPUTS:
%   basis   eigenbasis (scalar from rheome.eigen.modes, or Dirac from rheome.eigen.dirac_frame)
%   vertex  seed vertex index
%   G       [K x N] kernel evaluated on its native axis (time for 'ts', freq for 'js')
%   domain  'ts' | 'js'
%   nKeep   ('js' only) keep the first nKeep time samples (default N)
%   dir     (Dirac only) [3x1] dipole direction (default [1;0;0])
%
% OUTPUT:
%   U       [nV x nT] scalar field, or [3nV x nT] vector field (Dirac)
%
% See also: rheome.eigen.dirac_frame, rheome.forward.reconstruct, rheome.filters.wave, rheome.show.filmstrip, rheome.show.vectors
%
% Author: Diellor Basha, 2026

    isDirac = isfield(basis, 'nVert');
    if isDirac
        if nargin < 6 || isempty(dir), dir = [1;0;0]; end
        seed = zeros(4*basis.nVert, 1);
        seed((vertex-1)*4 + (2:4)) = dir(:);          % oriented dipole in imaginary slots
        c0 = basis.Phi' * (basis.Mass * seed);        % [nModes x 1]
    else
        c0 = basis.Phi' * (basis.Mass(:, vertex));    % scalar delta -> [K x 1]
    end

    switch lower(domain)
        case 'ts'
            coeffs = G .* c0;
        case 'js'
            Gt = real(ifft(G, [], 2));                % joint-spectral -> eigen-time response
            if nargin >= 5 && ~isempty(nKeep), Gt = Gt(:, 1:min(nKeep, size(Gt,2))); end
            coeffs = Gt .* c0;
        otherwise
            error('filters:impulse:domain', 'domain must be ''ts'' or ''js''.');
    end

    if isDirac
        U = rheome.forward.reconstruct(coeffs, basis);       % [3nV x nT] vector field
    else
        U = basis.Phi * coeffs;                       % [nV x nT] scalar field
    end
end

% Author: Diellor Basha, 2026
