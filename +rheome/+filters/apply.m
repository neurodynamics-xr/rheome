function Y = apply(basis, X, gains)
% FILTERS.APPLY  Spectral filtering on the surface:  Y = Phi * (g(lambda) .* (Phi' M X)).
%
%   Y = rheome.filters.apply(basis, X, gains)
%
% Projects the field(s) X onto the eigenbasis, scales each mode coefficient by the
% spectral gain g(lambda), and reconstructs. This is g(L) applied to X -- analysis
% and synthesis in one step -- the operation underneath every atom.
%
% INPUTS:
%   basis  struct from rheome.eigen.modes (.Phi, .Lambda, .Mass, .nV)
%   X      [nV x nS] one or more vertex fields (e.g. a delta, or coordinate maps)
%   gains  spectral gains, any of:
%            [K x 1]           a single filter
%            [K x Nf]          a filterbank (Nf filters)  -> extra output dimension
%            function handle   g(lambda) evaluated on basis.Lambda
%
% OUTPUT:
%   Y      [nV x nS]        if a single filter
%          [nV x nS x Nf]   if a filterbank (nS>1)
%          [nV x Nf]        if a filterbank applied to a single field (nS==1)
%
% See also: rheome.filters.localize, rheome.eigen.modes
%
% Author: Diellor Basha, 2026

    if isa(gains, 'function_handle'), gains = gains(basis.Lambda); end
    if size(gains, 1) ~= numel(basis.Lambda)
        error('filters:apply:gains', 'gains has %d rows but the basis has %d modes.', ...
            size(gains, 1), numel(basis.Lambda));
    end

    C  = basis.Phi' * (basis.Mass * X);   % [K x nS] mode coefficients
    nS = size(X, 2);
    Nf = size(gains, 2);

    if Nf == 1
        Y = basis.Phi * (gains .* C);     % [nV x nS]
    else
        Y = zeros(basis.nV, nS, Nf);
        for f = 1:Nf
            Y(:, :, f) = basis.Phi * (gains(:, f) .* C);
        end
        if nS == 1
            Y = reshape(Y, basis.nV, Nf); % [nV x Nf] convenience for single-field banks
        end
    end
end

% Author: Diellor Basha, 2026
