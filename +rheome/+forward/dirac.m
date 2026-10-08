function Gm = dirac(Gain, dbasis)
% FORWARD.DIRAC  Project an unconstrained leadfield into the Dirac eigenbasis.
%
%   Gm = rheome.forward.dirac(Gain, dbasis)
%
% Spectral change-of-basis of the forward model: each per-vertex gain 3-vector is
% embedded as a pure-imaginary quaternion psi = [0, gx, gy, gz] and projected onto
% the B-orthonormal Dirac eigenvectors,
%       Gm = Psi' * (B * Phi)         [nCh x nModes],   B = kron(Mass, I4).
% Faithful port of Brainstorm's bst_dirac TRANSFORM (done whole-brain here; the
% cortex LBO is block-diagonal across hemispheres, so this equals the per-hemi loop).
%
% INPUTS:
%   Gain    [nCh x 3nV] unconstrained leadfield (columns [gx,gy,gz] per vertex)
%   dbasis  Dirac eigenbasis from rheome.eigen.dirac_frame (.Phi, .Mass, .nVert, .nModes)
%
% OUTPUT:
%   Gm      [nCh x nModes] mode-forward (the "Dirac eigenmode leadfield")
%
% See also: rheome.eigen.dirac_frame, rheome.forward.reconstruct, rheome.inverse.dirac
%
% Author: Diellor Basha, 2026

    G   = double(Gain);
    nV  = dbasis.nVert;
    nCh = size(G, 1);
    if size(G, 2) ~= 3*nV
        error('forward:dirac:size', 'Gain has %d columns but the basis expects 3*nV = %d.', size(G,2), 3*nV);
    end

    % Embed the gain as a pure-imaginary quaternion field Psi [4nV x nCh] (w rows = 0).
    Psi = zeros(4*nV, nCh);
    Psi(2:4:end, :) = G(:, 1:3:end).';   % x -> i
    Psi(3:4:end, :) = G(:, 2:3:end).';   % y -> j
    Psi(4:4:end, :) = G(:, 3:3:end).';   % z -> k

    Gm = Psi' * (dbasis.Mass * dbasis.Phi);   % [nCh x nModes]
end

% Author: Diellor Basha, 2026
