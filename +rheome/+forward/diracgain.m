function [Gmode, Psi] = diracgain(G, dbasis)
% FORWARD.DIRACGAIN  The SYNTHESIS gain: one sensor pattern per Dirac eigenmode, Gmode = G * Phi.
%
%   [Gmode, Psi] = rheome.forward.diracgain(G, dbasis)
%
% ⭐ WHAT THIS IS FOR. The leadfield answers "what sensor pattern does a source delta produce".
% Composed with the Dirac basis it answers "what sensor pattern does an EIGENMODE produce", which is
% the operator a mode-space simulation needs: seeding becomes a choice of coefficients C [nModes x nT]
% and the sensors are one GEMM, B = Gmode*C, with nothing of size [3nV x nT] ever built.
%
% ⚠⚠ THIS IS NOT rheome.forward.dirac, AND THE DIFFERENCE IS SILENT. rheome.forward.dirac is the ANALYSIS-convention
% gain -- what rheome.inverse.dirac consumes -- and carries a mass weighting: measured on a reference subject the
% two differ by a factor of about 1.2e5, which is 1/(mean vertex area) = 9.1e4 to within the spread of
% the areas. Substituting rheome.forward.dirac here type-checks, runs, returns tesla and reproduces the right
% spatial pattern 100+ dB too small. Nothing in the output looks wrong. ⭐ Use this to SYNTHESISE and
% rheome.forward.dirac to build an inverse; `rheome.operators.registry` now carries both rows with this caveat.
%
% INPUTS
%   G       [nCh x 3nV] vertex leadfield for the SAME vertices the basis covers (rheome.forward.leadfield
%           with GlobalVertices set)
%   dbasis  Dirac basis restricted to those vertices: .Phi [4nV x nModes], .nVert, .nModes
%
% OUTPUTS
%   Gmode  [nCh x nModes] tesla per unit coefficient   Psi  [3nV x nModes] the synthesis operator
%
% See also: rheome.forward.leadfield, rheome.forward.dirac, rheome.forward.reconstruct, rheome.forward.simulate
%
% Author: Diellor Basha, 2026

    nV = dbasis.nVert;  nM = dbasis.nModes;
    if size(G,2) ~= 3*nV
        error('forward:diracgain:shape', 'G has %d columns, expected 3*nVert = %d.', size(G,2), 3*nV);
    end
    % the current lives in the IMAGINARY part of the quaternion: rows 2,3,4 of each vertex block
    Psi = zeros(3*nV, nM);
    Psi(1:3:end, :) = dbasis.Phi((1:nV)*4 - 2, :);
    Psi(2:3:end, :) = dbasis.Phi((1:nV)*4 - 1, :);
    Psi(3:3:end, :) = dbasis.Phi((1:nV)*4,     :);
    Gmode = G * Psi;
end

% Author: Diellor Basha, 2026
