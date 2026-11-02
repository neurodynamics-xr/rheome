function [Z, s] = modecoords(P, mass, X)
% CORES.MODECOORDS  A vertex operator in mass-orthonormal mode coordinates, Z = P*(mass.*X)./s (nxr Cortex.coords).
%
%   [Z, s] = rheome.cores.modecoords(P, mass, X)   % P [K x V] modes as rows, mass [V x 1] (mm^2), X [V x C]
%
% s = sqrt(sum(P.^2 .* mass', 2)) is each mode's norm (K x 1), as nxr Cortex.s.
%
% Ports nxr-analysis src/nxr/vortex.py Cortex.__init__, coords (commit 903fdc2).
%
% See also: rheome.cores.operator
%
% Author: Diellor Basha, 2026
    mass = mass(:);
    s = sqrt(sum(P.^2 .* mass.', 2));
    Z = (P * (mass .* X)) ./ s;
end
