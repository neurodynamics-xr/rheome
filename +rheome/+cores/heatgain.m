function g = heatgain(lam, sigmaMM)
% CORES.HEATGAIN  Heat-kernel gain exp(-lambda*sigma^2/2) on Laplace-Beltrami eigenvalues (nxr core_operator).
%
%   g = rheome.cores.heatgain(lam)            % sigma = 10 mm (nxr CORE_SIGMA)
%   g = rheome.cores.heatgain(lam, sigmaMM)   % lam in 1/mm^2
%
% Ports nxr-analysis src/nxr/vortex.py core_operator (commit 903fdc2).
%
% See also: rheome.cores.operator, rheome.cores.hatgain
%
% Author: Diellor Basha, 2026
    if nargin < 2 || isempty(sigmaMM), sigmaMM = 10; end
    g = exp(-lam(:) * sigmaMM^2 / 2);
end
