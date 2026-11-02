function g = hatgain(lam, areaTotal, level)
% CORES.HATGAIN  Level hat gain lambda*s2*exp(-lambda*s2/2), scaled to peak 1 (nxr hat_operator).
%
%   g = rheome.cores.hatgain(lam, areaTotal, level)
%
% Peak at lambda_peak = 4*pi*2^(level-1/2) / (areaTotal/2), the geometric centre of the ladder level's Weyl band, s2 = 2/lambda_peak.
% areaTotal is the summed vertex area (mm^2) of the mesh; the formula assumes two hemispheres (so half of it is one
% hemisphere), as nxr does. lam in 1/mm^2. This is a new gain, distinct from rheome.filters Mexican hats.
%
% Ports nxr-analysis src/nxr/vortex_cores.py hat_operator (commit 903fdc2).
%
% See also: rheome.cores.operator, rheome.cores.heatgain
%
% Author: Diellor Basha, 2026
    lamPeak = 4*pi * 2^(level - 0.5) / (areaTotal / 2);
    s2 = 2 / lamPeak;
    g = lam(:) * s2 .* exp(-lam(:) * s2 / 2);
    g = g / max(g);
end
