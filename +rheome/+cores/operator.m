function H = operator(P, s, Z, g)
% CORES.OPERATOR  The core operator P'*((g/s).*Z) as nxr builds it, in single.
%
%   H = rheome.cores.operator(P, s, Z, g)   % P [K x V], s [K x 1] mode norms, Z [K x C] curl in modes, g [K x 1] gain
%   H = rheome.cores.operator(P, s, Z, rheome.cores.heatgain(lam))                    % curl smoothed, sigma 10 mm
%   H = rheome.cores.operator(P, s, Z, rheome.cores.hatgain(lam, sum(mass), level))   % curl at a ladder level
%
% H [V x C] is single, as nxr casts it; everything after is single (compare at ~1e-6 relative).
%
% Ports nxr-analysis src/nxr/vortex.py core_operator and src/nxr/vortex_cores.py hat_operator (commit 903fdc2).
%
% See also: rheome.cores.heatgain, rheome.cores.hatgain, rheome.cores.modecoords, rheome.cores.threshold
%
% Author: Diellor Basha, 2026
    H = single(P.' * ((g(:) ./ s(:)) .* Z));
end
