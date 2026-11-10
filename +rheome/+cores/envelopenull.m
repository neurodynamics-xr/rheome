function Y = envelopenull(xr, good, per, white)
% CORES.ENVELOPENULL  The lag-free surrogate rescaled frame by frame to the data's own envelope (lag-free, but it bursts when the data bursts).
%
%   Y = rheome.cores.envelopenull(xr, good, per, white)
%
% xr [C x m] real band signal (single or double), good [1 x m] logical, per frames per cycle, white [C x m] noise draw or a
% RandStream (see rheome.cores.lagfreenull). xn = lagfreenull(xr, good, white); sm(v) = sqrt(boxsmooth(mean over channels v^2, per))
% (zero-padded, one cycle); Y = single(xn .* (sm(xr) ./ max(sm(xn), 1e-30))), with sm(xn) taken on the single-precision xn. The
% rescaling happens before the smoothing, so Y's own smoothed envelope matches the data's to ~5 %, not exactly.
%
% Ports nxr-analysis src/nxr/vortex.py envelope_matched_null (line 794, commit 30b9622).
%
% See also: rheome.cores.lagfreenull, rheome.cores.boxsmooth
%
% Author: Diellor Basha, 2026
    arguments
        xr (:,:) {mustBeNumeric, mustBeReal}
        good logical
        per (1,1) double {mustBeInteger, mustBePositive}
        white
    end
    xn = rheome.cores.lagfreenull(xr, good, white);
    sm = @(v) sqrt(rheome.cores.boxsmooth(mean(double(v).^2, 1), per));
    Y = single(double(xn) .* (sm(xr) ./ max(sm(xn), 1e-30)));
end
