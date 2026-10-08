function E = scalogram(obj)
% SCALOGRAM  Energy per spatial scale at every frame -- the cortical scalogram.
%
%   E = scalogram(fp)      -> [nScale x nFrames]
%
% ⭐ THE VIEW THE HYPOTHESIS IS READ FROM. With scale on one axis and time on the other,
% a pattern breaking into smaller vortices appears directly as energy MIGRATING UP the
% scale axis -- not as a shift in a summary statistic that has to be trusted. The centroid
% is a one-number projection of this; the image is the thing itself.
%
% Computed in the coefficient domain, one pass per scale over C [Ks x nT], so no vertex
% field is formed: ~37 M operations for a 7-scale, 6500-frame page.
%
% ⚠ RAW ENERGY, NOT NORMALISED. Column-normalising (each frame's scales summing to 1) is
% what separates a change of SCALE from a change of POWER, and is usually what you want to
% look at -- but it is a display choice, so it is made by the caller and not baked in here.
%
% See also: scaleEnergy, centroid, rheome.flowbrowser/drawScalogram
%
% Author: Diellor Basha, 2026

    C  = double(obj.Coefficients);                 % [Ks x nT]
    nS = obj.NumScales;
    E  = zeros(nS, size(C, 2));
    for g = 1:nS
        E(g, :) = sum(abs(obj.ScaleGains(:, g) .* C).^2, 1);
    end
end

% Author: Diellor Basha, 2026
