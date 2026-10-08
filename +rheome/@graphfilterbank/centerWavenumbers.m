function k = centerWavenumbers(obj)
% CENTERWAVENUMBERS  Gain-weighted centroid of sqrt(lambda) per member.
%
%   k = centerWavenumbers(gfb)      [1 x M], rad/m
%
% The centerFrequencies analogue -- but a WAVENUMBER, not a frequency: for a geometric
% Laplacian lambda has units 1/m^2, so sqrt(lambda) is rad/m and "frequency" would
% import the wrong dimension.
%
% A DISPLAY summary only: it is contaminated by where the spectrum was truncated, so it
% is not the scale a member is matched to. Report widths(gfb) for that.
%
% See also: wavelengths, widths, qfactor
%
% Author: Diellor Basha, 2026

    % The integration floor follows rheome.filters.frame: the smallest NONZERO eigenvalue when
    % a spectrum or explicit range was supplied, else the LPFactor fallback.
    if ~isempty(obj.Lmin_), lo = max(obj.Lmin_, eps);
    else,                   lo = max(obj.Lmax_/obj.LPFactor, eps);
    end
    lg = linspace(lo, obj.Lmax_, 512)';
    M  = obj.NumMembers;
    k  = zeros(1, M);
    for m = 1:M
        gm = abs(obj.G_{m}(lg));
        if sum(gm) > 0, k(m) = sqrt(sum(lg .* gm) / sum(gm)); end
    end
end

% Author: Diellor Basha, 2026
