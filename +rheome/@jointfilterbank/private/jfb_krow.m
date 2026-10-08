function r = jfb_krow(k, kWanted)
% JFB_KROW  Map a wavenumber to its ROW on a spectralResponse image.
%
%   r = jfb_krow(k, kWanted)
%
% sqrt(lambda) is non-uniformly spaced, so spectralResponse plots rows by INDEX with the
% wavenumbers as tick labels. Anything overlaid on that image must be mapped the same way,
% or it lands somewhere else entirely.
%
% ⚠ NOT interp1. A symmetric domain has DEGENERATE MULTIPLETS -- on a sphere degree l has
% multiplicity 2l+1, so sqrt(lambda) repeats and interp1 refuses non-unique sample points.
% This takes the CENTRE of the matching multiplet instead, which is also where the energy
% visually sits.
%
% Author: Diellor Basha, 2026

    k = double(k(:));
    r = zeros(size(kWanted));
    for i = 1:numel(kWanted)
        d = abs(k - double(kWanted(i)));
        r(i) = mean(find(d <= min(d) + 1e-9*max(1, max(k))));
    end
end

% Author: Diellor Basha, 2026
