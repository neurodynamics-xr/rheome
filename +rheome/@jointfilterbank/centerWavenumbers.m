function k = centerWavenumbers(obj)
% CENTERWAVENUMBERS  Energy-weighted centroid of sqrt(lambda) per member.
%
%   k = centerWavenumbers(jfb)      [1 x Nf], rad/m
%
% A WAVENUMBER, not a frequency: for a geometric Laplacian lambda has units 1/m^2, so
% sqrt(lambda) is rad/m. Marginalises the member's full 2-D gain, so it is defined for
% non-separable kernels too.
%
% See also: centerFrequencies
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  om = obj.Omega;
    kk  = sqrt(lam);
    k   = zeros(1, obj.NumMembers);
    for m = 1:obj.NumMembers
        P = abs(jfb_member(obj, m, lam, om)).^2;
        s = sum(P(:));
        if s > 0, k(m) = sum(kk .* sum(P, 2)) / s; end
    end
end

% Author: Diellor Basha, 2026
