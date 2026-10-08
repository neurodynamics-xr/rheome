function bw = powerbw(obj)
% POWERBW  Half-power WAVENUMBER band of each member.
%
%   bw = powerbw(gfb)     [M x 2], columns [k_lo k_hi] in rad/m
%
% "Bandwidth" is used generically on any axis (cf. spatial bandwidth in optics); the
% axis here is the wavenumber, not a frequency.
%
% See also: qfactor, centerWavenumbers
%
% Author: Diellor Basha, 2026

    lam = linspace(0, obj.Lmax_, 4096)';
    H   = graphfilters(obj, 'Lambda', lam);
    k   = sqrt(lam);
    M   = obj.NumMembers;
    bw  = nan(M, 2);
    for m = 1:M
        p = H(:,m).^2;
        [pk, ipk] = max(p);
        if pk <= 0, continue; end
        lo = find(p(1:ipk) >= pk/2, 1, 'first');
        hi = ipk - 1 + find(p(ipk:end) >= pk/2, 1, 'last');
        bw(m,:) = [k(lo), k(hi)];
    end
end

% Author: Diellor Basha, 2026
