function f = centerFrequencies(obj)
% CENTERFREQUENCIES  Energy-weighted centroid of temporal frequency per member.
%
%   f = centerFrequencies(jfb)      [1 x Nf], Hz
%
% ⭐ This is the one cwtfilterbank name that transfers unchanged AND means what it says:
% on the TIME axis omega really is a frequency, so Hz is correct here even though
% "frequency" was wrong for lambda.
%
% See also: centerWavenumbers
%
% Author: Diellor Basha, 2026

    lam = obj.Lambda;  om = obj.Omega;  ff = obj.Frequencies;
    f = zeros(1, obj.NumMembers);
    for m = 1:obj.NumMembers
        P = abs(jfb_member(obj, m, lam, om)).^2;
        s = sum(P(:));
        if s > 0, f(m) = sum(ff .* sum(P, 1)) / s; end
    end
end

% Author: Diellor Basha, 2026
