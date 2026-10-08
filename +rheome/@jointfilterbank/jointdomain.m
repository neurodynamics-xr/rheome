function d = jointdomain(obj)
% JOINTDOMAIN  This bank's joint plane as a domain descriptor: (eigenmode x frequency).
%
%   d = jointdomain(jfb)
%   [ok, why] = rheome.domain.same(jointdomain(a), jointdomain(b))
%
% ⭐ The bank has always known it lives on a product -- `axes` returns both axes and iscompatible
% checks them -- but as a class-local struct, so "the same plane" meant something different here than
% in +pipeline or in the tile store. This returns the shared type instead, with .half and .boundary
% carried as CONVENTIONS so rheome.domain.same compares them.
%
% ⚠⚠ THIS DOES NOT REPLACE iscompatible, and the reason is structural rather than a shortcut. A domain
% "carries no coordinates: two fields are on the same domain when the ids match" (+domain/Contents.m),
% so a descriptor holds COUNTS and conventions and cannot hold eigenvalues. iscompatible compares the
% lambda and f VALUES with a tolerance -- catching the same mode count over different anatomy, which no
% id comparison can see. Use rheome.domain.same for the counts and the conventions, iscompatible for the
% values, and do not expect one to subsume the other.
%
% See also: axes, iscompatible, rheome.domain.product, rheome.domain.same
%
% Author: Diellor Basha, 2026

    ax = axes(obj);
    dM = rheome.domain.index(numel(ax.lambda), Name="eigenmode");
    dF = rheome.domain.index(numel(ax.f),      Name="frequency", Units="Hz");
    d  = rheome.domain.product(dM, dF, Name="eigenmode x frequency", ...
             Conventions=struct('half', char(ax.half), 'boundary', char(ax.boundary), ...
                                'fs', ax.fs));
end

% Author: Diellor Basha, 2026
