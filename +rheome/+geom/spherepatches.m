function P = spherepatches(Sph, opts)
% GEOM.SPHEREPATCHES  Shared tiles on the registered sphere, and where the polar gauge breaks.
%
%   P = rheome.geom.spherepatches(Sph)                       % ico3: 642 patches per hemisphere
%   P = rheome.geom.spherepatches(Sph, Level=2, Axis=[1 0 0], MaxSpread=15)
%
% Every vertex goes to the nearest vertex of an ico-Level sphere (rheome.geom.icosphere) on its
% REGISTERED sphere point, so patch k is the same anatomical place in every subject: a sum over
% a patch is comparable across a cohort without resampling anything. Then, per patch, how far the
% rheome.geom.sphereframe gauge turns across it:
%
%   spread = max over member vertices of |angle between the member's meridian, parallel-transported
%            along the great circle to the patch centre, and the centre's own meridian|
%
% It is ~0 on the equator and grows as 2*rho*cot(colatitude) towards a pole (rho = patch radius),
% so it is THE number that says whether a pooled (north, west) vector means anything. A patch is
% .excluded when spread > MaxSpread or it holds a pole vertex; read those in the second chart
% (Axis=[1 0 0]), where they lie near the equator.
%
% INPUTS:  Sph [nV x 3] registered-sphere points (any radius)
%   Level      ico subdivision level, default 3 (642 patches; ~8 deg apart, ~14 mm on a 100 mm sphere)
%   Axis       pole axis of the gauge being judged, default [0 0 1]
%   MaxSpread  degrees, default 15
%
% OUTPUT (struct P):
%   .patch [nV x 1]  patch of each vertex        .n [nP x 1] vertices per patch
%   .centre [nP x 3] unit centres                .colat .lon [nP x 1] of the centres, degrees
%   .spread [nP x 1] degrees (Inf for an empty patch)   .hasPole [nP x 1]   .excluded [nP x 1]
%
% See also: rheome.geom.sphereframe, rheome.geom.icosphere, rheome.geom.tilemean
%
% Author: Diellor Basha, 2026

    arguments
        Sph (:,3) double
        opts.Level (1,1) double {mustBeInteger, mustBeNonnegative} = 3
        opts.Axis (1,3) double = [0 0 1]
        opts.MaxSpread (1,1) double = 15
    end
    ax = opts.Axis / norm(opts.Axis);
    u = Sph ./ vecnorm(Sph,2,2);
    [C, ~] = rheome.geom.icosphere(opts.Level);  C = C ./ vecnorm(C,2,2);  nP = size(C,1);
    [~, k] = max(u*C', [], 2);
    n = accumarray(k, 1, [nP 1]);

    mer = @(x) (ax - (x*ax').*x) ./ max(vecnorm(ax - (x*ax').*x, 2, 2), eps);
    p = C(k,:);  ei = mer(u);  ep = mer(p);
    % transport ei from u to p along the great circle (Rodrigues about u x p)
    w = cross(u, p, 2);  sn = vecnorm(w,2,2);  cs = sum(u.*p,2);
    kk = w ./ max(sn, eps);
    Re = ei.*cs + cross(kk, ei, 2).*sn + kk.*sum(kk.*ei,2).*(1-cs);
    d = abs(atan2d(sum(p.*cross(ep, Re, 2),2), sum(ep.*Re,2)));
    [~, north] = max(u*ax');  [~, south] = min(u*ax');
    hasPole = false(nP,1);  hasPole(k([north south])) = true;
    spread = accumarray(k, d, [nP 1], @max, Inf);
    colat = acosd(max(-1, min(1, C*ax')));
    if abs(ax(3)) > 0.9, a1 = [1 0 0]; else, a1 = [0 0 1]; end
    a1 = a1 - (a1*ax')*ax;  a1 = a1/norm(a1);  a2 = cross(ax, a1);
    P = struct('patch', k, 'n', n, 'centre', C, 'colat', colat, 'lon', atan2d(C*a2', C*a1'), ...
               'spread', spread, 'hasPole', hasPole, ...
               'excluded', hasPole | ~(spread <= opts.MaxSpread), 'axis', ax, 'level', opts.Level);
end

% Author: Diellor Basha, 2026
