function P = movingpeak(atoms, cums, opts)
% FLOW.MOVINGPEAK  A scalar peak carried along a path at a constant speed: the planted ground truth.
%
%   P = rheome.flow.movingpeak(atoms, cums, SpeedMS=0.25, SampleRate=300)
%   P = rheome.flow.movingpeak(atoms, cums, SpeedMS=0, StillSec=2)          % the stationary control
%
% The scalar sibling of rheome.flow.movingvortex. `atoms(:,k)` is the map with its peak at waypoint k --
% typically impulse(gfb, m, path(k)) of a Laplace-Beltrami mexhat, which has the positive core and
% negative ring of a curl map around a vortex -- and `cums(k)` is waypoint k's arclength along the
% path. At time t the peak sits at arclength s = SpeedMS * t, and the map is the linear blend of the
% two atoms either side of s, so the core moves continuously rather than hopping a waypoint at a time.
%
% ⭐ BUILD THE PATH ON rheome.geom.edgegraph, NOT ON rheome.geom.geodesic. Shortest paths on edgegraph give `cums`
% on a ruler that is symmetric and +1.4% against the analytic geodesic. Choosing waypoints by the heat
% method -- as rheome.flow.movingvortex does -- switches between near-equal routes around a sulcus, and the
% planted peak teleports (measured: 46.7 mm steps on a 3 mm spacing; see plant_track_omega.m).
%
% INPUTS
%   atoms      [nV x K] one map per waypoint, peak normalised to 1 by the caller if wanted
%   cums       [K x 1] arclength of each waypoint (m), increasing, cums(1) = 0
%   SpeedMS    arclength speed (m/s); 0 holds the peak at waypoint 1 for StillSec
%   SampleRate frames per second (300 = 30.9 per cycle at a 9.71 Hz alpha peak)
%   StillSec   record length when SpeedMS = 0 (2)
%
% OUTPUT (struct P)
%   .X [nV x nT] the maps     .t [1 x nT] s     .s [nT x 1] arclength of the peak (m)
%   .k [nT x 1] nearest waypoint to the peak   .k0, .w  the blend: waypoint k0 with weight 1-w
%
% See also: rheome.flow.movingvortex, rheome.geom.edgegraph, rheome.graphfilterbank/impulse, plant_tiletrack_omega
%
% Author: Diellor Basha, 2026

    arguments
        atoms double
        cums (:,1) double
        opts.SpeedMS (1,1) double {mustBeNonnegative} = 0.25
        opts.SampleRate (1,1) double {mustBePositive} = 300
        opts.StillSec (1,1) double {mustBePositive} = 2
    end
    K = size(atoms, 2);
    if numel(cums) ~= K, error('flow:movingpeak:size', 'cums has %d entries for %d atoms.', numel(cums), K); end
    if K < 2 || any(diff(cums) <= 0)
        error('flow:movingpeak:path', 'Need >= 2 waypoints with strictly increasing arclength.');
    end
    v = opts.SpeedMS;  fs = opts.SampleRate;
    if v == 0, nT = round(opts.StillSec * fs);  s = zeros(nT, 1);
    else,      nT = floor(cums(end) / v * fs) + 1;  s = (0:nT-1)' / fs * v;
    end
    k0 = discretize(s, [cums(1:end-1); inf]);  k0(isnan(k0)) = K - 1;  k0 = min(k0, K - 1);
    w = (s - cums(k0)) ./ (cums(k0 + 1) - cums(k0));
    P.X = atoms(:, k0) .* (1 - w') + atoms(:, k0 + 1) .* w';
    P.t = (0:nT-1) / fs;
    P.s = s;  P.k0 = k0;  P.w = w;
    P.k = k0 + (w > 0.5);
end

% Author: Diellor Basha, 2026
