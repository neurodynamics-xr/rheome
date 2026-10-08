function R = onsetspread(Y, G, fs, seed, tPeak, opts)
% DETECT.ONSETSPREAD  How a burst's ONSET spreads across tiles: slowness from half-rise times.
%
%   R = rheome.detect.onsetspread(Y, G, fs, seed, tPeak)
%   R = rheome.detect.onsetspread(Y, G, fs, seed, tPeak, Radius=0.12, Pre=1, Post=0.5, Centres=xyz)
%
% Y is a tile envelope [nTile x nT] at fs (rheome.geom.tilemean of |J^|), G is rheome.geom.tiles WITH a Ruler (.D),
% and (seed, tPeak) is one burst: its tile and the time (s) it peaks there. Each tile within Radius of
% the seed that takes part in the burst gets an onset time, and the onsets are regressed on distance.
%
% ⭐ ONSET IS THE HALF-RISE TIME, NOT A THRESHOLD CROSSING. A fixed threshold is crossed EARLIER by a
% louder tile even when every tile rises at the same instant, which manufactures propagation from
% strong tiles to weak ones. t50 is when a tile's envelope reaches halfway from its own pre-peak floor
% to its own peak inside the window: invariant to each tile's amplitude, so a simultaneous rise reads
% as simultaneous whatever the amplitudes.
%
% ⭐ SLOWNESS, NOT SPEED. The slope of t50 on geodesic distance from the seed is in s/m: 0 is a
% simultaneous rise (standing, or faster than the timing resolves), positive spreads OUTWARD from the
% seed, negative converges on it. Speed is 1/slowness and is infinite exactly where nothing propagates,
% so it is the wrong quantity to average. The origin is the burst's PEAK tile, not its earliest one:
% choosing the earliest tile as origin would make every slope positive by selection.
%
% A planar fit t50 = a + s . x, with x the participating tile centres in their own two principal
% directions (Centres required), gives a slowness VECTOR: |s| and a direction, for a front that
% crosses the region rather than spreading from its peak.
%
% INPUTS
%   Y        [nTile x nT] envelope     G  rheome.geom.tiles with Ruler (.D; .centre for Centres)
%   fs       frames per second         seed, tPeak  the event
%   Radius   tiles within this geodesic distance of the seed (0.12 m)
%   Pre/Post window before / after tPeak searched for each tile's rise and peak (1, 0.5 s)
%   Smooth   moving mean applied first, symmetric so it shifts no tile (0.1 s)
%   MinRise  a tile participates if its rise is >= MinRise x the seed's (0.5)
%   MinTiles fewer participating tiles -> no fit (5)
%   Centres  [nTile x 3] tile centre coordinates, for the planar fit ([] = skip)
%
% OUTPUT (struct R)
%   .tiles .t50 (s, relative to tPeak) .dist (m) .rise   participating tiles
%   .slowness (s/m, Theil-Sen)  .slownessLS  .r2  .spreadS (range of t50)  .n
%   .planar  slowness vector magnitude (s/m)  .planarDir  unit direction  .planarR2
%   .ok      false if too few tiles took part
%
% See also: rheome.detect.spindle, rheome.geom.tiles, alpha_burstspread_omega
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        G struct
        fs (1,1) double {mustBePositive}
        seed (1,1) double {mustBeInteger, mustBePositive}
        tPeak (1,1) double
        opts.Radius (1,1) double {mustBePositive} = 0.12
        opts.Pre (1,1) double {mustBePositive} = 1
        opts.Post (1,1) double {mustBePositive} = 0.5
        opts.Smooth (1,1) double {mustBeNonnegative} = 0.1
        opts.MinRise (1,1) double {mustBeNonnegative} = 0.5
        opts.MinTiles (1,1) double {mustBePositive} = 5
        opts.Centres = []
    end
    R = struct('tiles',[],'t50',[],'dist',[],'rise',[],'slowness',NaN,'slownessLS',NaN,'r2',NaN, ...
        'spreadS',NaN,'n',0,'planar',NaN,'planarDir',[NaN NaN],'planarR2',NaN,'ok',false);
    i0 = round((tPeak - opts.Pre) * fs) + 1;  i1 = round((tPeak + opts.Post) * fs) + 1;
    if i0 < 1 || i1 > size(Y, 2), return; end
    k = max(1, round(opts.Smooth * fs));
    cand = find(G.D(seed, :) <= opts.Radius);
    Ys = movmean(Y(cand, i0:i1), k, 2);
    tt = ((i0:i1) - 1) / fs - tPeak;
    t50 = nan(numel(cand), 1);  rise = zeros(numel(cand), 1);
    for c = 1:numel(cand)
        s = Ys(c, :);  [pk, ip] = max(s);
        if ip < 2, continue; end
        [base, ib] = min(s(1:ip));  rise(c) = pk - base;
        lv = base + 0.5 * rise(c);
        j = find(s(ib:ip) >= lv, 1) + ib - 1;                    % first reach of half-rise
        if isempty(j) || j <= 1, continue; end
        t50(c) = tt(j-1) + (lv - s(j-1)) / (s(j) - s(j-1)) * (tt(j) - tt(j-1));
    end
    rs = rise(cand == seed);
    part = isfinite(t50) & rise >= opts.MinRise * rs;
    if sum(part) < opts.MinTiles, return; end
    tl = cand(part);  t = t50(part);  d = G.D(seed, tl)';
    R.tiles = tl(:);  R.t50 = t;  R.dist = d;  R.rise = rise(part);  R.n = numel(tl);
    R.spreadS = max(t) - min(t);
    % Theil-Sen slope over pairs at different distances
    [a, b] = find(triu(true(numel(d)), 1));
    dd = d(b) - d(a);  ok = abs(dd) > 1e-6;
    R.slowness = median((t(b(ok)) - t(a(ok))) ./ dd(ok));
    p = polyfit(d, t, 1);  R.slownessLS = p(1);
    R.r2 = 1 - sum((t - polyval(p, d)).^2) / max(sum((t - mean(t)).^2), realmin);
    if ~isempty(opts.Centres)
        X = opts.Centres(tl, :);  X = X - mean(X, 1);
        [~, ~, V] = svd(X, 'econ');  Z = X * V(:, 1:2);
        A = [ones(numel(t), 1) Z];  cf = A \ t;
        R.planar = norm(cf(2:3));  R.planarDir = (cf(2:3) / max(norm(cf(2:3)), realmin))';
        R.planarR2 = 1 - sum((t - A * cf).^2) / max(sum((t - mean(t)).^2), realmin);
    end
    R.ok = true;
end

% Author: Diellor Basha, 2026
