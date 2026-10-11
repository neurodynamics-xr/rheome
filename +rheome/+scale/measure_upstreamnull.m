function [T, X] = measure_upstreamnull(name, S, opts)
% SCALE.MEASURE_UPSTREAMNULL  Gate N calibration in one head: does the lagged-source null read planted non-travel as non-travel and keep planted travel?
%
%   [T, X] = rheome.scale.measure_upstreamnull(name)
%   [T, X] = rheome.scale.measure_upstreamnull(name, S, Gradient=f, SnrDB=[0 10], Hz=[2 5], SpeedMS=[0.5 1 2])
%
% Why: U* = U - U_surrogate (rheome.scale.measure_upstream, gate N) is a claim about travel only if the null
% (rheome.forward.laggedrest) cancels activity that does NOT travel and leaves activity that does. This plants
% both in the participant's own head and resting recording, and reads every plant with the frozen readout v1.
%
% PLANTS (normal current, area-weighted, through S.G, ADDED to the recording's selected rows S.iSel):
%   travel   a plane wave a(t) cos(2 pi f t - sigma 2 pi x / lambda) over the whole cortex, x the vertex position
%            along Axis, lambda = SpeedMS / Hz, sigma = +1 (travel along +Axis) or -1
%   lagged   per hemisphere two discs of each DiscMM radius that do not move, one in the anterior and one in the
%            posterior quarter along Axis (random centres, Locations draws): a(t) cos(2 pi f t) in front and
%            Ratio a(t) cos(2 pi f t - delta) behind, delta in LagRad; dir = -sign(delta) (the back leads:
%            the lag looks like travel along +Axis)
%   none     the recording alone, its null drawn with NullSeeds seeds (the null's own seed spread)
%   a(t) = exp(0.3 z), z a 2 s moving mean of white noise at unit SD (slow amplitude drift, as at rest).
%   SnrDB  = the plant's mean sensor power over the recording's own 1-8 Hz sensor power.
%
% OUTPUT
%   X  one row per plant: kind, hz, speed_ms, lambda_mm, dir, lag_rad, disc_mm, snr_db, seed, U, U_surrogate, U_star
%   T  rheome.scale.rows("upstreamnull"), per SnrDB (band "snr<dB>"): null_ustar_absmax (largest |U*| of a
%      lagged plant: the calibrated threshold), null_ustar_absmedian, lagged_dev_median (median |U - 0.5| of the
%      lagged plants: what the null had to cancel), seed_ustar_sd (the "none" seeds); and per lambda (band
%      "snr<dB>_lambda<mm>"): travel_power (share of travel plants with |U*| above the threshold and U* of
%      U - 0.5's sign), travel_dir_sign (share whose U - 0.5 has the sign sigma: only meaningful when g
%      increases along Axis)
%
% ⚠ THE READOUT ALIASES PAST HALF A CYCLE between tile centres (~120 mm): lambda < ~240 mm can read backwards
%   whatever the null does; travel_dir_sign shows it, travel_power does not (it compares U* with U only).
% ⚠ A plane wave and two discs are the simplest movers and non-movers: they calibrate the null on the head and
%   the readout, not on real resting activity's geometry.
%
% See also: rheome.scale.measure_upstream, rheome.forward.laggedrest, rheome.flow.upstreamshare
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Gradient char = ''
        opts.SignCheck (1,1) logical = true
        opts.SnrDB (1,:) double = [0 10]
        opts.Hz (1,:) double = [2 5]
        opts.SpeedMS (1,:) double = [0.5 1 2]
        opts.LagRad (1,:) double = [-pi/2 -pi/4 pi/4 pi/2]
        opts.Ratio (1,1) double = 0.8
        opts.DiscMM (1,:) double = [20 40]
        opts.Locations (1,1) double {mustBeInteger, mustBePositive} = 2
        opts.NullSeeds (1,1) double {mustBeInteger, mustBePositive} = 3
        opts.Axis (1,3) double = [1 0 0]                 % SCS: +x anterior
        opts.MinSpanS (1,1) double = 180
        opts.OwnMin (1,1) double = 0.25
        opts.Seed (1,1) double = 1
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    nT = size(F, 2);  t = (0:nT-1) / fs;
    N = double(S.Sf.VertNormals);  N = N ./ vecnorm(N, 2, 2);
    Gn = N(:,1)' .* S.G(:, 1:3:end) + N(:,2)' .* S.G(:, 2:3:end) + N(:,3)' .* S.G(:, 3:3:end);   % [nSel x nV]
    w = zeros(S.nV, 1);  for h = ["L" "R"], w(S.B.(h).gv) = full(sum(S.B.(h).lbo.Mass, 2)); end
    x = double(S.Sf.Vertices) * (opts.Axis(:) / norm(opts.Axis));     % metres along Axis
    C = fft(F, [], 2);  f = (0:nT-1) * fs / nT;  C(:, ~(f >= 1 & f <= 8 | f >= fs-8 & f <= fs-1)) = 0;
    pBand = mean(real(ifft(C, [], 2)).^2, 'all');  clear C
    rng(opts.Seed, 'twister');
    a = movmean(randn(1, nT), round(2 * fs));  a = exp(0.3 * a / std(a));

    P = struct('kind', {}, 'hz', {}, 'speed_ms', {}, 'lambda_mm', {}, 'dir', {}, 'lag_rad', {}, 'disc_mm', {}, 'seed', {}, 'gc', {}, 'gs', {});
    for hz = opts.Hz
        for c = opts.SpeedMS
            lam = c / hz;
            for sg = [1 -1]
                k = sg * 2 * pi / lam;
                P(end+1) = struct('kind', "travel", 'hz', hz, 'speed_ms', c, 'lambda_mm', 1000 * lam, 'dir', sg, ...
                    'lag_rad', NaN, 'disc_mm', NaN, 'seed', 0, 'gc', Gn * (w .* cos(k * x)), 'gs', Gn * (w .* sin(k * x))); %#ok<AGROW>
            end
        end
        for r = opts.DiscMM
            for loc = 1:opts.Locations
                [Df, Db] = i_discs(S, x, r / 1000);
                for d = opts.LagRad
                    P(end+1) = struct('kind', "lagged", 'hz', hz, 'speed_ms', NaN, 'lambda_mm', NaN, 'dir', -sign(d), ...
                        'lag_rad', d, 'disc_mm', r, 'seed', loc, 'gc', Gn * (w .* (Df + opts.Ratio * cos(d) * Db)), ...
                        'gs', Gn * (w .* (opts.Ratio * sin(d) * Db))); %#ok<AGROW>
                end
            end
        end
    end

    up = @(Y, seed) i_read(name, S, Y, seed, opts);
    X = table();
    for s = 1:opts.NullSeeds
        X = [X; i_row("none", NaN, NaN, NaN, 0, NaN, NaN, NaN, s, up(F, opts.Seed + 100 * s))]; %#ok<AGROW>
    end
    for snr = opts.SnrDB
        for i = 1:numel(P)
            p = P(i);  Y = a .* (p.gc * cos(2*pi*p.hz*t) + p.gs * sin(2*pi*p.hz*t));
            Y = F + Y * sqrt(10^(snr/10) * pBand / mean(Y.^2, 'all'));
            X = [X; i_row(p.kind, p.hz, p.speed_ms, p.lambda_mm, p.dir, p.lag_rad, p.disc_mm, snr, p.seed, up(Y, opts.Seed))]; %#ok<AGROW>
            fprintf('[upstreamnull %s] %d/%d %s %g Hz %g dB: U %.3f U* %.3f\n', name, i, numel(P), p.kind, p.hz, snr, X.U(end), X.U_star(end));
        end
    end
    T = i_summary(X, opts.SnrDB);
end

function u = i_read(name, S, Y, seed, opts)
    Tu = rheome.scale.measure_upstream(name, S, Data=Y, Passes=["periodic" "surrogate"], Halves=false, ...
        Gradient=opts.Gradient, SignCheck=opts.SignCheck, MinSpanS=opts.MinSpanS, OwnMin=opts.OwnMin, Seed=seed);
    v = @(m) [Tu.value(Tu.metric == m); NaN];
    u = [v("U") v("U_surrogate") v("U_star")];  u = u(1, :);           % NaN when the recording is not read
end

function [Df, Db] = i_discs(S, x, r)
% one anterior-quarter and one posterior-quarter disc per hemisphere, at random centres
    V = double(S.Sf.Vertices);  Df = zeros(S.nV, 1);  Db = Df;
    for h = ["L" "R"]
        gv = S.B.(h).gv;  q = quantile(x(gv), [0.25 0.75]);
        fr = gv(x(gv) >= q(2));  bk = gv(x(gv) <= q(1));
        cf = fr(randi(numel(fr)));  cb = bk(randi(numel(bk)));
        Df(gv(vecnorm(V(gv,:) - V(cf,:), 2, 2) <= r)) = 1;  Db(gv(vecnorm(V(gv,:) - V(cb,:), 2, 2) <= r)) = 1;
    end
end

function R = i_row(kind, hz, c, lam, d, lag, disc, snr, seed, u)
    R = table(string(kind), hz, c, lam, d, lag, disc, snr, seed, u(1), u(2), u(3), 'VariableNames', ...
        {'kind','hz','speed_ms','lambda_mm','dir','lag_rad','disc_mm','snr_db','seed','U','U_surrogate','U_star'});
end

function T = i_summary(X, snrs)
    T = rheome.scale.rows("upstreamnull", "seed_ustar_sd", std(X.U_star(X.kind == "none")), "share");
    for snr = snrs
        b = sprintf("snr%g", snr);  L = X(X.kind == "lagged" & X.snr_db == snr, :);
        thr = max(abs(L.U_star));
        T = [T; rheome.scale.rows("upstreamnull", ["null_ustar_absmax" "null_ustar_absmedian" "lagged_dev_median"], ...
            [thr median(abs(L.U_star)) median(abs(L.U - 0.5))], "share", b)]; %#ok<AGROW>
        Tr = X(X.kind == "travel" & X.snr_db == snr, :);
        for lam = unique(Tr.lambda_mm)'
            r = Tr(Tr.lambda_mm == lam, :);
            pw = mean(abs(r.U_star) > thr & sign(r.U_star) == sign(r.U - 0.5));
            T = [T; rheome.scale.rows("upstreamnull", ["travel_power" "travel_dir_sign"], ...
                [pw mean(sign(r.U - 0.5) == r.dir)], "fraction", sprintf("%s_lambda%g", b, round(lam)))]; %#ok<AGROW>
        end
    end
end

% Author: Diellor Basha, 2026
