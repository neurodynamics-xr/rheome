function [T, X, Best] = measure_grouptrack(name, S, opts)
% SCALE.MEASURE_GROUPTRACK  Viterbi tile tracking of the alpha envelope, real against the mode-shift surrogate.
%
%   [T, X] = rheome.scale.measure_grouptrack(name)
%   [T, X] = rheome.scale.measure_grouptrack(name, S, Depths=7, JLevels=[0 2])
%
% The single-subject group-track analysis with the report printing and the blob-size aside removed; the same seed,
% the same surrogate, the same calibration. Per hemisphere: the 8-16 Hz analytic envelope of the
% minimum-norm current as 1000 LB mode coefficients at 300 frames/s (rheome.flow.envelopemodes, cached in
% rheome.load.root -- node-local on the cluster); a SURROGATE with every mode's time course circularly
% shifted by its own random offset (>= one 2 s window); per tile depth and dyadic frame rate,
% rheome.detect.tilepath (track-before-detect, Viterbi on the tile graph) in each 2 s window; the baseline
% is the BaseQ quantile of the surrogate, the score threshold is set on HALF the surrogate windows
% at Alpha false paths per second and the false rate tested on the other half.
%
% Returns rheome.scale.rows (analysis "grouptrack"), band = "<hemi>/d<depth>/<fps>fps", metrics:
%   tile_mm           median tile diameter                                    mm
%   false_rate_test   surrogate passes per second on the held-out half        1/s
%   frac_pass         fraction of 2 s windows with a path above threshold     fraction
%   dur_s             median duration of the passing paths                    s
%   net_mm            median net displacement per 2 s window (passing paths)  mm
%   speed_ms          median net speed (net / duration) ⭐                    m/s
%   straight          median straightness (net / path length)                 ratio
%   frac_beyond_tile  fraction of passing paths whose net exceeds one tile    fraction
%   sur_n_pass, sur_net_mm, sur_speed_ms   the same for the held-out surrogate passes
% plus "n_windows" (band ""). X is the per-configuration table (grouptrack.csv).
% Best (third output): per configuration, the highest-scoring REAL path -- .hemi .depth .jj .rate
% .window .score .passed .track (rheome.detect.tilepath track, frames within its window) .G (the tiles),
% which rheome.scale.measure_eventsensors maps onto the sensor data.
%
% ⚠ netMM is quantised to the tiling: displacement below one tile is not measured, by design.
% ⚠ Everything is through the instrument -- a trajectory is the trajectory of the RENDERED blob.
%
% See also: rheome.detect.tilepath, rheome.detect.trackstats, rheome.geom.tiles, rheome.flow.envelopemodes, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depths double = 6:8
        opts.JLevels double = [0 2 4 6]
        opts.Alpha (1,1) double = 0.1
        opts.BaseQ (1,1) double = 0.9
        opts.WindowS (1,1) double = 2
        opts.Rate (1,1) double = 300
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 17
    end
    rng(opts.Seed);
    if isempty(S), B = rheome.load.bases(name); else, B = S.B; end
    E = rheome.flow.envelopemodes(name, Band=[8 16], Rate=opts.Rate);
    FR = opts.Rate;  W = opts.WindowS;  res = [];  nWin = NaN;  Best = struct([]);
    for hh = opts.Hemis
        H = B.(char(hh));  Sh = H.S;  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;
        av = full(sum(Mm, 2));  C = double(E.(char(hh)).C);  nFr = size(C, 2);  nWin = floor(nFr/(W*FR));
        g = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=H.lbo.L, M=Mm, MaxDepth=max(opts.Depths));
        sh = randi([W*FR, nFr - W*FR], size(C, 1), 1);
        Cs = zeros(size(C));  for k = 1:size(C, 1), Cs(k, :) = circshift(C(k, :), sh(k)); end
        for d = opts.Depths
            G = rheome.geom.tiles(Tr, Sh, d, Ruler=g);  P = double(G.P);  Wd = (P' * (av .* Phi)) ./ (P' * av);
            Y0 = Wd * C;  Ys0 = Wd * Cs;
            for jj = opts.JLevels
                [Yr, rate, nF] = i_winblock(Y0, FR, jj, W*FR, nWin);  Ysr = i_winblock(Ys0, FR, jj, W*FR, nWin);
                b0 = quantile(reshape(Ysr(:, 1:nWin*nF), [], 1), opts.BaseQ);
                ss = zeros(nWin, 1);  sst = cell(nWin, 1);
                for w = 1:nWin
                    pth = rheome.detect.tilepath(Ysr(:, (w-1)*nF + (1:nF)), G, Baseline=b0, SampleRate=rate);
                    if pth.nTracks, ss(w) = pth.tracks(1).score; sst{w} = i_stats(pth, G, rate); end
                end
                cal = 1:2:nWin;  tst = 2:2:nWin;
                thr = i_rate_thr(ss(cal), opts.Alpha * W);
                falseTest = sum(ss(tst) > thr) / (numel(tst) * W);
                rs = zeros(nWin, 1);  rst = cell(nWin, 1);
                for w = 1:nWin
                    pth = rheome.detect.tilepath(Yr(:, (w-1)*nF + (1:nF)), G, Baseline=b0, SampleRate=rate);
                    if pth.nTracks, rs(w) = pth.tracks(1).score; rst{w} = i_stats(pth, G, rate); end
                end
                pr = rs > thr;  ps = ss(tst) > thr;
                [bs, bw] = max(rs);
                if bs > 0
                    pth = rheome.detect.tilepath(Yr(:, (bw-1)*nF + (1:nF)), G, Baseline=b0, SampleRate=rate);
                    Best = [Best, struct('hemi', hh, 'depth', d, 'jj', jj, 'rate', rate, 'window', bw, ...
                        'score', bs, 'passed', bs > thr, 'track', pth.tracks(1), 'G', G)]; %#ok<AGROW>
                end
                R = vertcat(rst{pr});  Sx = vertcat(sst{tst(ps)});
                if isempty(R), R = i_stats([], G, rate); end
                if isempty(Sx), Sx = i_stats([], G, rate); end
                res = [res; double(hh == "R") d median(G.diameterMM) rate b0 thr falseTest ...
                       mean(pr) median(R.durationS) median(R.netMM) median(R.speedMS) median(R.straight, 'omitnan') ...
                       mean(R.netMM > median(G.diameterMM)) height(Sx) median(Sx.netMM) median(Sx.speedMS)]; %#ok<AGROW>
                fprintf('[grouptrack %s] %s d%d %7.3g fps: pass %.2f, false %.2f/s, net %.0f mm, speed %.3f m/s\n', ...
                    name, hh, d, rate, mean(pr), falseTest, median(R.netMM), median(R.speedMS));
            end
        end
        clear C Cs Y0 Ys0
    end
    X = array2table(res, 'VariableNames', {'hemiR','depth','tileMM','rate','baseline','thr','falseTestS', ...
        'fracPass','durS','netMM','speedMS','straight','fracBeyondTile','nSurPass','surNetMM','surSpeedMS'});
    T = rheome.scale.rows("grouptrack", "n_windows", nWin, "windows");
    m = ["tile_mm" "false_rate_test" "frac_pass" "dur_s" "net_mm" "speed_ms" "straight" ...
         "frac_beyond_tile" "sur_n_pass" "sur_net_mm" "sur_speed_ms"];
    u = ["mm" "1/s" "fraction" "s" "mm" "m/s" "ratio" "fraction" "windows" "mm" "m/s"];
    col = ["tileMM" "falseTestS" "fracPass" "durS" "netMM" "speedMS" "straight" "fracBeyondTile" ...
           "nSurPass" "surNetMM" "surSpeedMS"];
    hn = ["L" "R"];
    for i = 1:height(X)
        b = sprintf('%s/d%d/%gfps', hn(X.hemiR(i) + 1), X.depth(i), X.rate(i));
        T = [T; rheome.scale.rows("grouptrack", m, X{i, cellstr(col)}, u, b)]; %#ok<AGROW>
    end
end

% Dyadic block means inside each window of nW frames, windows concatenated: [nTile x nF*nWin].
function [Y, rate, nF] = i_winblock(Y0, FR, jj, nW, nWin)
    bs = 2^jj;  nF = floor(nW/bs);  rate = FR/bs;  nT = size(Y0, 1);
    Z = reshape(Y0(:, 1:nW*nWin), nT, nW, nWin);  Z = Z(:, 1:nF*bs, :);
    Y = reshape(mean(reshape(Z, nT, bs, nF, nWin), 2), nT, nF*nWin);
end

function t = i_rate_thr(sc, allowPerWindow)
    sc = sort(sc(:), 'descend');  k = floor(allowPerWindow * numel(sc));
    if k >= numel(sc), t = 0; else, t = sc(k+1); end
end

function s = i_stats(pth, G, rate)
    if isempty(pth)
        s = array2table(zeros(0, 9), 'VariableNames', {'track','nFrames','durationS','netMM','pathMM', ...
            'straight','hops','meanValue','speedMS'});
        return
    end
    pth.tracks = pth.tracks(1);  s = rheome.detect.trackstats(pth, G, rate);
end

% Author: Diellor Basha, 2026
