function [T, X, G] = measure_periodicflow(name, S, opts)
% SCALE.MEASURE_PERIODICFLOW  Apparent alpha flow on the TOTAL and the PERIODIC envelope, over many tiles.
%
%   [T, X] = rheome.scale.measure_periodicflow(name)
%   [T, X] = rheome.scale.measure_periodicflow(name, rheome.scale.sensors(name), NumTiles=12)
%   [T, X] = rheome.scale.measure_periodicflow(name, S, Centres=c)     % explicit tile centres (samples)
%   [T, X, G] = rheome.scale.measure_periodicflow(...)                  % + the flow in the group gauge
%
% The single-subject apparent-flow and periodic-flow analyses with the figures removed, run on
% NumTiles 2 s tiles spread EVENLY across the recording instead of the single strongest one. The
% report periodicflow_article ends by asking for exactly this: speed against amplitude ACROSS
% tiles, where the samples are independent, because 599 frames of one tile are ~25 samples.
% Per tile, identical to the scripts:
%   Welch PSD of the whole recording (8 s Hann, 50%), aperiodic fit with knee on [1 45] Hz minus
%   60/120 Hz -> rheome.spectral.decompose on the tile's FFT (+-2 s margin) with the unit conversion
%   derived from the data (the ratio of MEANS of |C|^2 and the PSD over the fit band), periodic
%   gains applied with Hermitian symmetry -> for TOTAL and PERIODIC: 8-16 Hz IIR band-pass,
%   Hilbert AT THE SENSORS, decimate 600 -> 300 Hz, left-hemisphere minimum-norm current
%   (S.Res, SnrFixed 3, amplitude), |J^| (rheome.flow.activation), Horn-Schunck Alpha=1 (rheome.flow.apparent).
%
% Returns rheome.scale.rows (analysis "periodicflow"); band = "total" | "periodic" | "" (spectrum):
%   periodic_fraction, _q25, _q75   PSD periodic share of 8-16 Hz power, over channels  fraction
%   aperiodic_share                  1 - periodic_fraction: the 1/f share of the band     fraction
%   aperiodic_exponent, fit_r2       chi of the knee fit, and its log-space R^2
%   tile_inband_periodic_median      the split's own in-band share, median over tiles     fraction
%   tile_bg_above_obs_median         % bins where the background exceeds |C|^2 (0/100 = unit error)
%   n_tiles, frames_per_tile         n_tiles counts ACCEPTED tiles only
%   n_tiles_rejected                 tiles dropped by the stationarity screen (below)
%   clean_start_s, clean_end_s       the span tiles were placed in (rheome.scale.cleanspan)
% per band (total, periodic):
%   speed_median, speed_p95          median over tiles of each tile's median / p95 |v|     m/s
%   curl_p95, div_p95                median over tiles of p95 |curl v|, |div v|           1/s
%   r_within_median                  median over tiles of corr(frame mean envelope, frame median speed)
%   neff_median                      median Bartlett effective n of that correlation      samples
%   r_within_pooled, z_within        n_eff-weighted Fisher pool of the within-tile r, and its Stouffer z
%   r_across, p_across               corr ACROSS tiles of tile mean envelope vs tile median speed ⭐
%   env_spatial_r_median             (periodic only) spatial corr of the two tile-mean envelopes
% X is the per-tile table (one row per tile x band), written as flowtiles.csv by rheome.scale.run; its
% power_ratio and rejected columns say which tiles the summaries used.
%
% G is the velocity IN THE GROUP GAUGE, written as gaugeflow.csv: one row per chart x band x patch
% (rheome.geom.spherepatches, ico-PatchLevel on the registered sphere, so patch k is the same place in
% every subject), pooled over every frame of every accepted tile and every vertex of the patch:
%   n_vertices, n_samples                vertices in the patch, vertex-frames pooled
%   vn_mean, vw_mean                     mean north / west velocity (rheome.geom.sphereframe)    m/s
%   tnn, tnw, tww                        mean v (x) v in (north, west): the orientation tensor, (m/s)^2
%   env_mean                             mean activation |J^| (rheome.flow.activation units)
%   colat_deg, lon_deg, spread_deg, has_pole, excluded   the patch, and whether the gauge holds there
% chart "z" has its poles at the sphere's +-z (the group gauge); chart "x" at +-x, to read the polar
% patches z excludes. Left hemisphere only, as the flow. Empty if the surface has no Reg.Sphere.
%
% ⚠⚠ TILES ARE SCREENED FOR STATIONARITY, AND A BAD TILE NO LONGER COSTS THE OTHER ELEVEN. The split
% compares each tile's |C|^2 with a background fitted on the WHOLE recording, so it is only valid
% where the tile looks like the recording. Three guards, in order:
%   1. rheome.scale.cleanspan trims artefact blocks off both edges and the Welch PSD and the tile centres
%      use only [clean_start, clean_end]. Before 2026-10-07 tile 12 ended on the last sample, and on
%      ~3% of a resting cohort that sample sat inside a ~3.5 s common-mode end-of-recording artefact.
%      ⚠ The Welch PSD also leaves out every segment touching an INTERIOR bad block
%      (rheome.scale.cleanwelch). Without that, a 50x artefact at 20-25 s of one subject made the PSD 20x
%      a typical stretch and all 12 clean tiles read power_ratio 0.03-0.07 (another subject: 0.15-0.33),
%      so both failed with notiles. With no interior bad block the PSD is the pwelch call, bit for bit.
%   2. a tile whose window (margins included) touches a remaining bad block is rejected.
%   3. power_ratio = (data-derived |C|^2/PSD scale) / (fs*nT/2, the analytic one). It is 0.88-1.13
%      on one subject's clean tiles, and over EVERY 6 s window it spans 0.44-1.82 on a reference subject
%      (595 windows) and 0.87-1.67 on the first subject's clean span. Outside [1/MaxPowerRatio, MaxPowerRatio]
%      (default 3) the tile is rejected. Its old tile 12 read ~3e4 -- what tripped
%      spectral:decompose:scale as a "UNIT MISMATCH" -- and its old tile 1 read 3.1.
% Rejected tiles stay in X as NaN rows with the reason.
%
% ⚠ BARTLETT n_eff = n / (1 + 2 sum_{k=1}^{n/5} (1-k/n) rho_x(k) rho_y(k)). At 300 Hz both series
% have lag-1 autocorrelation ~0.99, so the nominal n = 599 overstates the evidence ~20-fold.
% ⚠ a reference subject's strongest tile read r = -0.425 (total) and -0.162 (periodic), periodic fraction
% 0.693, speed 15.1 / 13.2 mm/s; Centres= reproduces that tile.
%
% See also: rheome.flow.apparent, rheome.flow.activation, rheome.spectral.decompose, rheome.spectral.aperiodic, rheome.scale.run,
%           rheome.scale.cleanspan, rheome.scale.cleanwelch
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.NumTiles (1,1) double {mustBeInteger, mustBePositive} = 12
        opts.Centres double = []
        opts.Band (1,2) double = [8 16]
        opts.Rate (1,1) double = 300
        opts.TileS (1,1) double = 2
        opts.MarginS (1,1) double = 2
        opts.Alpha (1,1) double = 1
        opts.SegS (1,1) double = 8
        opts.FitRange (1,2) double = [1 45]
        opts.MaxPowerRatio (1,1) double {mustBeGreaterThan(opts.MaxPowerRatio, 1)} = 3
        opts.EdgeS (1,1) double = 10
        opts.PatchLevel (1,1) double = 3
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    band = opts.Band;  FITR = opts.FitRange;
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    dec = round(fs / opts.Rate);  rate = fs / dec;
    SL = S.B.L.S;  gvL = double(S.B.L.gv(:))';
    rowsL = reshape((gvL - 1) * 3 + (1:3)', [], 1);
    K = S.Res.ImagingKernel(rowsL, :);

    % 0. where the recording is clean enough to tile (edges trimmed, artefact blocks marked)
    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);
    cs = CS.first:CS.last;

    % 1. the PSD split of the clean span
    nf = 2^nextpow2(opts.SegS * fs);
    if numel(cs) < nf, error('scale:periodicflow:short', 'clean span of %.1f s is shorter than one %g s Welch segment', numel(cs)/fs, nf/fs); end
    [Pfit, ffit, W] = rheome.scale.cleanwelch(F, fs, cs, CS.mask, nf);
    if W.nUsed < W.nSeg
        fprintf('[periodicflow %s] Welch PSD on %d of %d segments (%d touch an artefact block)\n', name, W.nUsed, W.nSeg, W.nSeg - W.nUsed);
    end
    mf = ffit >= FITR(1) & ffit <= FITR(2);  for h = [60 120], mf = mf & ~(ffit > h-2 & ffit < h+2); end
    ap = rheome.spectral.aperiodic(Pfit(mf, :), ffit(mf), struct('knee', true));
    fm = ffit(mf);  mb = fm >= band(1) & fm < band(2);  Pk = Pfit(mf, :);
    pf = sum(max(Pk(mb, :) - ap.Pap(mb, :), 0), 1) ./ sum(Pk(mb, :), 1);

    % 2. tiles, each with a margin so the frequency-domain gains do not ring into it
    Lm = round(opts.MarginS * fs);  Lt = round(opts.TileS * fs);  N = size(F, 2);
    lo = CS.first - 1 + Lm + floor(Lt/2) + 1;  hi = CS.last - Lm - ceil(Lt/2);
    if hi < lo, error('scale:periodicflow:short', 'clean span of %.1f s (of %.1f s) holds no %g s tile with %g s margins', numel(cs)/fs, N/fs, opts.TileS, opts.MarginS); end
    c = opts.Centres(:)';
    if isempty(c), c = round(linspace(lo, hi, min(opts.NumTiles, hi - lo + 1))); end
    if any(c < lo | c > hi), error('scale:periodicflow:centre', 'tile centres must lie in [%d %d]', lo, hi); end
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', band(1), ...
                    'HalfPowerFrequency2', band(2), 'SampleRate', fs);
    nm = ["total" "periodic"];
    X = table();
    GZ = i_gaugeprep(S, gvL, SL, opts.PatchLevel);      % [] without a registration sphere
    if ~isempty(GZ), acc = zeros(SL.nV, 7, 2, 2); end   % vertex x [n vn vw nn nw ww env] x band x chart
    for ti = 1:numel(c)
        seg = (c(ti) - floor(Lt/2) - Lm) : (c(ti) + ceil(Lt/2) - 1 + Lm);
        Xw = F(:, seg);  nT = numel(seg);
        Fw = fft(Xw, [], 2);  fw = (0:nT-1) * fs / nT;  pos = 2:floor(nT/2) + 1;
        Pint = interp1(ffit, Pfit, fw(pos)', 'linear', 'extrap')';
        inFit = fw(pos)' >= FITR(1) & fw(pos)' <= FITR(2);
        sc = mean(abs(Fw(:, pos(inFit))).^2, 'all') / mean(Pint(:, inFit), 'all');
        pr = sc / (fs * nT / 2);           % 1 when the tile has the recording's fit-band power
        why = "";
        if any(CS.mask(seg)), why = "artefact block in window";
        elseif pr > opts.MaxPowerRatio || pr < 1/opts.MaxPowerRatio
            why = string(sprintf("power ratio %.3g outside [1/%g %g]", pr, opts.MaxPowerRatio, opts.MaxPowerRatio));
        end
        if why == ""
            try
                D = rheome.spectral.decompose(Fw(:, pos), fw(pos)', struct('fitPower', (Pfit(mf, :) * sc)', ...
                                        'fitF', ffit(mf), 'knee', true));
            catch e
                if e.identifier ~= "spectral:decompose:scale", rethrow(e); end
                why = "decompose: " + string(e.message);
            end
        end
        if why ~= ""
            X = [X; i_rejected(ti, c(ti)/fs, nm, pr, why)]; %#ok<AGROW>
            fprintf('[periodicflow %s] tile %2d/%d at %6.1f s: REJECTED, %s\n', name, ti, numel(c), c(ti)/fs, why);
            continue
        end
        bgAbove = 100 * mean(rheome.spectral.evaluate(D.ap, fw(pos)')' > abs(Fw(:, pos)).^2, 'all');
        Hp = zeros(size(Fw));  Hp(:, pos) = D.hper;
        for j = pos, m = nT - j + 2; if m <= nT && m > j, Hp(:, m) = Hp(:, j); end, end
        Xper = real(ifft(Fw .* Hp, [], 2));
        ba = filtfilt(bp, Xw.');  bp2 = filtfilt(bp, Xper.');
        shareBand = sum(bp2(:).^2) / sum(ba(:).^2);
        keep = (Lm + 1):(nT - Lm);  idx = keep(1:dec:end);
        Am = cell(1, 2);
        for v = 1:2
            if v == 1, Xv = Xw; else, Xv = Xper; end
            Fa = hilbert(filtfilt(bp, Xv.')).';                 % ⚠ time down columns, then back
            A = rheome.flow.activation(K * Fa(:, idx), Rate=rate, Band=band);
            o = rheome.flow.apparent(A, SL, Rate=rate, Band=band, Alpha=opts.Alpha);
            mE = mean(A(:, 1:end-1), 1);  mS = median(o.speed, 1);
            r = corr(mE(:), mS(:));  ne = i_neff(mE(:), mS(:));
            Am{v} = mean(A, 2);
            if ~isempty(GZ)
                for ch = 1:2, acc(:,:,v,ch) = acc(:,:,v,ch) + i_frameSums(o.velocity * rate, A(:, 1:end-1), GZ.fr{ch}); end
            end
            X = [X; table(ti, c(ti)/fs, nm(v), o.nT, mean(mE), median(o.speed(:)), prctile(o.speed(:), 95), ...
                 prctile(abs(o.vorticity(:)), 95), prctile(abs(o.divergence(:)), 95), r, ne, ...
                 shareBand, bgAbove, D.partition, NaN, o.seconds, pr, "", ...
                 'VariableNames', i_cols())]; %#ok<AGROW>
        end
        X.env_spatial_r(end) = corr(Am{1}, Am{2});
        fprintf('[periodicflow %s] tile %2d/%d at %6.1f s: r total %+.3f, periodic %+.3f, share %.2f\n', ...
            name, ti, numel(c), c(ti)/fs, X.r_within(end-1), X.r_within(end), shareBand);
    end

    G = table();
    if ~isempty(GZ), G = i_gaugeTable(acc, GZ, nm); end
    ok = X.rejected == "";
    if ~any(ok), error('scale:periodicflow:notiles', 'all %d tiles were rejected by the stationarity screen', numel(c)); end
    nAcc = sum(ok & X.band == "total");
    T = rheome.scale.rows("periodicflow", ...
        ["periodic_fraction" "periodic_fraction_q25" "periodic_fraction_q75" "aperiodic_share" ...
         "aperiodic_exponent" "fit_r2" "tile_inband_periodic_median" "tile_bg_above_obs_median" ...
         "n_tiles" "frames_per_tile" "n_tiles_rejected" "clean_start_s" "clean_end_s"], ...
        [median(pf) prctile(pf, 25) prctile(pf, 75) 1 - median(pf) median(ap.exponent) median(ap.r2) ...
         median(X.inband_periodic(ok & X.band == "total")) median(X.bg_above_obs_pct(ok & X.band == "total")) ...
         nAcc median(X.frames(ok)) numel(c) - nAcc (CS.first - 1)/fs CS.last/fs], ...
        ["fraction" "fraction" "fraction" "fraction" "chi" "R2" "fraction" "percent" "tiles" "frames" ...
         "tiles" "s" "s"]);
    for b = nm
        x = X(ok & X.band == b, :);
        w = max(x.neff - 3, 0);  z = atanh(max(min(x.r_within, 0.999999), -0.999999));
        rp = tanh(sum(w .* z) / max(sum(w), eps));  zs = sum(sqrt(w) .* z) / sqrt(height(x));
        if height(x) >= 3, [ra, pa] = corr(x.env_mean, x.speed_median); else, ra = NaN; pa = NaN; end
        m = ["speed_median" "speed_p95" "curl_p95" "div_p95" "r_within_median" "neff_median" ...
             "r_within_pooled" "z_within" "r_across" "p_across"];
        v = [median(x.speed_median) median(x.speed_p95) median(x.curl_p95) median(x.div_p95) ...
             median(x.r_within) median(x.neff) rp zs ra pa];
        u = ["m/s" "m/s" "1/s" "1/s" "r" "samples" "r" "z" "r" "p"];
        if b == "periodic"
            m(end+1) = "env_spatial_r_median";  v(end+1) = median(x.env_spatial_r);  u(end+1) = "r"; %#ok<AGROW>
        end
        T = [T; rheome.scale.rows("periodicflow", m, v, u, b)]; %#ok<AGROW>
    end
end

function GZ = i_gaugeprep(S, gvL, SL, level)
% the two charts' frames and patches on the left hemisphere (rheome.geom.sphereframe, rheome.geom.spherepatches)
    GZ = [];
    if ~isfield(S, 'Sf') || ~isfield(S.Sf, 'Sphere') || isempty(S.Sf.Sphere)
        fprintf('[periodicflow] no Reg.Sphere on the surface: gaugeflow skipped\n');  return
    end
    Sp = S.Sf.Sphere(gvL, :);  ax = [0 0 1; 1 0 0];
    for ch = 1:2
        GZ.fr{ch} = rheome.geom.sphereframe(SL.Vertices, SL.Faces, SL.VertNormals, Sp, Axis=ax(ch,:));
        GZ.P{ch}  = rheome.geom.spherepatches(Sp, Level=level, Axis=ax(ch,:));
    end
end

function a = i_frameSums(vel, A, fr)
% per vertex, over frames: [n, sum vn, sum vw, sum vn^2, sum vn*vw, sum vw^2, sum env]; the
% frame's singular vertices contribute nothing
    vx = vel(1:3:end, :);  vy = vel(2:3:end, :);  vz = vel(3:3:end, :);
    e1 = fr.e1;  e2 = fr.e2;  ok = ~fr.singular;  e1(~ok,:) = 0;  e2(~ok,:) = 0;
    vn = vx.*e1(:,1) + vy.*e1(:,2) + vz.*e1(:,3);
    vw = vx.*e2(:,1) + vy.*e2(:,2) + vz.*e2(:,3);
    nT = size(vel, 2);
    a = [ok*nT, sum(vn,2), sum(vw,2), sum(vn.^2,2), sum(vn.*vw,2), sum(vw.^2,2), ok.*sum(A,2)];
end

function G = i_gaugeTable(acc, GZ, nm)
% vertex sums -> patch means, one row per chart x band x patch
    G = table();  chart = ["z" "x"];
    for ch = 1:2
        P = GZ.P{ch};  nP = numel(P.n);
        for b = 1:2
            s = zeros(nP, 7);
            for j = 1:7, s(:,j) = accumarray(P.patch, acc(:,j,b,ch), [nP 1]); end
            m = s(:,2:7) ./ s(:,1);
            G = [G; table(repmat(chart(ch),nP,1), repmat(nm(b),nP,1), (1:nP)', P.n, s(:,1), ...
                 m(:,1), m(:,2), m(:,3), m(:,4), m(:,5), m(:,6), P.colat, P.lon, P.spread, ...
                 P.hasPole, P.excluded, 'VariableNames', {'chart','band','patch','n_vertices', ...
                 'n_samples','vn_mean','vw_mean','tnn','tnw','tww','env_mean','colat_deg', ...
                 'lon_deg','spread_deg','has_pole','excluded'})]; %#ok<AGROW>
        end
    end
end

function v = i_cols()
    v = {'tile','centre_s','band','frames','env_mean','speed_median','speed_p95', ...
         'curl_p95','div_p95','r_within','neff','inband_periodic','bg_above_obs_pct', ...
         'partition_err','env_spatial_r','flow_seconds','power_ratio','rejected'};
end

function R = i_rejected(ti, cs, nm, pr, why)
% the two NaN rows (total, periodic) a rejected tile leaves in the per-tile table
    n = numel(nm);  z = NaN(n, 1);
    R = table(repmat(ti, n, 1), repmat(cs, n, 1), nm(:), z, z, z, z, z, z, z, z, z, z, z, z, z, ...
              repmat(pr, n, 1), repmat(string(why), n, 1), 'VariableNames', i_cols());
end

function n = i_neff(x, y)
% Bartlett's effective sample size of corr(x, y) for two autocorrelated series, lags up to n/5
    n0 = numel(x);  M = max(1, floor(n0/5));
    rx = i_acf(x, M);  ry = i_acf(y, M);  k = (1:M)';
    n = n0 / (1 + 2 * sum((1 - k/n0) .* rx .* ry));
    n = min(max(n, 2), n0);
end

function r = i_acf(x, M)
    x = x - mean(x);  d = sum(x.^2);  r = zeros(M, 1);
    for k = 1:M, r(k) = sum(x(1:end-k) .* x(1+k:end)) / d; end
end

% Author: Diellor Basha, 2026
