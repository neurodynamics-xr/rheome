function [T, X] = measure_inject(name, S, opts)
% SCALE.MEASURE_INJECT  Known movers injected into this participant's own recording, tracked back (MS1 G7, Fig. 8).
%
%   [T, X] = rheome.scale.measure_inject(name)
%   [T, X] = rheome.scale.measure_inject(name, S, Depths=[3 7 8], JLevels=[2 4], KAmp=2, NPer=12)
%
% The group version of alpha_inject_omega.m (nxr-cortical-flow-matlab; sub-0002, left hemisphere, v18
% Fig. 8 and section 5.2), same source, same tracker, same null, per hemisphere:
%   source   a Laplace-Beltrami heat-kernel blob (SigmaMM 15) along the normal, 10 Hz carrier, forwarded
%            through the participant's leadfield; movers are waypoint blobs blended at the sensors by
%            rheome.flow.movingpeak (exact: the forward is linear)
%   cases    still, swap (two blobs SwapMM apart trading amplitude: nothing travels), and movers at
%            Speeds 0.02 / 0.05 / 0.1 m/s along a cortical shortest path of v*W; NPer random 2 s windows
%            each, present for the window +-1 s with 0.25 s cosine ramps
%   amplitude KAmp x the PARTICIPANT'S OWN median per-window alpha envelope peak (over a unit blob's)
%   readout  each injected window band-passed (8-16 Hz) and Hilbert-transformed on an 8 s segment exactly
%            as the data were (rheome.flow.envelopemodes), minimum norm, LB modes, tile means; its paired
%            control is the same segment without the injection
%   tracker  rheome.detect.tilepath, top 3 paths; the baseline and threshold from a mode-shift surrogate of
%            the recorded envelope (rheome.detect.tilepathnull, Alpha false paths/s, odd windows calibrate,
%            even windows test); FOUND = a passing path with >= half its frames within NearMM of the truth
%   real     every 2 s window of the recording through the same tracker and threshold (case "real")
%
% X (inject.csv): one row per hemisphere x case x injection x depth x rate, with the injected and the
% control readouts (pass found near netMM speedMS straight durS, c_*), the window's `half` of the
% recording (injections 1..NPer/2 come from the first half; real windows by position) for the ICC, and
% the surrogate's held-out false rate `nullFalseS`.
% T (rheome.scale.rows, analysis "inject"), band "<case>/d<depth>/<fps>fps", medians over both hemispheres:
%   hit_rate      fraction of injections FOUND                                   fraction
%   speed_ratio   median recovered net speed / true speed, found movers ⭐       ratio
%   net_mm        median net displacement of the found path                      mm
%   control_found fraction of paired controls with a passing path near the site  fraction
%   false_rate    surrogate passes per second on held-out windows (band d/fps)   1/s
%   real_frac_pass, real_net_mm   real alpha: windows passing the null, their median net displacement
% plus hit_rate and speed_ratio per half (band suffix "_h1" / "_h2").
%
% ⚠ Depth 3 tiles are lobar (8 per hemisphere): there a 0.1 m/s mover crosses at most one tile boundary,
%   so its speed ratio is quantised. It is in the design as the independence limit (MS1 section 9.5).
% ⚠ A window holds two objects -- the injection and the recording's own alpha -- hence the top 3 paths.
%
% See also: rheome.scale.measure_grouptrack, rheome.scale.measure_trackfactorial, rheome.detect.tilepath,
%           rheome.detect.tilepathnull, rheome.flow.movingpeak, rheome.flow.envelopemodes
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depths double = [3 7 8]
        opts.JLevels double = [2 4]
        opts.KAmp (1,1) double = 2
        opts.NPer (1,1) double {mustBeInteger, mustBePositive} = 12
        opts.Cases string = ["still" "swap" "move02" "move05" "move10"]
        opts.Speeds double = [0 0 0.02 0.05 0.1]
        opts.SigmaMM (1,1) double = 15
        opts.SwapMM (1,1) double = 40
        opts.NearMM (1,1) double = 50
        opts.Alpha (1,1) double = 0.1
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 23
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    rng(opts.Seed);
    FR = 300;  W = 2;  MARGIN = 1;  PAD = 2;  F0 = 10;
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    E = rheome.flow.envelopemodes(name, Band=[8 16], Rate=FR);
    dec = round(fs/FR);  nSeg = (W + 2*MARGIN + 2*PAD)*fs;  tSeg = (0:nSeg-1)/fs - (MARGIN + PAD);
    env = i_env(tSeg, W, MARGIN, 0.25);  cIdx = (MARGIN + PAD)*fs + (1:dec:W*fs);
    uW = min(max(tSeg/W, 0), 1);
    vn = {'pass','found','near','netMM','speedMS','straight','durS'};
    X = table();
    for hh = opts.Hemis
        H = S.B.(char(hh));  Sh = H.S;  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  lam = H.lbo.Lambda(:);
        av = full(sum(Mm, 2));  gv = double(H.gv(:));  nV = numel(gv);  nrm = double(Sh.VertNormals);
        rows = reshape((gv' - 1)*3 + (1:3)', [], 1);  K = S.Res.ImagingKernel(rows, :);  GL = S.G(:, rows);
        hk = exp(-lam*(opts.SigmaMM*1e-3)^2/2);
        C = double(E.(char(hh)).C);  nFr = size(C, 2);  nWin = floor(nFr/(W*FR));
        ge = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=H.lbo.L, M=Mm, MaxDepth=max(opts.Depths));
        sh = randi([W*FR, nFr - W*FR], size(C, 1), 1);
        Cs = zeros(size(C));  for k = 1:size(C, 1), Cs(k, :) = circshift(C(k, :), sh(k)); end
        nd = numel(opts.Depths);  nj = numel(opts.JLevels);
        Gs = cell(1, nd);  Wd = Gs;  Dcv = Gs;  NUL = cell(nd, nj);
        for di = 1:nd
            Gs{di} = rheome.geom.tiles(Tr, Sh, opts.Depths(di), Ruler=ge);  P = double(Gs{di}.P);
            Wd{di} = (P'*(av.*Phi))./(P'*av);  Dcv{di} = distances(ge, Gs{di}.centre);
            for ji = 1:nj
                [Ys, rate, nF] = i_winblock(Wd{di}*Cs, FR, opts.JLevels(ji), W*FR, nWin);
                NUL{di, ji} = rheome.detect.tilepathnull(Ys, Gs{di}, nF, rate, Alpha=opts.Alpha, WindowS=W);
                % real alpha: every window against the same threshold (best path, as grouptrack)
                Yr = i_winblock(Wd{di}*C, FR, opts.JLevels(ji), W*FR, nWin);
                for w = 1:nWin
                    sr = i_path(Yr(:, (w-1)*nF + (1:nF)), Gs{di}, NUL{di, ji}, rate, [], [], opts.NearMM);
                    X = [X; i_row(hh, "real", NaN, w, w, 1 + (w > nWin/2), opts.Depths(di), Gs{di}, rate, ...
                                  NUL{di, ji}, sr, nan(1, 7), vn)]; %#ok<AGROW>
                end
            end
        end
        clear Cs
        % amplitude reference: the real envelope's typical per-window peak, and a unit blob's
        ref = zeros(nWin, 1);
        for w = 1:nWin, ref(w) = max(Phi*mean(C(:, (w-1)*W*FR + (1:W*FR)), 2)); end
        s1 = zeros(20, 1);  tt = (0:W*fs-1)/fs;
        for q = 1:20
            g0 = GL*i_blob(randi(nV), Phi, Mm, hk, nrm);
            A = rheome.flow.activation(K*hilbert(filtfilt(bp, (g0*cos(2*pi*F0*tt)).')).');
            s1(q) = max(mean(A, 2));
        end
        amp = opts.KAmp*median(ref)/median(s1);
        for ci = 1:numel(opts.Cases)
            for n = 1:opts.NPer
                hf = 1 + (n > opts.NPer/2);  lims = [3, floor(nWin/2); floor(nWin/2) + 1, nWin - 2];
                w = randi(lims(hf, :));  s0 = (w-1)*W*fs - (MARGIN + PAD)*fs + (1:nSeg);
                v0 = randi(nV);  car = cos(2*pi*F0*tSeg + 2*pi*rand);  v = opts.Speeds(ci);
                switch opts.Cases(ci)
                    case "still"
                        bsig = (GL*i_blob(v0, Phi, Mm, hk, nrm))*(env.*car);  vT = repmat(v0, 1, numel(cIdx));
                    case "swap"
                        d = distances(ge, v0);  [~, v1] = min(abs(d - opts.SwapMM*1e-3));
                        car2 = cos(2*pi*F0*tSeg + 2*pi*rand);
                        bsig = (GL*i_blob(v0, Phi, Mm, hk, nrm))*(env.*cos(pi/2*uW).*car) + ...
                               (GL*i_blob(v1, Phi, Mm, hk, nrm))*(env.*sin(pi/2*uW).*car2);
                        vT = repmat(v0, 1, numel(cIdx));  vT(uW(cIdx) > 0.5) = v1;
                    otherwise
                        Lp = v*W;  cand = [];  tries = 0;
                        while isempty(cand)
                            tries = tries + 1;
                            if tries > 200, error('scale:inject:path', 'No %.0f mm shortest path on the %s hemisphere.', 1e3*Lp, hh); end
                            v0 = randi(nV);  d = distances(ge, v0);  cand = find(abs(d - Lp) < 2e-3);
                        end
                        [pv, ~, eg] = shortestpath(ge, v0, cand(randi(numel(cand))));
                        cu = [0; cumsum(ge.Edges.Weight(eg))];  keep = [true; diff(cu) > 0];  pv = pv(keep);  cu = cu(keep);
                        Ga = zeros(size(GL, 1), numel(pv));
                        for k = 1:numel(pv), Ga(:, k) = GL*i_blob(pv(k), Phi, Mm, hk, nrm); end
                        PK = rheome.flow.movingpeak(Ga, cu, SpeedMS=v, SampleRate=fs);
                        nTr = size(PK.X, 2);  inW = find(tSeg >= 0, 1);
                        idx = min(max((1:nSeg) - inW + 1, 1), nTr);                 % hold at the ends
                        bsig = PK.X(:, idx).*(env.*car);  vT = pv(PK.k(idx(cIdx)))';
                end
                Cin = i_envmodes(F(:, s0) + amp*bsig, bp, cIdx, K, Phi, Mm);
                Cct = i_envmodes(F(:, s0), bp, cIdx, K, Phi, Mm);
                for di = 1:nd
                    for ji = 1:nj
                        [Yi, rate, nF] = i_winblock(Wd{di}*Cin, FR, opts.JLevels(ji), W*FR, 1);
                        Yc = i_winblock(Wd{di}*Cct, FR, opts.JLevels(ji), W*FR, 1);
                        vTb = vT((0:nF-1)*2^opts.JLevels(ji) + ceil(2^opts.JLevels(ji)/2));
                        si = i_path(Yi, Gs{di}, NUL{di, ji}, rate, Dcv{di}, vTb, opts.NearMM);
                        sc = i_path(Yc, Gs{di}, NUL{di, ji}, rate, Dcv{di}, vTb, opts.NearMM);
                        X = [X; i_row(hh, opts.Cases(ci), v, n, w, hf, opts.Depths(di), Gs{di}, rate, ...
                                      NUL{di, ji}, si, sc, vn)]; %#ok<AGROW>
                    end
                end
            end
            fprintf('[inject %s] %s %-7s: %d injections\n', name, hh, opts.Cases(ci), opts.NPer);
        end
        clear C
    end
    T = i_metrics(X, opts);
end

function r = i_row(hh, cs, v, n, w, hf, d, G, rate, Nl, si, sc, vn)
    r = [table(hh, cs, v, n, w, hf, d, median(G.diameterMM), rate, Nl.falseTestS, ...
               'VariableNames', {'hemi','case','speed','n','win','half','depth','tileMM','rate','nullFalseS'}), ...
         array2table([si sc], 'VariableNames', [vn strcat('c_', vn)])];
end

function T = i_metrics(X, opts)
    T = table();  r = @(m, v, u, b) rheome.scale.rows("inject", m, v, u, b);
    for d = opts.Depths
        for rt = unique(X.rate(X.depth == d))'
            k = X.depth == d & abs(X.rate - rt) < 1e-9;  b = sprintf('d%d/%gfps', d, rt);
            q = X(k & X.case == "real", :);
            T = [T; r(["false_rate" "real_frac_pass" "real_net_mm"], ...
                      [median(X.nullFalseS(k)) mean(q.pass) median(q.netMM(q.pass == 1), 'omitnan')], ...
                      ["1/s" "fraction" "mm"], b)]; %#ok<AGROW>
            for c = opts.Cases
                q = X(k & X.case == c, :);  f = q(q.found == 1, :);  bc = c + "/" + b;
                sr = f.speedMS ./ f.speed;  sr(f.speed == 0) = NaN;
                T = [T; r(["hit_rate" "speed_ratio" "net_mm" "control_found"], ...
                          [mean(q.found) median(sr, 'omitnan') median(f.netMM, 'omitnan') mean(q.c_found)], ...
                          ["fraction" "ratio" "mm" "fraction"], bc)]; %#ok<AGROW>
                for hf = 1:2
                    qh = q(q.half == hf, :);  fh = qh(qh.found == 1, :);  s = fh.speedMS ./ fh.speed;  s(fh.speed == 0) = NaN;
                    T = [T; r(["hit_rate" "speed_ratio"], [mean(qh.found) median(s, 'omitnan')], ["fraction" "ratio"], bc + "_h" + hf)]; %#ok<AGROW>
                end
            end
        end
    end
end

% Heat-kernel blob at vertex v along the normal, interleaved [x1;y1;z1;x2;...] to match the gain rows.
function j = i_blob(v, Phi, Mm, hk, nrm)
    d = zeros(size(Phi, 1), 1);  d(v) = 1;
    a = Phi*(hk.*(Phi'*(Mm*d)));  a = a/max(a);
    j = reshape((a.*nrm)', [], 1);
end

% 1 inside [-margin, W+margin] with cosine ramps of r seconds at both outer ends.
function e = i_env(t, W, margin, r)
    e = ones(size(t));  a = -margin;  b = W + margin;
    e(t < a | t > b) = 0;
    up = t >= a & t < a + r;  e(up) = 0.5 - 0.5*cos(pi*(t(up) - a)/r);
    dn = t > b - r & t <= b;  e(dn) = 0.5 - 0.5*cos(pi*(b - t(dn))/r);
end

% The envelope of a sensor segment as LB mode coefficients over the central window.
function Cw = i_envmodes(Fs, bp, cIdx, K, Phi, Mm)
    Fa = hilbert(filtfilt(bp, Fs.')).';
    Cw = Phi'*(Mm*rheome.flow.activation(K*Fa(:, cIdx)));
end

function [Y, rate, nF] = i_winblock(Y0, FR, jj, nW, nWin)
    bs = 2^jj;  nF = floor(nW/bs);  rate = FR/bs;  nT = size(Y0, 1);
    Z = reshape(Y0(:, 1:nW*nWin), nT, nW, nWin);  Z = Z(:, 1:nF*bs, :);
    Y = reshape(mean(reshape(Z, nT, bs, nF, nWin), 2), nT, nF*nWin);
end

% Top 3 disjoint paths; the injection is FOUND if a passing one stays near it (>= half its frames within
% nearMM). Returns [pass found near netMM speedMS straight durS] of the found path, else the best one.
% Without a truth (Dcv empty: real windows) only the best path is scored.
function r = i_path(Y, G, Nl, rate, Dcv, vTb, nearMM)
    r = [0 0 NaN NaN NaN NaN NaN];
    p = rheome.detect.tilepath(Y, G, Baseline=Nl.baseline, SampleRate=rate, NumPaths=3 - 2*isempty(Dcv));
    if p.nTracks == 0, return; end
    ok = [p.tracks.score] > Nl.threshold;  q = 1;  fnd = false;
    if ~isempty(Dcv)
        nr = zeros(p.nTracks, 1);
        for i = 1:p.nTracks
            x = p.tracks(i);  vv = vTb(x.frames(:));
            nr(i) = mean(Dcv(sub2ind(size(Dcv), x.tiles(:), vv(:))) <= nearMM*1e-3);
        end
        fnd = ok(:) & nr >= 0.5;
        if any(fnd), [~, q] = max(nr.*fnd); end
    else
        nr = NaN;
    end
    one = p;  one.tracks = p.tracks(q);  one.nTracks = 1;  st = rheome.detect.trackstats(one, G, rate);
    r = [any(ok) any(fnd) nr(q) st.netMM st.speedMS st.straight st.durationS];
end

% Author: Diellor Basha, 2026
