function [T, X] = measure_trackfactorial(name, S, opts)
% SCALE.MEASURE_TRACKFACTORIAL  The 0.52 diagnosis as a 2^5 factorial on this participant (MS1 G8, Fig. 9).
%
%   [T, X] = rheome.scale.measure_trackfactorial(name)                       % head level "own" only
%   [T, X] = rheome.scale.measure_trackfactorial(name, S, RefHead="subject01") % both head levels
%
% MS1 v18 section 5.3: exact sphere movers tracked through subject01's leadfield read 0.52 of their speed;
% injections into sub-0002's recording (Fig. 8, rheome.scale.measure_inject) read 0.93-1.0. The two tests
% differ in five ways. Here each is a factor at two levels, crossed, in the same participant:
%   ruler       "fold"   net displacement on the folded cortex (rheome.detect.trackstats, edge-graph geodesic)
%                        over the true folded net speed (cortical geodesic from the true start to end vertex)
%               "proj"   the path's tile centres projected onto the EXACT trajectory on the registration
%                        sphere, slope of arc-length position against time, over v      (readout factor)
%   signal      "carrier" 10 Hz carrier on the source; sensors band-passed 8-16 Hz (IIR order 8, filtfilt on
%                        a 6 s segment) and Hilbert-transformed, as the data are
%               "smooth"  the source amplitude itself, no carrier and no filter
%   background  "real"   the participant's own resting record, a random 2 s window (its analytic 8-16 Hz signal)
%               "synth"  LB-smooth (15 mm heat envelope) complex AR(1) source noise with alpha's correlation
%                        time (1/9.71 s), forwarded; carrier level: on a 10 Hz carrier. Scaled so its median
%                        per-window envelope peak equals the real background's
%   matching    "top3"   a passing path among the top 3 that stays within NearMM of the truth (Fig. 8's rule)
%               "best"   the single best path only (Fig. 9's rule)                    (readout factor)
%   head        "own"    the participant's cortex, registration sphere, gain and minimum norm
%               "ref"    RefHead's (subject01): the mover lives on its sphere/cortex, goes through its gain and
%                        inverse on the channels it shares with the participant (by name); the participant's
%                        real background enters through those channels
% Baseline: the DIRECT arm (no instrument: the source map itself + the synthetic background, own head),
% whose speed read 0.97-0.98 on subject01's sphere.
%
% MOVER (the sphere arm's): a geodesic Gaussian (SigmaMM 25, peak 1) on the participant's FreeSurfer
% registration sphere (Reg.Sphere, R ~ 100 mm), along a random great circle at Speeds 0.05-1 m/s for
% min(PathMM, v*W), at a random offset in the W = 2 s window and absent outside it; sphere vertex i IS
% hemisphere cortex vertex i, so it drives a normal current there. Amplitude: KAmp x the real background's
% median per-window envelope peak, over a unit mover's own peak at the same signal level.
% TRACKER: rheome.detect.tilepath on the CORTEX tiles (rheome.geom.tiles, depths Depths), maps pooled
% through the hemisphere's LB modes as in Fig. 8, at the frame rate chosen a priori per (speed, depth)
% so the mover takes about half a tile per frame (Fig. 9's rule); baseline and threshold from NNull
% mover-free windows of the same cell (rheome.detect.tilepathnull, Alpha 0.1/s, odd calibrate / even
% test). HIT: passes, and >= half the path's frames that hold the mover lie within NearMM of its centre
% (sphere geodesic, tile centre to true centre). Speed ratios are for hits only.
%
% X (trackfactorial.csv): one row per hemisphere x cell x speed x rep x depth x (ruler, matching):
%   head signal background matching ruler arm hemi speed rep depth rate tileMM nullFalseS hit speedRatio
% T (rheome.scale.rows, analysis "trackfactorial"): hit_rate and speed_ratio (median over hits) per cell,
% band "<ruler>/<signal>/<background>/<matching>/<head>", pooled over speeds, depths and hemispheres, and
% per depth (band suffix "/d<depth>"); the direct arm as band "direct/<ruler>/<matching>"; n_head_levels.
% The group model (plan G8): speed_ratio ~ the five factors + two-way interactions + (1 | participant).
%
% ⚠ Fixed across cells, and so NOT diagnosed here: blob width (25 mm, Fig. 9's; Fig. 8 used a 15 mm heat
%   kernel), the a-priori frame rate (Fig. 8 read 18.75 fps), and tracking on cortex tiles (Fig. 9 tracked
%   on sphere tiles). The ruler factor includes the sphere-to-cortex distortion by construction.
% ⚠ No RefHead (empty, and RHEOME_REFHEAD unset or not a cache): the head factor runs "own" only and
%   n_head_levels = 1; the group model then has four factors.
%
% See also: rheome.scale.measure_inject, rheome.detect.tilepath, rheome.detect.tilepathnull,
%           sphere_tilepath_validate (nxr-cortical-flow-matlab)
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.RefHead string = ""
        opts.Speeds double = [0.05 0.1 0.25 0.5 1]
        opts.Depths double = [3 6 7 8]
        opts.JLevels double = [0 2 4 6]
        opts.NRep (1,1) double {mustBeInteger, mustBePositive} = 6
        opts.NNull (1,1) double {mustBeInteger, mustBePositive} = 24
        opts.SigmaMM (1,1) double = 25
        opts.PathMM (1,1) double = 120
        opts.KAmp (1,1) double = 2
        opts.NearMM (1,1) double = 50
        opts.Alpha (1,1) double = 0.1
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 41
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    ref = opts.RefHead;  if ref == "", ref = string(getenv('RHEOME_REFHEAD')); end
    if ref ~= "" && ~isfile(fullfile(rheome.load.root(), ref, 'bases.mat')), ref = ""; end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
    chn = string({st.chan.Channel(S.iSel).Name});  clear st
    heads = "own";  if ref ~= "", heads = ["own" "ref"]; end
    X = table();
    for hd = heads
        if hd == "own"
            H0 = struct('B', S.B, 'Sf', S.Sf, 'G', S.G, 'K', S.Res.ImagingKernel, 'iF', 1:numel(chn));
        else
            H0 = i_refhead(ref, chn);
        end
        for hh = opts.Hemis
            rng(opts.Seed + 1000*(hd == "ref") + 100*(hh == "R"));
            X = [X; i_hemi(H0, F(H0.iF, :), fs, hd, hh, opts)]; %#ok<AGROW>
            fprintf('[trackfactorial %s] head %s, %s done\n', name, hd, hh);
        end
    end
    T = i_metrics(X, numel(heads));
end

% The reference head on the participant's channels: gain rows matched by name, its own noise covariance
% on those rows, its own minimum norm (SnrFixed 3).
function H0 = i_refhead(ref, chn)
    Sr = rheome.scale.sensors(char(ref));  st = rheome.load.study(char(ref));
    cr = string({st.chan.Channel(Sr.iSel).Name});  clear st
    [~, ir, iF] = intersect(cr, chn, 'stable');
    if numel(ir) < 0.5*numel(chn), error('scale:trackfactorial:channels', ...
            'RefHead %s shares %d of %d channels with the participant.', ref, numel(ir), numel(chn)); end
    G = Sr.G(ir, :);  ncm = struct('NoiseCov', Sr.ncm.NoiseCov(ir, ir));
    Res = rheome.inverse.mne(G, ncm, struct('ChannelTypes', {Sr.chT(ir)}, 'InverseMeasure', 'amplitude', ...
                                            'nVert', Sr.nV, 'SnrFixed', 3));
    H0 = struct('B', Sr.B, 'Sf', Sr.Sf, 'G', G, 'K', Res.ImagingKernel, 'iF', iF(:)');
end

function X = i_hemi(H0, F, fs, hd, hh, opts)
    FR = 300;  W = 2;  PAD = 2;  F0 = 10;  nW = W*FR;
    Hm = H0.B.(char(hh));  Sh = Hm.S;  Phi = Hm.lbo.Phi;  Mm = Hm.lbo.Mass;  lam = Hm.lbo.Lambda(:);
    av = full(sum(Mm, 2));  gv = double(Hm.gv(:));  nV = numel(gv);
    if isempty(H0.Sf.Sphere), error('scale:trackfactorial:noSphere', 'The surface has no Reg.Sphere.'); end
    Vs = H0.Sf.Sphere(gv, :);  R = mean(vecnorm(Vs, 2, 2));  Vu = Vs ./ vecnorm(Vs, 2, 2);
    nrm = double(Sh.VertNormals);  nrm = nrm ./ vecnorm(nrm, 2, 2);
    rows = reshape((gv' - 1)*3 + (1:3)', [], 1);  K = H0.K(rows, :);  GL = H0.G(:, rows);
    Gn = GL(:, 1:3:end).*nrm(:, 1)' + GL(:, 2:3:end).*nrm(:, 2)' + GL(:, 3:3:end).*nrm(:, 3)';
    hN = exp(-lam*0.015^2/2);  GP = Gn*(Phi.*hN');  rho = exp(-1/(fs*(1/9.71)));
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    nSeg = (W + 2*PAD)*fs;  tSeg = (0:nSeg-1)/fs - PAD;  cIdx = PAD*fs + round((0:nW-1)*fs/FR) + 1;
    nWr = floor(size(F, 2)/(W*fs));
    ge = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=Hm.lbo.L, M=Mm, MaxDepth=max(opts.Depths));
    nd = numel(opts.Depths);  Gs = cell(1, nd);  Wd = Gs;  Cu = Gs;  jR = zeros(numel(opts.Speeds), nd);
    for di = 1:nd
        Gs{di} = rheome.geom.tiles(Tr, Sh, opts.Depths(di), Ruler=ge);  P = double(Gs{di}.P);
        Wd{di} = (P'*(av.*Phi))./(P'*av);  Cu{di} = Vu(Gs{di}.centre, :);
        dT = 1e-3*median(Gs{di}.diameterMM);
        for vi = 1:numel(opts.Speeds)          % a priori: about half a tile per frame (Fig. 9's rule)
            [~, jR(vi, di)] = min(abs(log2(FR./2.^opts.JLevels) - log2(opts.Speeds(vi)/(0.5*dT))));
        end
    end
    % the instrument: analytic sensors at the window's frames -> |J| -> LB modes; direct: the source map
    modes = @(A) Phi'*(Mm*A);
    inst = @(Za) modes(rheome.flow.activation(K*Za));
    realw = @() randi([2, nWr - 3]);
    seg = @(w) F(:, (w-1)*W*fs - PAD*fs + (1:nSeg));
    % backgrounds: analytic window [nCh x nW] at signal level sg, given a real window index or synthetic
    function Za = bgA(kind, sg, w)
        if kind == "real"
            Fs = seg(w);  Za = hilbert(filtfilt(bp, Fs.')).';  Za = Za(:, cIdx);  return
        end
        Z = i_ar1(size(Phi, 2), nSeg, rho);
        if sg == "carrier", Za = real(GP*(Z.*exp(1i*2*pi*F0*tSeg)));  Za = hilbert(filtfilt(bp, Za.')).';  Za = Za(:, cIdx);
        else, Za = GP*Z(:, cIdx); end
    end
    % scales: real background's median per-window envelope peak r0; synthetic matched to it; unit movers
    r0 = median(arrayfun(@(q) max(mean(rheome.flow.activation(K*bgA("real", "carrier", realw())), 2)), 1:10));
    sc.synth.carrier = r0/median(arrayfun(@(q) max(mean(rheome.flow.activation(K*bgA("synth", "carrier", 0)), 2)), 1:6));
    sc.synth.smooth  = r0/median(arrayfun(@(q) max(mean(rheome.flow.activation(K*bgA("synth", "smooth", 0)), 2)), 1:6));
    sc.real.carrier = 1;  sc.real.smooth = 1;
    blob = @(p) exp(-(R*acos(max(min(Vu*p(:), 1), -1))).^2/(2*(opts.SigmaMM*1e-3)^2));
    u = zeros(8, 2);  tt = (0:nSeg-1)/fs;
    for q = 1:8
        b = Gn*blob(Vu(randi(nV), :));
        u(q, 1) = max(mean(rheome.flow.activation(K*i_hilb(bp, b*cos(2*pi*F0*tt), cIdx)), 2));
        u(q, 2) = max(rheome.flow.activation(K*b, Method="norm"));
    end
    amp.carrier = opts.KAmp*r0/median(u(:, 1));  amp.smooth = opts.KAmp*r0/median(u(:, 2));
    cells = table(["carrier"; "carrier"; "smooth"; "smooth"; "smooth"], ["real"; "synth"; "real"; "synth"; "direct"], ...
                  'VariableNames', {'signal', 'background'});
    if hd == "ref", cells = cells(cells.background ~= "direct", :); end
    X = table();
    for c = 1:height(cells)
        sg = cells.signal(c);  bk = cells.background(c);  isDirect = bk == "direct";
        % the map of one window: mover sensors b300 [nCh x nW] (or source X300 [nV x nW]) + its background
        if isDirect
            mapw = @(X300) modes(abs(opts.KAmp*X300 + i_nsrc(Phi, hN, rho, nSeg, cIdx)));
        elseif sg == "carrier"
            mapw = @(b300, car) inst(i_hilb(bp, i_bgraw(bk, seg, realw, GP, rho, nSeg, F0, tSeg, sc.(bk).carrier) ...
                                    + amp.carrier*i_upsample(b300, cIdx, nSeg).*car, cIdx));
        else
            mapw = @(b300) inst(sc.(bk).smooth*bgA(bk, "smooth", realw()) + amp.smooth*b300);
        end
        % null: mover-free windows, per depth and the rates the rule picks
        NUL = cell(nd, numel(opts.JLevels));  Cn = cell(1, opts.NNull);
        for k = 1:opts.NNull
            if isDirect, Cn{k} = mapw(zeros(nV, nW));
            elseif sg == "carrier", Cn{k} = mapw(zeros(size(Gn, 1), nW), zeros(1, nSeg));
            else, Cn{k} = mapw(zeros(size(Gn, 1), nW)); end
        end
        for di = 1:nd
            for jj = unique(jR(:, di))'
                Yn = cell2mat(cellfun(@(C) i_block(Wd{di}*C, opts.JLevels(jj)), Cn, 'uni', 0));
                NUL{di, jj} = rheome.detect.tilepathnull(Yn, Gs{di}, size(Yn, 2)/opts.NNull, FR/2^opts.JLevels(jj), ...
                                                         Alpha=opts.Alpha, WindowS=W);
            end
        end
        Cn = []; %#ok<NASGU>
        for vi = 1:numel(opts.Speeds)
            v = opts.Speeds(vi);  Lp = min(opts.PathMM*1e-3, v*W);  nTr = max(2, round(Lp/v*FR));
            for r = 1:opts.NRep
                [at, c0, b0] = i_greatcircle(R);  off = randi(nW - nTr + 1) - 1;
                sT = nan(1, nW);  sT(off + (1:nTr)) = v*(0:nTr-1)/FR;  pres = ~isnan(sT);
                pT = nan(nW, 3);  pT(pres, :) = at(sT(pres));
                X300 = zeros(nV, nW);  for t = find(pres), X300(:, t) = blob(pT(t, :)); end
                [~, vS] = max(Vu*pT(find(pres, 1), :)');  [~, vE] = max(Vu*pT(find(pres, 1, 'last'), :)');
                vFold = distances(ge, vS, vE)/((nTr - 1)/FR);            % the true NET speed on the folded cortex
                dense = at((0:2e-4:Lp)');
                if isDirect, C = mapw(X300);
                elseif sg == "carrier", C = mapw(Gn*X300, cos(2*pi*F0*tSeg + 2*pi*rand));
                else, C = mapw(Gn*X300); end
                for di = 1:nd
                    jj = jR(vi, di);  bs = 2^opts.JLevels(jj);  rate = FR/bs;  Y = i_block(Wd{di}*C, opts.JLevels(jj));
                    fm = (0:size(Y, 2)-1)*bs + ceil(bs/2);  Nl = NUL{di, jj};
                    p = rheome.detect.tilepath(Y, Gs{di}, Baseline=Nl.baseline, SampleRate=rate, NumPaths=3);
                    sc4 = i_score(p, Nl.threshold, pT(fm, :), Gs{di}, Cu{di}, R, rate, v, vFold, dense, opts.NearMM);
                    base = table(string(hd), sg, bk, "instrument", string(hh), v, r, opts.Depths(di), rate, ...
                                 median(Gs{di}.diameterMM), Nl.falseTestS, 'VariableNames', ...
                                 {'head','signal','background','arm','hemi','speed','rep','depth','rate','tileMM','nullFalseS'});
                    if isDirect, base.arm = "direct";  base.background = "synth"; end
                    X = [X; [repmat(base, 4, 1), sc4]]; %#ok<AGROW>
                end
            end
        end
    end
end

% [matching ruler hit speedRatio] x 4: best / top3 x fold / proj
function s = i_score(p, thr, pT, G, Cu, R, rate, v, vFold, dense, nearMM)
    m = ["best"; "best"; "top3"; "top3"];  ru = ["fold"; "proj"; "fold"; "proj"];
    s = table(m, ru, zeros(4, 1), nan(4, 1), 'VariableNames', {'matching','ruler','hit','speedRatio'});
    if p.nTracks == 0, return; end
    pres = ~any(isnan(pT), 2);  near = zeros(p.nTracks, 1);
    for q = 1:p.nTracks
        f = p.tracks(q).frames(:);  k = p.tracks(q).tiles(:);  fp = pres(f);
        if any(fp)
            near(q) = mean(R*acos(max(min(sum(Cu(k(fp), :).*pT(f(fp), :), 2), 1), -1)) <= nearMM*1e-3);
        end
    end
    ok = [p.tracks.score]' > thr & near >= 0.5;
    qq = [1, 0];  if any(ok), [~, qq(2)] = max(near.*ok); end
    for i = 1:2
        q = qq(i);  if q == 0 || ~ok(q), continue; end
        one = p;  one.tracks = p.tracks(q);  one.nTracks = 1;  st = rheome.detect.trackstats(one, G, rate);
        f = p.tracks(q).frames(:);  k = p.tracks(q).tiles(:);  fp = pres(f);
        [~, j] = max(Cu(k(fp), :)*dense', [], 2);  tf = f(fp)/rate;  pr = NaN;
        if numel(unique(tf)) >= 3, c = polyfit(tf, (j - 1)*2e-4, 1);  pr = c(1)/v; end
        s.hit(2*i - 1:2*i) = 1;  s.speedRatio(2*i - 1) = st.speedMS/vFold;  s.speedRatio(2*i) = pr;
    end
end

function T = i_metrics(X, nHeads)
    T = rheome.scale.rows("trackfactorial", "n_head_levels", nHeads, "levels");
    r = @(m, v, u, b) rheome.scale.rows("trackfactorial", m, v, u, b);
    I = X(X.arm == "instrument", :);
    [g, key] = findgroups(I(:, {'ruler','signal','background','matching','head'}));
    for i = 1:height(key)
        q = I(g == i, :);  b = strjoin(string(key{i, :}), "/");
        T = [T; r(["hit_rate" "speed_ratio"], [mean(q.hit) median(q.speedRatio(q.hit == 1), 'omitnan')], ["fraction" "ratio"], b)]; %#ok<AGROW>
        for d = unique(q.depth)'
            qd = q(q.depth == d, :);
            T = [T; r(["hit_rate" "speed_ratio"], [mean(qd.hit) median(qd.speedRatio(qd.hit == 1), 'omitnan')], ...
                      ["fraction" "ratio"], b + "/d" + d)]; %#ok<AGROW>
        end
    end
    D = X(X.arm == "direct", :);
    [g, key] = findgroups(D(:, {'ruler','matching'}));
    for i = 1:height(key)
        q = D(g == i, :);
        T = [T; r(["hit_rate" "speed_ratio"], [mean(q.hit) median(q.speedRatio(q.hit == 1), 'omitnan')], ...
                  ["fraction" "ratio"], "direct/" + strjoin(string(key{i, :}), "/"))]; %#ok<AGROW>
    end
end

% a random great circle on the sphere of radius R: at(s) [n x 3] unit vectors, s arc length (m)
function [at, c, b] = i_greatcircle(R)
    c = randn(1, 3);  c = c/norm(c);  b = cross(c, [0 0 1]);  if norm(b) < 0.1, b = cross(c, [1 0 0]); end
    b = b/norm(b);  at = @(s) cos(s(:)/R).*c + sin(s(:)/R).*b;
end

% complex AR(1) per mode, unit variance, correlation rho per sample
function Z = i_ar1(n, nT, rho)
    Z = complex(randn(n, nT), randn(n, nT))/sqrt(2);
    for t = 2:nT, Z(:, t) = rho*Z(:, t-1) + sqrt(1 - rho^2)*Z(:, t); end
end

% synthetic source-space background for the direct arm: |.| is taken with the mover, unit median window peak
function N = i_nsrc(Phi, hN, rho, nSeg, cIdx)
    Z = i_ar1(size(Phi, 2), nSeg, rho);  N = Phi*(hN.*Z(:, cIdx));  N = N/median(max(abs(N), [], 1));
end

% the raw (real-valued) background segment at the sensors, for the carrier level
function B = i_bgraw(bk, seg, realw, GP, rho, nSeg, F0, tSeg, s)
    if bk == "real", B = seg(realw());
    else, B = s*real(GP*(i_ar1(size(GP, 2), nSeg, rho).*exp(1i*2*pi*F0*tSeg))); end
end

function Za = i_hilb(bp, B, cIdx)
    Za = hilbert(filtfilt(bp, B.')).';  Za = Za(:, cIdx);
end

% window frames at 300/s -> the segment's samples (linear in time; zero outside the window)
function B = i_upsample(b300, cIdx, nSeg)
    B = zeros(size(b300, 1), nSeg);  B(:, cIdx(1):cIdx(end)) = interp1(cIdx, b300.', cIdx(1):cIdx(end)).';
end

function Y = i_block(Y0, j)
    bs = 2^j;  nF = floor(size(Y0, 2)/bs);
    Y = reshape(mean(reshape(Y0(:, 1:nF*bs), size(Y0, 1), bs, nF), 2), size(Y0, 1), nF);
end

% Author: Diellor Basha, 2026
