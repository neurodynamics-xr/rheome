function [T, X] = measure_aperiodic(name, S, opts)
% SCALE.MEASURE_APERIODIC  A moving change of the aperiodic (1/f) background planted through this participant's
% MEG: when does it read as a wave, and does the periodic split remove it? (MS1 G10, Fig. 14; section 10.1)
%
%   [T, X] = rheome.scale.measure_aperiodic(name)
%   [T, X] = rheome.scale.measure_aperiodic(name, S, SigmaMM=[30 60 100], Reps=8, NullWindows=32)
%
% The group version of nxr-cortical-flow-matlab aperiodic_wave_omega.m (180bbaf; sub-0002, a 30 mm patch, the
% existence case), at lobar scale and with the threshold and the split's selectivity per participant:
%   PLANT    a Gaussian patch (sigma SigmaMM, geodesic) whose centre moves along a PathMM shortest path at v m/s,
%            centred in a WindowS window with cosine ramps (v = 0: the patch sits at mid-path for the whole window,
%            a regional difference). With g the patch gain (1 at its core) and N0, N0', N1 independent smooth 1/f
%            fields (rheome.spectral.synthesise per Laplace-Beltrami mode, heat-smoothed over 15 mm):
%              offset   X = sqrt(OffsetX - 1) g N0'                   core power x OffsetX, exponent unchanged
%              exponent X = (sqrt(1-g) - 1) N0 + sqrt(g) N1           core exponent chi0 -> chi0 + ExponentChange
%              control  X = g c (cos 2pi 12t + cos 2pi 24t)           a coherent moving RHYTHM, 0 dB in Band: the
%                                                                     positive control the split must keep
%   ARMS     emptyroom  X forwarded (normal currents through the gain) onto a synthetic cortical 1/f background
%                       (exponent chi0 = median resting sensor exponent, 2-40 Hz; level = resting aperiodic power
%                       less the empty room's) + the participant's OWN EMPTY ROOM
%            rest       the participant's OWN RESTING recording + the forwarded offset increment (exponent mode needs
%                       a source background to replace, so it is emptyroom only)
%   SPLIT    rheome.spectral.decompose at the sensors on the 6 s window's FFT, the background fitted (knee, 1-45 Hz)
%            on a PSD over FitWindowsS seconds around the window: 6 = the window itself (the fit sees the plant at
%            full weight), 30, Inf = the whole record (the plant diluted to 6/L of the fit)
%   READOUTS Band (8-16 Hz) FFT-mask analytic signal -> MNE read along normals smoothed over 20 mm -> envelope at
%            25 frames/s, then on the total and on each periodic split:
%              viterbi   rheome.detect.tilepath at each of Depths (3, 7) against rheome.detect.tilepathnull on the
%                        plant-free windows (Alpha false paths/s; odd windows calibrate); "on the patch" = the best
%                        passing path within sqrt(2) sigma of the moving centre for >= half the frames it exists
%              bands     rheome.differential.helmholtzbands of the mean apparent velocity over the crossing: the
%                        scalar-potential and stream-function extrema in the bands nearest BandMM (130, 183 mm),
%                        distance to the mid-path (mm); plant and its paired null
%              flow      rheome.flow.apparent (Horn-Schunck, Alpha 0.01, the validated value): the signed along-chord
%                        speed within 2 sigma of the centre (m/s); plant and paired null
%   NULL     every plant window is a background window + the plant; the same window without it is its paired null.
%            The HELD-OUT null on-patch rate of a cell: each plant geometry applied to every even (test) null window.
%
% Returns T (rheome.scale.rows, analysis "aperiodic"):
%   threshold_db   per arm x sigma x speed x signal x depth (band "<arm>_s<sigma>_v<v>_<signal>_d<depth>"): the
%                  smallest in-band dB change (10 log10(1 + plant/background), Band, while the patch is there) of
%                  the offset ladder at which the on-patch rate exceeds the participant's held-out null 95th
%                  percentile; NaN = never exceeded
%   selectivity    per arm x fit window (band "<arm>_<signal>"): median split reduction of the offset plants over
%                  median reduction of the coherent control; reduction = 1 - in-band power of the plant after the
%                  split's gains / before (the gains fitted on plant + background). 1 = not selective
%   reduction      the two medians (band "<arm>_<signal>_offset" | "_control")
% and X, one row per cell (arm, mode, level, sigma, speed, signal), medians over Reps: in-band plant dB and dB
% change, reduction, along-chord speed and null, per depth the passing and on-patch rates (/s) and the held-out
% null 95th percentile (/s), per BandMM the phi / psi extremum errors and their nulls.
%
% ⚠ One hemisphere is planted (Hemi); the emptyroom arm has no sources elsewhere (the rest arm has the real ones).
% ⚠ With an empty room on fewer channels, every arm uses the shared channels and its own MNE kernel (SnrFixed 3).
%   A participant without an empty room runs the rest arm only.
% ⚠ The fit-window PSD mixes the window's own Welch PSD at weight 6/L with the background's over the L span (the
%   rest record; or the empty room + the cortical background's mean PSD) -- the plant exists only in the window,
%   so this is the record's Welch PSD up to the segments straddling the window's edges.
%
% See also: rheome.spectral.decompose, rheome.detect.tilepath, rheome.differential.helmholtzbands,
%           rheome.scale.measure_movingvortex
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Arms string = ["emptyroom" "rest"]
        opts.SigmaMM double = [30 60 100]
        opts.OffsetX double = [2 10 30 100]
        opts.ExponentChange double = -0.5
        opts.SpeedsMS double = [0 0.05 0.15 0.5]
        opts.FitWindowsS double = [6 30 Inf]
        opts.Depths double = [3 7]
        opts.BandMM double = [130 183]
        opts.Band (1,2) double = [8 16]
        opts.Reps (1,1) double {mustBeInteger, mustBePositive} = 8
        opts.NullWindows (1,1) double {mustBeInteger, mustBePositive} = 32
        opts.Hemi (1,1) string = "L"
        opts.PathMM (1,1) double = 100
        opts.WindowS (1,1) double = 2
        opts.MarginS (1,1) double = 2
        opts.Alpha (1,1) double = 0.1
        opts.Seed (1,1) double = 41
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    FITR = [1 45];  FR = 100;  FE = 25;  W = opts.WindowS;  SIG = opts.SigmaMM/1e3;  PATH = opts.PathMM/1e3;

    %% the cortex, the instrument, the records
    Hm = S.B.(char(opts.Hemi));  Sh = Hm.S;  lbo = Hm.lbo;  Phi = lbo.Phi;  lam = lbo.Lambda(:);  Mm = lbo.Mass;
    gv = double(Hm.gv(:));  nV = size(Sh.Vertices, 1);  av = full(sum(Mm, 2));
    Nl = Sh.VertNormals ./ max(vecnorm(Sh.VertNormals, 2, 2), eps);  Sh.VertNormals = Nl;
    st = rheome.load.study(name);  fs = st.rec.sfreq;  REST = double(st.rec.F(S.iSel, :));  clear st
    E = [];  iC = (1:numel(S.iSel))';
    if ismember("emptyroom", opts.Arms)
        try, [E, iC] = i_emptyroom(name, S);
        catch e, fprintf('[aperiodic %s] no empty-room arm: %s\n', name, e.message);
        end
    end
    arms = opts.Arms;  if isempty(E), arms = setdiff(arms, "emptyroom", 'stable'); end
    Res = S.Res;
    if numel(iC) < numel(S.iSel)                         % the shared channels get their own kernel
        Res = rheome.inverse.mne(S.G(iC, :), struct('NoiseCov', S.ncm.NoiseCov(iC, iC)), ...
            struct('ChannelTypes', {S.chT(iC)}, 'InverseMeasure', 'amplitude', 'nVert', S.nV, 'SnrFixed', 3));
    end
    REST = REST(iC, :);  ER = [];
    if ~isempty(E)
        ER = E.F;  if E.fs ~= fs, [pq, qq] = rat(fs/E.fs);  ER = resample(ER.', pq, qq).'; end
    end
    Ns = Phi*(exp(-lam*0.020^2/2) .* (Phi.'*(Mm*Nl)));  Ns = Ns ./ vecnorm(Ns, 2, 2);   % normals smoothed over 20 mm
    Kr = Res.ImagingKernel;  Kn = Ns(:,1).*Kr(3*gv-2,:) + Ns(:,2).*Kr(3*gv-1,:) + Ns(:,3).*Kr(3*gv,:);  clear Res Kr
    Gs = S.G(iC, :);  Gn = Gs(:,3*gv-2).*Nl(:,1).' + Gs(:,3*gv-1).*Nl(:,2).' + Gs(:,3*gv).*Nl(:,3).';  clear Gs
    nT = round((W + 2*opts.MarginS)*fs);  keep = round(opts.MarginS*fs) + (1:round(W*fs));  Lw = nT/fs;
    nseg = round(2*fs);  wfun = hann(nseg);
    [Pr, fw] = pwelch(REST.', wfun, nseg/2, nseg, fs);
    m = fw >= 2 & fw <= 40;  apr = rheome.spectral.aperiodic(Pr(m,:), fw(m));  chi0 = median(apr.exponent);
    mb = fw >= FITR(1) & fw <= FITR(2);  apk = rheome.spectral.aperiodic(Pr(mb,:), fw(mb), struct('knee', true));
    target = trapz(fw(mb), sum(apk.Pap, 2));
    if ~isempty(E)
        Pe = pwelch(ER.', wfun, nseg/2, nseg, fs);  pe = trapz(fw(mb), sum(Pe(mb,:), 2));
        target = max(target - pe, 0.1*pe);
    end
    hN = exp(-lam*0.015^2/2);
    field = @(chi) Phi*(hN .* rheome.spectral.synthesise(nT, fs, chi, ...
                               struct('nS', numel(lam), 'power', ones(numel(lam),1), 'band', FITR)));
    scl = sqrt(target / i_inband(Gn*field(chi0), fs, FITR));
    fprintf('[aperiodic %s] %s: %d vertices, %d channels, %g Hz, chi0 %.2f, arms %s\n', name, opts.Hemi, nV, ...
        numel(iC), fs, chi0, join(arms, '+'));

    ge = rheome.geom.edgegraph(Sh);  Tr = rheome.geom.tree(Sh, L=lbo.L, M=Mm, MaxDepth=max(opts.Depths));
    nd = numel(opts.Depths);  G = cell(1, nd);  Wt = G;
    for k = 1:nd
        G{k} = rheome.geom.tiles(Tr, Sh, opts.Depths(k), Ruler=ge);  Pt = double(G{k}.P);
        Wt{k} = (Pt.'.*av.') ./ (Pt.'*av);
    end
    C = struct('Kn', Kn, 'S', Sh, 'lbo', lbo, 'fs', fs, 'FR', FR, 'FE', FE, 'keep', keep, 'nV', nV, ...
               'band', opts.Band, 'BandMM', opts.BandMM);  C.Wt = Wt;
    nFr = round(W*FE);

    %% the paths, shared by the arms
    rng(opts.Seed);  paths = cell(opts.Reps, 1);
    for r = 1:opts.Reps
        cand = [];  tries = 0;
        while isempty(cand)
            s0 = randi(nV);  d0 = distances(ge, s0);  cand = find(abs(d0 - PATH) < 2e-3);  tries = tries + 1;
            if tries > 200, error('scale:aperiodic:path', 'no %g mm geodesic on this hemisphere', opts.PathMM); end
        end
        [pv, ~, eg] = shortestpath(ge, s0, cand(randi(numel(cand))));
        c = [0; cumsum(ge.Edges.Weight(eg))];  ok = [true; diff(c) > 0];
        paths{r} = struct('v', pv(ok), 'c', c(ok), 'D', distances(ge, pv(ok)));
    end

    sigs = ["total" "p" + i_wlab(opts.FitWindowsS)];
    plants = [struct('mode', "offset", 'level', num2cell(opts.OffsetX)), ...
              struct('mode', "exponent", 'level', num2cell(opts.ExponentChange)), struct('mode', "control", 'level', 0)];
    rows = {};
    for arm = arms
        rng(opts.Seed + find(arm == ["emptyroom" "rest"]));
        if arm == "rest", rec = REST; else, rec = ER; end
        nSlot = floor(size(rec, 2)/fs/Lw);
        if nSlot < opts.NullWindows
            error('scale:aperiodic:short', '%s: %d disjoint %g s windows, need %d', arm, nSlot, Lw, opts.NullWindows);
        end
        ix0 = round(Lw*(sort(randperm(nSlot, opts.NullWindows)) - 1)*fs);   % disjoint windows
        BG = cell(opts.NullWindows, 1);  N0c = BG;  Pctx = 0;
        for w = 1:opts.NullWindows
            ix = ix0(w) + (1:nT);
            if arm == "rest", BG{w} = REST(:, ix);
            else
                N0c{w} = scl*field(chi0);  Bc = Gn*N0c{w};  BG{w} = Bc + ER(:, ix);
                Pctx = Pctx + pwelch(Bc.', wfun, nseg/2, nseg, fs)/opts.NullWindows;   % the cortex's stationary PSD
            end
        end
        fit = @(Y, w) i_fitpsd(Y, rec, ix0(w), nT, fs, wfun, nseg, opts.FitWindowsS, Pctx);

        % --- the null chains: every window, every signal; Viterbi calibrated per signal x depth
        NC = cell(opts.NullWindows, numel(sigs));
        for w = 1:opts.NullWindows
            Ys = i_split(BG{w}, fs, fit(BG{w}, w), fw, FITR);
            for si = 1:numel(sigs), NC{w, si} = i_chain(Ys{si}, C, w <= opts.Reps); end
        end
        cal = cell(numel(sigs), nd);  NP = cell(numel(sigs), nd);  test = 2:2:opts.NullWindows;
        for si = 1:numel(sigs)
            for k = 1:nd
                Yn = cell2mat(cellfun(@(c) c.tile{k}, NC(:, si).', 'uni', 0));
                cal{si,k} = rheome.detect.tilepathnull(Yn, G{k}, nFr, FE, Alpha=opts.Alpha, WindowS=W);
                NP{si,k} = arrayfun(@(w) rheome.detect.tilepath(NC{w,si}.tile{k}, G{k}, ...
                           Baseline=cal{si,k}.baseline, SampleRate=FE), test, 'uni', 0);   % held-out null best paths
            end
        end
        fprintf('[aperiodic %s] %s: null calibrated on %d windows\n', name, arm, opts.NullWindows);

        % --- the plants
        pl = plants;  if arm == "rest", pl = plants([plants.mode] ~= "exponent"); end
        for r = 1:opts.Reps
            pp = paths{r};  dMid = pp.D(ceil(numel(pp.v)/2), :).';
            for s = 1:numel(SIG)
                for v = opts.SpeedsMS
                    [g, cpos, tang, pres] = i_patch(pp, v, SIG(s), nT, fs, keep, W, PATH, Sh.Vertices, FE);
                    pc = max(g, [], 1) > 0.5*max(g, [], 'all');      % the samples while the patch is there
                    for ip = 1:numel(pl)
                        switch pl(ip).mode
                            case "offset",   Xp = sqrt(pl(ip).level - 1) * g .* (scl*field(chi0));
                            case "exponent", Xp = (sqrt(1-g) - 1).*N0c{r} + sqrt(g).*(scl*field(chi0 + pl(ip).level));
                            case "control"
                                tt = (0:nT-1)/fs;  Xp = g .* (cos(2*pi*12*tt) + cos(2*pi*24*tt));
                                Xp = Xp * sqrt(i_bpow(BG{r}, fs, opts.Band, pc) / i_bpow(Gn*Xp, fs, opts.Band, pc));
                        end
                        Bp = Gn*Xp;  Y = BG{r} + Bp;
                        pdB = 10*log10(i_bpow(Bp, fs, opts.Band, pc) / i_bpow(BG{r}, fs, opts.Band, pc));
                        [Ys, Hp] = i_split(Y, fs, fit(Y, r), fw, FITR);
                        Fb = fft(Bp, [], 2);  e0 = i_inbandF(Fb, fs, opts.Band);
                        for si = 1:numel(sigs)
                            red = 0;  if si > 1, red = 1 - i_inbandF(Fb.*Hp{si-1}, fs, opts.Band)/e0; end
                            Q = i_chain(Ys{si}, C, true);  R = i_read(Q, NC{r,si}, cpos, tang, pres, SIG(s), dMid, C);
                            row = struct('arm', arm, 'mode', pl(ip).mode, 'level', pl(ip).level, ...
                                'sigmaMM', opts.SigmaMM(s), 'vMS', v, 'rep', r, 'signal', sigs(si), 'plantDB', pdB, ...
                                'changeDB', 10*log10(1 + 10^(pdB/10)), 'reduction', red, ...
                                'alongMS', R.along, 'alongNullMS', R.alongNull);
                            for k = 1:nd
                                dk = "d" + opts.Depths(k);
                                p = rheome.detect.tilepath(Q.tile{k}, G{k}, Baseline=cal{si,k}.baseline, SampleRate=FE);
                                [row.("pass_" + dk), row.("on_" + dk)] = ...
                                    i_onpatch(p, cal{si,k}.threshold, G{k}, cpos, pres, Sh, SIG(s));
                                [~, nullOn] = cellfun(@(q) i_onpatch(q, cal{si,k}.threshold, G{k}, cpos, pres, Sh, ...
                                                      SIG(s)), NP{si,k});
                                row.("nullOnFrac_" + dk) = mean(nullOn);
                            end
                            for b = 1:numel(opts.BandMM)
                                lb = "_" + opts.BandMM(b);
                                row.("phiErrMM" + lb) = R.phi(b);  row.("psiErrMM" + lb) = R.psi(b);
                                row.("phiErrNullMM" + lb) = R.phiN(b);  row.("psiErrNullMM" + lb) = R.psiN(b);
                            end
                            rows{end+1, 1} = row; %#ok<AGROW>
                        end
                    end
                end
            end
            fprintf('[aperiodic %s] %s: rep %d/%d\n', name, arm, r, opts.Reps);
        end
    end
    [T, X] = i_metrics(struct2table(vertcat(rows{:}), 'AsArray', true), opts, W);
end

%% ------------------------------------------------------------------------------------------------------------------
function [T, X] = i_metrics(Xr, opts, W)
% per cell: medians over reps, rates per second, the held-out null 95th percentile; then threshold and selectivity
    [gi, X] = findgroups(Xr(:, {'arm','mode','level','sigmaMM','vMS','signal'}));
    med = @(c) splitapply(@(x) median(x, 'omitnan'), Xr.(c), gi);
    for c = ["plantDB" "changeDB" "reduction" "alongMS" "alongNullMS"], X.(c) = med(c); end
    X.reps = splitapply(@numel, Xr.rep, gi);
    rng(opts.Seed);
    for k = opts.Depths
        dk = "d" + k;
        X.("passRate_" + dk) = splitapply(@mean, Xr.("pass_" + dk), gi) / W;
        X.("onRate_" + dk) = splitapply(@mean, Xr.("on_" + dk), gi) / W;
        % the cell rate under the null: each rep a Bernoulli draw at its geometry's held-out on-patch fraction
        X.("null95_" + dk) = splitapply(@(p) prctile(mean(rand(4000, numel(p)) < p(:)', 2), 95), ...
                                        Xr.("nullOnFrac_" + dk), gi) / W;
    end
    for b = opts.BandMM
        for c = ["phiErrMM_" "psiErrMM_" "phiErrNullMM_" "psiErrNullMM_"], X.(c + b) = med(c + b); end
    end
    T = table();  r = @(m, v, u, b) rheome.scale.rows("aperiodic", m, v, u, b);
    off = X(X.mode == "offset", :);
    [gc, Cc] = findgroups(off(:, {'arm','sigmaMM','vMS','signal'}));
    for i = 1:height(Cc)
        o = sortrows(off(gc == i, :), 'changeDB');
        for k = opts.Depths
            dk = "d" + k;  hit = find(o.("onRate_" + dk) > o.("null95_" + dk), 1);  th = NaN;
            if ~isempty(hit), th = o.changeDB(hit); end
            T = [T; r("threshold_db", th, "dB", Cc.arm(i) + "_s" + Cc.sigmaMM(i) + "_v" + Cc.vMS(i) + "_" + ...
                      Cc.signal(i) + "_" + dk)]; %#ok<AGROW>
        end
    end
    for a = unique(X.arm)'
        for sg = unique(X.signal(X.signal ~= "total"))'
            k = X.arm == a & X.signal == sg;
            ro = median(X.reduction(k & X.mode == "offset"), 'omitnan');
            rc = median(X.reduction(k & X.mode == "control"), 'omitnan');
            T = [T; r(["selectivity" "reduction" "reduction"], [ro/rc ro rc], "ratio", ...
                      a + "_" + sg + ["" "_offset" "_control"])]; %#ok<AGROW>
        end
    end
end

function L = i_wlab(w)
    L = string(w) + "s";  L(isinf(w)) = "whole";
end

function P = i_fitpsd(Y, rec, ix0, nT, fs, wfun, nseg, Ls, Pctx)
% the fit PSD for each fit window L: the window's own Welch PSD at weight Lw/L, the background over the L span
    Pw = pwelch(Y.', wfun, nseg/2, nseg, fs);  Lw = nT/fs;  P = cell(1, numel(Ls));  nR = size(rec, 2);
    for j = 1:numel(Ls)
        a = min(Lw/Ls(j), 1);
        if a == 1, P{j} = Pw; continue; end
        if isinf(Ls(j)), span = 1:nR;
        else, n = min(round(Ls(j)*fs), nR);  s0 = min(max(ix0 + round(nT/2) - round(n/2), 0), nR - n);  span = s0 + (1:n);
        end
        P{j} = a*Pw + (1 - a)*(pwelch(rec(:, span).', wfun, nseg/2, nseg, fs) + Pctx);
    end
end

function [Ys, Hs] = i_split(Y, fs, Pfit, ff, FITR)
% the total and, per fit PSD, the periodic part at the sensors (Hermitian gains on the window's FFT)
    nT = size(Y, 2);  Fw = fft(Y, [], 2);  f = (0:nT-1)*fs/nT;  pos = 2:floor(nT/2)+1;  mf = ff >= FITR(1) & ff <= FITR(2);
    k = fs*nT/2;                                              % one-sided Welch PSD -> E|FFT bin|^2 of nT samples
    Ys = [{Y} cell(1, numel(Pfit))];  Hs = cell(1, numel(Pfit));
    mir = nT - pos + 2;  ok = mir > pos & mir <= nT;
    for j = 1:numel(Pfit)
        D = rheome.spectral.decompose(Fw(:, pos), f(pos)', struct('fitPower', (Pfit{j}(mf,:)*k)', 'fitF', ff(mf), 'knee', true));
        Hp = zeros(size(Fw));  Hp(:, pos) = D.hper;  Hp(:, mir(ok)) = D.hper(:, ok);
        Ys{j+1} = real(ifft(Fw.*Hp, [], 2));  Hs{j} = Hp;
    end
end

function p = i_inband(Y, fs, band)
    p = i_inbandF(fft(Y, [], 2), fs, band);
end

function p = i_inbandF(F, fs, band)               % summed over rows: in-band power by Parseval (positive bins x 2)
    nT = size(F, 2);  k = 1:floor(nT/2);  f = k*fs/nT;  in = f >= band(1) & f <= band(2);
    p = 2*sum(abs(F(:, k(in)+1)).^2, 'all') / nT^2;
end

function p = i_bpow(Y, fs, band, cols)          % band-passed power over the columns cols (FFT mask on the window)
    nT = size(Y, 2);  f = (0:nT-1)*fs/nT;  m = (f >= band(1) & f <= band(2)) | (f >= fs-band(2) & f <= fs-band(1));
    Yb = real(ifft(fft(Y, [], 2) .* m, [], 2));  p = sum(Yb(:, cols).^2, 'all');
end

function Q = i_chain(Y, C, full)
% analytic at the sensors, normal MNE current, envelope at FE, tiles per depth (and the apparent flow)
    nT = size(Y, 2);  f = (0:nT-1)*C.fs/nT;
    z = ifft(fft(Y, [], 2) .* (2*(f >= C.band(1) & f <= C.band(2))), [], 2);
    zf = C.Kn * z(:, C.keep(1):round(C.fs/C.FR):C.keep(end));
    bs = C.FR/C.FE;  nF = floor(size(zf, 2)/bs);
    A = squeeze(mean(reshape(abs(zf(:, 1:nF*bs)), C.nV, bs, nF), 2));
    Q.tile = cellfun(@(Wt) Wt*A, C.Wt, 'uni', 0);
    if full
        o = rheome.flow.apparent(A, C.S, Rate=C.FE, Band=C.band, Alpha=0.01);   % ⚠ Alpha 0.01, not the default 1
        Q.vel = single(o.velocity*C.FE);                                         % ⚠ .velocity is metres per FRAME
    end
end

function [g, cpos, tang, pres] = i_patch(pp, v, SIG, nT, fs, keep, W, PATH, V, FE)
% patch gain [nV x nT], and per envelope frame the centre, the path chord and whether the patch is there
    t = ((0:nT-1) - keep(1) + 1)/fs;
    if v == 0, Tc = W; s = PATH/2*ones(size(t)); else, Tc = PATH/v; s = (t - (W - Tc)/2)*v; end
    on = s >= 0 & s <= PATH;  if v == 0, on = t >= 0 & t <= W; end
    rmp = min(0.1, Tc/4);  t1 = t - max((W - Tc)/2, 0);  amp = double(on);
    if v > 0, amp = on .* min(1, min(t1, Tc - t1)/rmp); amp = 0.5 - 0.5*cos(pi*max(min(amp, 1), 0)); end
    if v == 0, amp(:) = 1; end                                     % a regional difference, present throughout
    sc = min(max(s, 0), pp.c(end));  k = discretize(sc, [pp.c(:); Inf]);  k = min(k, numel(pp.c) - 1);
    a = (sc - pp.c(k).') ./ (pp.c(k+1).' - pp.c(k).');
    g = zeros(size(pp.D, 2), nT);
    for j = find(amp > 0)
        d = (1 - a(j))*pp.D(k(j), :) + a(j)*pp.D(k(j)+1, :);  g(:, j) = amp(j)*exp(-d.^2/(2*SIG^2));
    end
    nF = round(W*FE);  fe = min(keep(1) + round((0:nF-1)*fs/FE) + round(fs/(2*FE)), nT);  kk = k(fe);
    cpos = (1 - a(fe)).'.*V(pp.v(kk), :) + a(fe).'.*V(pp.v(kk+1), :);
    % ⚠ the CHORD, not the local tangent: the path folds through sulcal walls; the rendered blob moves along the net
    tang = repmat((V(pp.v(end), :) - V(pp.v(1), :)) / norm(V(pp.v(end), :) - V(pp.v(1), :)), nF, 1);
    pres = amp(fe) > 0.5;
end

function R = i_read(Q, Q0, cpos, tang, pres, SIG, dMid, C)
    V = C.S.Vertices;  fr = find(pres(1:end-1));  if isempty(fr), fr = 1:size(Q.vel, 2); end
    R.along = i_along(Q.vel, fr, V, cpos, tang, SIG);  R.alongNull = i_along(Q0.vel, fr, V, cpos, tang, SIG);
    [R.phi, R.psi] = i_bands(mean(double(Q.vel(:, fr)), 2), C, dMid);
    [R.phiN, R.psiN] = i_bands(mean(double(Q0.vel(:, fr)), 2), C, dMid);
end

function al = i_along(vel, fr, V, cpos, tang, SIG)
    a = nan(numel(fr), 1);
    for j = 1:numel(fr)
        near = vecnorm(V - cpos(fr(j), :), 2, 2) <= 2*SIG;
        U = reshape(double(vel(:, fr(j))), 3, []).';  a(j) = mean(U(near, :)*tang(fr(j), :).');
    end
    al = median(a, 'omitnan');
end

function [ePhi, ePsi] = i_bands(J, C, dMid)
% the extremum of phi and of psi in the band nearest each BandMM, its distance to the mid-path (mm)
    Bd = rheome.differential.helmholtzbands(J, C.S, C.lbo, Maps=true);
    [~, mb] = min(abs(Bd.wavelengthMM(:) - C.BandMM(:)'), [], 1);
    ePhi = zeros(1, numel(mb));  ePsi = ePhi;
    for b = 1:numel(mb)
        ph = Bd.PhiBand(:, mb(b));  ps = Bd.PsiBand(:, mb(b));
        [~, i1] = max(abs(ph - median(ph)));  [~, i2] = max(abs(ps - median(ps)));
        ePhi(b) = 1e3*dMid(i1);  ePsi(b) = 1e3*dMid(i2);
    end
end

function [pass, on] = i_onpatch(pth, thr, G, cpos, pres, S, SIG)
% a passing best path, and whether it follows the moving centre (within sqrt(2) sigma, >= half its frames)
    pass = 0;  on = 0;
    if pth.nTracks == 0 || pth.tracks(1).score <= thr, return; end
    pass = 1;  x = pth.tracks(1);
    [~, vc] = min(pdist2(S.Vertices, cpos), [], 1);  tTru = G.tileOf(vc);
    f = x.frames(:);  ok = pres(f).';  if ~any(ok), return; end
    tl = x.tiles(:);  hit = tl(ok) == tTru(f(ok)) | G.D(sub2ind(size(G.D), tl(ok), tTru(f(ok)))) <= sqrt(2)*SIG;
    on = double(mean(hit) >= 0.5 && nnz(ok) >= 0.5*nnz(pres));
end

function [E, iE] = i_emptyroom(name, S)
% the participant's own empty room on the channels it shares, by name, with S.iSel (as measure_patterns)
    N = builtin('load', fullfile(rheome.load.root(), name, 'noise.mat'), 'nrec');  nrec = N.nrec;
    st = rheome.load.study(name);  nm = string(st.chan.Name(S.iSel));  clear st
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;
    [~, iE, jn] = intersect(nm, nn(ok), 'stable');  jo = find(ok);  jn = jo(jn);
    if isempty(iE), error('scale:aperiodic:noise', 'the empty room shares no good channel with the recording'); end
    E = struct('F', double(nrec.F(jn, :)), 'fs', nrec.sfreq);
end

% Author: Diellor Basha, 2026
