function [T, X] = measure_noisefloor(name, S, item, opts)
% SCALE.MEASURE_NOISEFLOOR  What the noise floor and the inverse manufacture, in this participant (MS1 G2, section 9.4).
%
%   [T, X] = rheome.scale.measure_noisefloor(name, S, "composition")
%   items: "composition" | "sizeruler" | "vortexscale" | "rotation" | "detection" | "diracangles"
%
% One item per call, so rheome.scale.run times each and an error in one does not cost the others. Every
% item runs on the participant's own leadfield, plain minimum norm (S.Res, SnrFixed 3, amplitude) and,
% where a background is needed, the participant's own empty room (noise.mat, channels matched BY NAME),
% on both hemispheres (Hemis). X is the per-unit table, T the rheome.scale.rows (analysis = item).
%
% composition   single-subject reference: vortex 0.036/0.608 irr/sol before, 0.338/0.253 after; spiral 0.347/0.236; empty
%               room 0.307/0.246; rest 0.334/0.262; gap at the focus 0.000; rows 0.427/0.127/0.348.
%               rheome.flow.seedvortex 140 mm, Phase pi/2 (vortex) and pi/4 (spiral), at Seeds random
%               vertices; read before and after leadfield + MNE (noise-free), and NumFrames strongest
%               8-12 Hz moments (real part of the sensor analytic signal) of rest and of the empty room.
%               Shares of |J|^2 (lumped area) in grad(Phi), N x grad(Psi), the normal part: over the
%               hemisphere ("whole") and over the top 10% of vertices by |J| ("focus"). Leadfield rows:
%               the same shares for each sensor's row. Metrics: irr_frac, sol_frac, normal_frac per
%               source x region; sol_gap_vortex_spiral and sol_gap_vortex_emptyroom at the focus.
% sizeruler     single-subject reference: slope +0.52 R2 0.83 on the cortex -> +0.02 / 0.002 at <= 10 dB, plateau ~100 mm;
%               location 35 / 43 / 89 mm at noiseless / 10 / 0 dB, 7 mm without the instrument.
%               Port of the size-ruler planting test on the plain MNE (the reference read it through the Dirac inverse; the
%               paper's estimator is now the plain MNE): heat-kernel blobs of 70-267 mm of NORMAL current
%               at Seeds depth-stratified vertices, noise from the participant's noise covariance at
%               SNR (sensor, dB), read by rheome.detect.blobscale (mexhat, 40-500 mm) at the |J| peak.
%               Metrics per condition: slope, r2, plateau_mm (median recovered size), loc_err_mm.
% vortexscale   single-subject reference: power-weighted wavelength 115-119 mm empty room, 173 mm rest alpha, 128 mm gamma;
%               ~2.5x more vortices at the floor. NumFrames evenly spaced frames of |K x analytic| in
%               alpha (8-12 Hz) and gamma (30-45 Hz), rest and empty room: the power-weighted LBO
%               wavelength 2*pi/sqrt(sum lambda c^2 / sum c^2) (constant mode out) and the count of
%               vortex-type critical points of the real-part field with |J| > 20% of its maximum.
% rotation      single-subject reference: iCoh +0.969 noiseless, +0.960 / +0.894 / +0.596 at 10 / 3 / 1 nAm; standing |iCoh|
%               <= 0.031. Port of the spin-recovery test: a 140 mm vortex and its pi/2-spun partner in
%               quadrature at sqrt(8*16) Hz (rotating) or the vortex alone (standing), Hann-windowed 2 s,
%               scaled to MomentNAm RSS, plus empty-room windows; the MNE estimate projected onto the
%               two planted patterns, 8-16 Hz, imaginary coherence between them.
% detection     single-subject reference: 0.09-0.56 nAm RSS across depth, a 41x range across families, rotor 140 mm 0.3 nAm.
%               rheome.forward.simulate (the threshold against the empty room's own 95th percentile)
%               per family at Seeds median-depth vertices, and a rheome.flow.seedrotor 140 mm rotor
%               through the vertex leadfield. Needs the Dirac basis (imported here when absent).
% diracangles   single-subject reference: rows a median 0.187 inside the 400-mode span, 4 of 270 principal angles with cos >
%               0.9, mean cos^2 0.041; the basis carries 0.026 of a vortex's measurable part. Leadfield
%               rows (hemisphere columns) against the x,y,z rows of the hemisphere's Dirac modes.
%
% See also: rheome.scale.run, rheome.scale.measure_plantfloors, rheome.differential.helmholtz,
%           rheome.detect.blobscale, rheome.forward.simulate
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S
        item (1,1) string {mustBeMember(item, ["composition" "sizeruler" "vortexscale" "rotation" "detection" "diracangles"])}
        opts.Hemis string = ["L" "R"]
        opts.Seeds (1,1) double {mustBeInteger, mustBePositive} = 8
        opts.NumFrames (1,1) double {mustBeInteger, mustBePositive} = 20
        opts.SNRdB double = [Inf 20 10 0]
        opts.MomentNAm double = [Inf 10 3 1]
        opts.Families string = ["oscillator" "resonator" "gabor" "travwave" "wave" "dampedwave" "diffusion"]
        opts.Seed (1,1) double = 37
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    rng(opts.Seed);
    X = feval(char("i_" + item), name, S, opts);
    T = feval(char("i_" + item + "_rows"), X, opts);
end

%% ---------- composition
function X = i_composition(name, S, o)
    [R, fs] = i_rest(name, S);  [E, iE] = i_emptyroom(name, S);
    X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  g = rheome.operators.gauge(C.Sh.Vertices, double(C.Sh.Faces), Method="diffusion");
        for k = 1:o.Seeds
            v = i_farvertex(C, g);
            for ty = ["vortex" "spiral"]
                ph = pi/2;  if ty == "spiral", ph = pi/4; end
                sv = rheome.flow.seedvortex(name, Vertex=v, WavelengthMM=140, Phase=ph, Hemi=h, Bases=S.B, Gauge=g, Check=false);
                X = [X; i_shares(sv.J, C, h, ty + "_planted")]; %#ok<AGROW>
                X = [X; i_shares(C.K * (C.G * sv.J), C, h, ty + "_meg")]; %#ok<AGROW>
            end
        end
        X = [X; i_shares(C.K * i_alphaframes(R, fs, o.NumFrames), C, h, "rest")]; %#ok<AGROW>
        X = [X; i_shares(C.K(:, iE) * i_alphaframes(E.F, E.fs, o.NumFrames), C, h, "emptyroom")]; %#ok<AGROW>
        X = [X; i_shares(C.G', C, h, "leadfield_rows")]; %#ok<AGROW>
    end
end

function X = i_shares(J, C, h, src)
% irrotational / solenoidal / normal shares of |J|^2 per column, over the hemisphere and at the focus
    Hh = rheome.differential.helmholtz(J, C.Sh);  a = C.area;
    m2 = @(V) V(1:3:end,:).^2 + V(2:3:end,:).^2 + V(3:3:end,:).^2;
    Jn = J(1:3:end,:).*C.n(:,1) + J(2:3:end,:).*C.n(:,2) + J(3:3:end,:).*C.n(:,3);
    tot = m2(J);  irr = m2(Hh.Virr);  sol = m2(Hh.Vsol);  nrm = Jn.^2;
    X = table();
    for rg = ["whole" "focus"]
        W = repmat(a, 1, size(J, 2));
        if rg == "focus", W = W .* (tot >= prctile(tot, 90, 1)); end
        den = sum(W .* tot, 1);
        X = [X; table(repmat(h, size(J,2), 1), repmat(src, size(J,2), 1), repmat(rg, size(J,2), 1), ...
             (sum(W.*irr,1)./den)', (sum(W.*sol,1)./den)', (sum(W.*nrm,1)./den)', ...
             'VariableNames', {'hemi','source','region','irrFrac','solFrac','normalFrac'})]; %#ok<AGROW>
    end
end

function T = i_composition_rows(X, ~)
    T = table();  r = @(m, v, b) rheome.scale.rows("composition", m, v, "fraction", b);
    for s = unique(X.source, 'stable')'
        for rg = ["whole" "focus"]
            k = X.source == s & X.region == rg;
            T = [T; r(["irr_frac" "sol_frac" "normal_frac"], median([X.irrFrac(k) X.solFrac(k) X.normalFrac(k)], 1), s + "_" + rg)]; %#ok<AGROW>
        end
    end
    f = @(s) median(X.solFrac(X.source == s & X.region == "focus"));
    T = [T; r(["sol_gap_vortex_spiral" "sol_gap_vortex_emptyroom"], ...
              [f("vortex_meg") - f("spiral_meg"), f("vortex_meg") - f("emptyroom")], "focus")];
end

%% ---------- size ruler
function X = i_sizeruler(~, S, o)
    LAM = [70 100 133 189 267];  NREAL = 5;
    NC = S.ncm.NoiseCov;  Lch = chol(NC + 1e-12*trace(NC)/size(NC,1)*eye(size(NC)), 'lower');
    X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  lbo = C.lbo;  P = lbo.Phi;  lam = lbo.Lambda(:);
        dmin = min(pdist2(C.Sh.Vertices, reshape(S.Loc, 3, [])'), [], 2);
        q = discretize(dmin, prctile(dmin, linspace(0, 100, o.Seeds + 1)));
        seeds = arrayfun(@(k) i_pick(find(q == k)), 1:o.Seeds);
        for l = LAM
            hk = exp(-((l*1e-3/(2*pi))^2/2) * lam);
            for v = seeds
                a = P * (hk .* (P' * (lbo.Mass * full(sparse(v, 1, 1, C.nV, 1)))));  a = a / max(abs(a));
                J = reshape((a .* C.n)', [], 1);  b = C.G * J;  d = 1e3 * distances(C.ge, v)';
                X = [X; i_ruler(a, "direct", Inf, 1, l, v, d, dmin, lbo, h)]; %#ok<AGROW>
                for snr = o.SNRdB
                    for rr = 1:(1 + (NREAL - 1) * isfinite(snr))
                        y = b;
                        if isfinite(snr), n = Lch * randn(size(b));  y = b + n * norm(b) / max(norm(n), realmin) * 10^(-snr/20); end
                        Jh = C.K * y;  amp = sqrt(Jh(1:3:end).^2 + Jh(2:3:end).^2 + Jh(3:3:end).^2);
                        X = [X; i_ruler(amp, "meg", snr, rr, l, v, d, dmin, lbo, h)]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

function row = i_ruler(amp, cond, snr, rr, l, v, d, dmin, lbo, h)
    [~, pk] = max(amp);
    bs = rheome.detect.blobscale(amp, lbo, Kernel="mexhat", Wavelengths=[40 500]*1e-3);
    row = table(h, cond, snr, rr, l, v, 1e3*dmin(v), 1e3*bs.wavelength(pk), d(pk), ...
        'VariableNames', {'hemi','condition','snrDB','real','plantedMM','vertex','depthMM','recoveredMM','locErrMM'});
end

function T = i_sizeruler_rows(X, o)
    T = table();  r = @(m, v, u, b) rheome.scale.rows("sizeruler", m, v, u, b);
    cond = ["direct" repmat("meg", 1, numel(o.SNRdB))];  snr = [Inf o.SNRdB];
    for c = 1:numel(cond)
        k = X.condition == cond(c) & X.snrDB == snr(c);
        b = cond(c);  if cond(c) == "meg", b = "meg_" + i_snr(snr(c)); end
        p = polyfit(X.plantedMM(k), X.recoveredMM(k), 1);  rr = corr(X.plantedMM(k), X.recoveredMM(k))^2;
        T = [T; r(["slope" "r2" "plateau_mm" "loc_err_mm"], [p(1) rr median(X.recoveredMM(k)) median(X.locErrMM(k))], ...
                  ["mm/mm" "r2" "mm" "mm"], b)]; %#ok<AGROW>
    end
end

%% ---------- vortex wavelength and count
function X = i_vortexscale(name, S, o)
    [R, fs] = i_rest(name, S);  [E, iE] = i_emptyroom(name, S);
    bands = struct('alpha', [8 12], 'gamma', [30 45]);
    X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  Sg = i_sg(C, name, h);  op = rheome.detect.operator(Sg);
        for src = ["rest" "emptyroom"]
            if src == "rest", Y = R; f = fs; K = C.K; else, Y = E.F; f = E.fs; K = C.K(:, iE); end
            for bn = string(fieldnames(bands))'
                [b, a] = butter(3, bands.(bn) / (f/2), 'bandpass');
                Z = hilbert(filtfilt(b, a, Y.')).';
                fr = round(linspace(5*f, size(Z, 2) - 5*f, o.NumFrames));
                Jc = K * Z(:, fr);  amp = sqrt(abs(Jc(1:3:end,:)).^2 + abs(Jc(2:3:end,:)).^2 + abs(Jc(3:3:end,:)).^2);
                c = C.lbo.Phi(:, 2:end)' * (C.lbo.Mass * amp);  lam = C.lbo.Lambda(2:end);
                wl = 1e3 * 2*pi ./ sqrt(sum(lam .* c.^2, 1) ./ sum(c.^2, 1));
                nv = zeros(1, numel(fr));
                for i = 1:numel(fr)
                    Jr = real(Jc(:, i));  cp = rheome.detect.criticalPoints(Jr, Sg, "vortex", op);
                    if isempty(cp.charge), continue, end
                    ai = sqrt(Jr(1:3:end).^2 + Jr(2:3:end).^2 + Jr(3:3:end).^2);
                    [~, nn] = min(pdist2(cp.pos, C.Sh.Vertices), [], 2);  nv(i) = sum(ai(nn) > 0.2 * max(ai));
                end
                X = [X; table(repmat(h, numel(fr), 1), repmat(src, numel(fr), 1), repmat(bn, numel(fr), 1), ...
                     fr(:)/f, wl(:), nv(:), 'VariableNames', {'hemi','source','band','tS','wavelengthMM','nVortex'})]; %#ok<AGROW>
            end
        end
    end
end

function T = i_vortexscale_rows(X, ~)
    T = table();
    for s = ["rest" "emptyroom"]
        for b = ["alpha" "gamma"]
            k = X.source == s & X.band == b;
            T = [T; rheome.scale.rows("vortexscale", ["wavelength_mm" "n_vortex"], ...
                     [median(X.wavelengthMM(k)) median(X.nVortex(k))], ["mm" "count"], s + "_" + b)]; %#ok<AGROW>
        end
    end
end

%% ---------- rotation (imaginary coherence)
function X = i_rotation(name, S, o)
    [E, iE] = i_emptyroom(name, S);  fs = E.fs;  nT = round(2*fs);  tt = (0:nT-1)/fs;  f0 = sqrt(8*16);
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    env = hann(nT)';  X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  g = rheome.operators.gauge(C.Sh.Vertices, double(C.Sh.Faces), Method="diffusion");
        G = C.G(iE, :);  K = C.K(:, iE);
        for k = 1:o.Seeds
            sv = rheome.flow.seedvortex(name, Vertex=i_farvertex(C, g), WavelengthMM=140, Hemi=h, Bases=S.B, Gauge=g, Check=false);
            JA3 = reshape(sv.J, 3, [])';  JB = reshape(rheome.filters.spin(JA3, g.normal, pi/2)', [], 1);
            PK = pinv([sv.J JB]) * K;  bA = G * sv.J;  bB = G * JB;  rss = max(vecnorm([sv.J JB], 2, 1));   % [2 x nCh]
            starts = randi(size(E.F, 2) - nT, 1, 20);
            for kind = ["rotating" "standing"]
                cB = sin(2*pi*f0*tt) .* env * (kind == "rotating");
                Bs = (bA * (cos(2*pi*f0*tt) .* env) + bB * cB) * (1e-9 / rss);
                for amp = o.MomentNAm
                    ic = nan(1, numel(starts));
                    for j = 1:numel(starts)
                        y = Bs;  if isfinite(amp), y = amp * Bs + E.F(:, starts(j) + (0:nT-1)); end
                        hz = hilbert(filtfilt(bp, (PK * y).')).';
                        c = mean(hz(1,:) .* conj(hz(2,:))) / sqrt(mean(abs(hz(1,:)).^2) * mean(abs(hz(2,:)).^2));
                        ic(j) = imag(c);  if ~isfinite(amp), break, end
                    end
                    X = [X; table(h, k, kind, amp, median(ic, 'omitnan'), ...
                         'VariableNames', {'hemi','seed','kind','momentNAm','iCoh'})]; %#ok<AGROW>
                end
            end
        end
    end
end

function T = i_rotation_rows(X, o)
    T = table();
    for amp = o.MomentNAm
        b = "noiseless";  if isfinite(amp), b = amp + "nAm"; end
        T = [T; rheome.scale.rows("rotation", ["icoh_rotating" "abs_icoh_standing"], ...
                 [median(X.iCoh(X.kind == "rotating" & X.momentNAm == amp)), ...
                  median(abs(X.iCoh(X.kind == "standing" & X.momentNAm == amp)))], "coh", b)]; %#ok<AGROW>
    end
end

%% ---------- detection thresholds
function X = i_detection(name, S, o)
    d = i_dirac(name);  [E, ~] = i_emptyroom(name, S);
    st = rheome.load.study(name);  st.rec = rmfield(st.rec, 'F');               % the forward reads flags only
    er = struct('chan', struct('Name', {E.nrec.ChannelName(:)'}, 'Type', {E.nrec.ChannelType(:)'}), 'rec', E.nrec);
    X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  ctx = struct('H', S.B.(char(h)), 'd', d, 'st', st, 'er', er);
        dmm = 1e3 * min(pdist2(C.Sh.Vertices, reshape(S.Loc, 3, [])'), [], 2);
        mid = find(abs(dmm - median(dmm)) < 5);  seeds = mid(randperm(numel(mid), min(o.Seeds, numel(mid))));
        for v = seeds(:)'
            for fam = [o.Families "rotor140"]
                if fam == "rotor140"
                    r = rheome.flow.seedrotor(name, Vertex=v, WavelengthMM=140, Hemi=h, Bases=S.B, Check=false);
                    s = rheome.forward.simulate(name, Family="oscillator", Band=[8 16], Vertex=v, Pattern=r.J, ...
                            NoiseTrials=40, MomentNAm=1, Hemi=h, Context=ctx, Verbose=false);
                else
                    s = rheome.forward.simulate(name, Family=fam, Band=[8 16], WavelengthMM=140, Vertex=v, ...
                            NoiseTrials=40, MomentNAm=1, Hemi=h, Context=ctx, Verbose=false);
                end
                X = [X; table(h, v, dmm(v), fam, s.threshMomentNAm, ...
                     'VariableNames', {'hemi','vertex','depthMM','family','threshNAm'})]; %#ok<AGROW>
            end
        end
    end
end

function T = i_detection_rows(X, ~)
    f = unique(X.family, 'stable');  m = arrayfun(@(q) median(X.threshNAm(X.family == q)), f);
    T = [rheome.scale.rows("detection", repmat("thresh_nam", numel(f), 1), m, "nAm", f); ...
         rheome.scale.rows("detection", "family_range", max(m) / min(m), "ratio", "")];
end

%% ---------- Dirac-basis principal angles
function X = i_diracangles(name, S, o)
    d = i_dirac(name);  X = table();
    for h = o.Hemis
        C = i_hemi(S, h);  gv = double(S.B.(char(h)).gv(:))';
        cols = find(d.Hemisphere == double(h == "R") + 1);
        D = full(d.Phi((gv' - 1)*4 + (2:4), cols));                      % x, y, z rows of the hemisphere's modes
        D = reshape(permute(reshape(D, numel(gv), 3, []), [2 1 3]), 3*numel(gv), []);   % -> [x1 y1 z1 x2 ...]
        Qd = orth(real(D));  Qg = orth(C.G');
        inside = vecnorm(Qd' * C.G', 2, 1).^2 ./ vecnorm(C.G', 2, 1).^2;
        cs = svd(Qg' * Qd);  cs = cs(1:min(size(Qg, 2), size(Qd, 2)));
        g = rheome.operators.gauge(C.Sh.Vertices, double(C.Sh.Faces), Method="diffusion");
        vm = zeros(1, o.Seeds);
        for k = 1:o.Seeds
            sv = rheome.flow.seedvortex(name, Vertex=i_farvertex(C, g), WavelengthMM=140, Hemi=h, Bases=S.B, Gauge=g, Check=false);
            m = Qg * (Qg' * sv.J);  vm(k) = norm(Qd' * m)^2 / max(norm(m)^2, realmin);
        end
        X = [X; table(h, median(inside), sum(cs > 0.9), numel(cs), mean(cs.^2), median(vm), ...
             'VariableNames', {'hemi','rowsInside','nCosAbove09','nAngles','meanCos2','vortexMeasurableInside'})]; %#ok<AGROW>
    end
end

function T = i_diracangles_rows(X, ~)
    T = rheome.scale.rows("diracangles", ["rows_inside" "n_cos_above_0p9" "n_angles" "mean_cos2" "vortex_measurable_inside"], ...
        median([X.rowsInside X.nCosAbove09 X.nAngles X.meanCos2 X.vortexMeasurableInside], 1), ...
        ["fraction" "count" "count" "cos2" "fraction"]);
end

%% ---------- shared
function C = i_hemi(S, h)
    H = S.B.(char(h));  C.Sh = H.S;  C.lbo = H.lbo;  C.nV = size(H.S.Vertices, 1);
    rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
    C.K = S.Res.ImagingKernel(rows, :);  C.G = S.G(:, rows);
    C.n = H.S.VertNormals ./ max(vecnorm(H.S.VertNormals, 2, 2), eps);
    C.area = full(sum(H.lbo.Mass, 2));  C.ge = rheome.geom.edgegraph(H.S);
end

function [F, fs] = i_rest(name, S)
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
end

function [E, iE] = i_emptyroom(name, S)
% the participant's own empty room on the channels it shares, by name, with S.iSel; iE indexes S's channels
    N = builtin('load', fullfile(rheome.load.root(), name, 'noise.mat'), 'nrec');  nrec = N.nrec;
    st = rheome.load.study(name);  nm = string(st.chan.Name(S.iSel));  clear st
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;
    [~, iE, jn] = intersect(nm, nn(ok), 'stable');  jo = find(ok);  jn = jo(jn);
    if isempty(iE), error('scale:noisefloor:noise', 'the empty room shares no good channel with the recording'); end
    E = struct('F', double(nrec.F(jn, :)), 'fs', nrec.sfreq, 'nrec', nrec);
end

function Y = i_alphaframes(F, fs, nF)
% the real part of the 8-12 Hz sensor analytic signal at the nF strongest moments (>= 1 s apart, 5 s edges)
    [b, a] = butter(3, [8 12] / (fs/2), 'bandpass');  Z = hilbert(filtfilt(b, a, F.')).';
    env = mean(abs(Z), 1);  e = round(5*fs);  env([1:e, end-e+1:end]) = 0;
    [~, fr] = findpeaks(env, 'MinPeakDistance', round(fs), 'SortStr', 'descend', 'NPeaks', nF);
    Y = real(Z(:, fr));
end

function v = i_farvertex(C, g)
% a random vertex >= 20 mm from the gauge's singular faces (the frame is meaningless in their one-ring)
    ok = true(C.nV, 1);  P = C.Sh.Vertices;
    if ~isempty(g.singular)
        Fc = double(C.Sh.Faces(g.singular, :));
        ok = min(pdist2(P, (P(Fc(:,1),:) + P(Fc(:,2),:) + P(Fc(:,3),:)) / 3), [], 2) >= 0.02;
    end
    k = find(ok);  v = k(randi(numel(k)));
end

function Sg = i_sg(C, name, h)
    Sg = C.Sh;  Sg.nV = C.nV;  Sg.Hemi = {(1:C.nV)'};  Sg.HemiLabel = {char("Cortex " + h)};
    Sg.Comment = char(name);  if ~isfield(Sg, 'SurfaceFile'), Sg.SurfaceFile = ''; end
end

function d = i_dirac(name)
    try, d = rheome.load.dirac(name);
    catch, rheome.import.dirac(name);  d = rheome.load.dirac(name);
    end
end

function v = i_pick(c)
    v = c(randi(numel(c)));
end

function s = i_snr(s)
    if isinf(s), s = "inf"; else, s = string(s) + "dB"; end
end

% Author: Diellor Basha, 2026
