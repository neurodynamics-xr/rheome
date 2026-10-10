function [T, X, strips] = measure_catalogue(name, S, item, opts)
% SCALE.MEASURE_CATALOGUE  The catalogue's patterns planted through this participant's own head, and three estimators reading them (MS1 G6, G3 rule 11, G16).
%
%   [T, X, strips] = rheome.scale.measure_catalogue(name, S, "catalogue")     % G6 + G3 rule 11 + G16 on the plants
%   [T, X]         = rheome.scale.measure_catalogue(name, S, "catalognulls")  % G16 false-propagation nulls
%
% Port of the single-subject catalogue validation (MS1 Fig. 15A-B: one reference gain, source noise, one hemisphere) to
% every participant, both hemispheres (Hemis), the participant's own leadfield and whitened minimum norm
% (S.Res: rheome.inverse.mne on the participant's noise covariance), and SENSOR noise from the
% participant's own resting record and empty room (spec: the group-analysis plan,
% G6, G3 rule 11, G16).
%
% catalogue     the 19 catalogue plants (i_catalogue: planar, source, sink, rotor, spiral, pair, saddle,
%               standing, nested, each meso and macro, and a drifting rotor) x Centres random centres per
%               hemisphere (the same centres for every pattern, so patterns are comparable), each read in arms
%                 direct       the plant itself, no instrument (the "same plants without the instrument" baseline,
%                              and the truth every other arm is judged against: truthOnMesh)
%                 meg/inf      normal current -> own gain -> whitened MNE, read along normals smoothed over
%                              NormalSigmaMM (20 mm), noise-free
%                 meg/rest     the same + a window of the participant's own resting record at SNRdB
%                 meg/emptyroom the same + a window of the participant's own empty room at SNRdB (channels by
%                              name; NaN rows when the participant has none)
%                 meglocal/inf planar_macro only (G3 rule 11): the meg/inf arm read along the LOCAL normals
%               Estimators (column estimator):
%                 framework    flow.phasegradient on the band-limited analytic field (speed, wavelength, direction,
%                              rotation period), detect.criticalPoints on the amplitude-weighted phase-velocity
%                              direction, detect.phasesingularity (chirality, drift track), differential.helmholtz,
%                              the standing index, the wavelet-scale wavelength, PAC and the fast-envelope speed
%                 bst_of       Brainstorm's bst_opticalflow (rheome.flow.bstopticalflow) on |the band-limited
%                              field| over one cycle, at each HornSchunck (the first is Brainstorm's default 0.01;
%                              column param); arms direct, meg/inf, meg/rest
%                 phasereg     rheome.flow.phaseregression (circular-linear regression of phase on geodesic distance,
%                              PatchMM patches) at every vertex within 1.3 E of the centre; arms direct, meg/inf, meg/rest
%               The comparators are judged on the same quantities: phase speed, direction error (against the
%               direct arm's time-mean phase velocity), critical-point type, sign of div/curl at the core
%               (chirality for the rotational classes), and "propagation declared".
% catalognulls  G16 false-propagation nulls, NullsPerHemi per hemisphere per type, two Gaussian patches of normal
%               current (NullExtentMM) at 10 Hz, NullSepMM apart (geodesic):
%                 coherent     zero lag (stationary coherent dipoles)
%                 lagged       NullLag apart (two phase-lagged generators, Zhigalov & Jensen 2023)
%               arms direct, meg/inf, meg/rest; all three estimators. Nothing propagates, so "propagation declared"
%               is a false propagation.
%
% ⭐ PROPAGATION DECLARED, ONE RULE FOR EVERY ESTIMATOR (pre-declared): the amplitude-weighted directional
%   coherence of the estimator's time-mean velocity field over the region, |sum a v^| / sum a, is at least
%   PropCoherence (0.5) AND its median speed lies in PropSpeed (0.1-10 m/s).
% ⭐ SNR IS MEASURED IN THE READOUT BAND (pattern frequency x [0.7 1.3]): the signal and the noise window are
%   both band-passed there and their channel-summed powers set the noise gain; the noise is added broadband.
% ⭐ RATIOS: column ratio = recovered / ref, ref the catalogue's NOMINAL value where it has one (speed,
%   wavelength, period, drift, envelope speed), otherwise truthOnMesh (the direct arm's own reading).
% ⚠ THE CHART IS NOT THE SINGLE-SUBJECT VALIDATION'S. That took the azimuth about the centre on the reference subject's FreeSurfer registration sphere;
%   the nsp staging carries no sphere, so here theta is the azimuth in the centre's tangent plane on the cortex
%   SMOOTHED over NormalSigmaMM (LBO heat kernel on the coordinates); r is the geodesic distance on the folded
%   cortex, as there. The drifting rotor's core moves along the geodesic from the centre towards e1, not a great
%   circle. Every arm is judged against the direct arm on the same chart, so the chart cancels.
% ⚠ The sampling rate is the record's (600 Hz) decimated by an integer to >= SPC samples per cycle (nested: 10 per
%   fast cycle), so sensor-noise windows are the participant's own samples, resampled (anti-aliased) by that integer.
%
% X (one row per hemi x centre x pattern x arm x noise x estimator x param x quantity): nominal, truthOnMesh,
% recovered, ratio. T = rheome.scale.rows (analysis = item): per pattern x arm/noise x estimator x quantity the
% median ratio (metric ratio_median) or the fraction correct (metric preserved_frac); for the nulls the
% false-propagation rate per estimator x type x arm. strips: centre 1 of the first hemisphere, six phases of the
% middle cycle per pattern and arm (Graphics: Fig. 15 strips for the illustrating participant).
%
% See also: rheome.scale.run, rheome.flow.phaseregression, rheome.flow.bstopticalflow, rheome.flow.phasegradient,
%           rheome.detect.criticalPoints, rheome.scale.measure_plantfloors
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S
        item (1,1) string {mustBeMember(item, ["catalogue" "catalognulls"])}
        opts.Hemis string = ["L" "R"]
        opts.Centres (1,1) double {mustBeInteger, mustBePositive} = 3
        opts.Patterns string = string.empty
        opts.SNRdB (1,1) double = 10
        opts.NormalSigmaMM (1,1) double = 20
        opts.SPC (1,1) double = 16
        opts.HornSchunck double = [0.01 0.001 0.1 1]
        opts.PatchMM (1,1) double = 20
        opts.NullsPerHemi (1,1) double {mustBeInteger, mustBePositive} = 3
        opts.NullSepMM (1,2) double = [40 80]
        opts.NullExtentMM (1,1) double = 10
        opts.NullLag (1,1) double = pi/2
        opts.PropSpeed (1,2) double = [0.1 10]
        opts.PropCoherence (1,1) double = 0.5
        opts.Seed (1,1) double = 43
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    rng(opts.Seed);
    N = i_noise(name, S);
    strips = struct();
    if item == "catalogue", [X, strips] = i_catalogue_run(S, N, opts);  T = i_catalogue_rows(X);
    else, X = i_nulls_run(S, N, opts);  T = i_nulls_rows(X); end
end

%% ---------- the catalogue
function P = i_catalogue()
% name class scale f(Hz) lambda(mm) E(mm) sep(mm) fast(Hz) mod drift(m/s) -- the validation's table, verbatim
    P = cell2struct({
     'planar_meso'    'planar'   'meso'  20  16.5 25 NaN NaN NaN 0
     'planar_macro'   'planar'   'macro' 10  200  100 NaN NaN NaN 0
     'source_meso'    'source'   'meso'  10  57   30 NaN NaN NaN 0
     'source_macro'   'source'   'macro' 10  200  100 NaN NaN NaN 0
     'sink_meso'      'sink'     'meso'  10  57   30 NaN NaN NaN 0
     'sink_macro'     'sink'     'macro' 10  200  100 NaN NaN NaN 0
     'rotor_meso'     'rotor'    'meso'  10  NaN  15 NaN NaN NaN 0
     'rotor_macro'    'rotor'    'macro' 13  NaN  50 NaN NaN NaN 0
     'spiral_meso'    'spiral'   'meso'  12  20   15 NaN NaN NaN 0
     'spiral_macro'   'spiral'   'macro' 12  200  50 NaN NaN NaN 0
     'pair_meso'      'pair'     'meso'  12  NaN  20 15  NaN NaN 0
     'pair_macro'     'pair'     'macro' 12  NaN  60 60  NaN NaN 0
     'saddle_meso'    'saddle'   'meso'  2   20   20 NaN NaN NaN 0
     'saddle_macro'   'saddle'   'macro' 2   100  60 NaN NaN NaN 0
     'standing_meso'  'standing' 'meso'  10  20   25 NaN NaN NaN 0
     'standing_macro' 'standing' 'macro' 10  80   100 NaN NaN NaN 0
     'nested_meso'    'nested'   'meso'  5   30   25 NaN 40  0.5 0
     'nested_macro'   'nested'   'macro' 5   200  100 NaN 40  0.5 0
     'rotor_drift'    'drift'    'macro' 10  NaN  50 NaN NaN NaN 0.05
    }, {'name','class','scale','f','lam','E','sep','ff','mod','drift'}, 2);
    for i = 1:numel(P), P(i).name = string(P(i).name); P(i).class = string(P(i).class); P(i).scale = string(P(i).scale); end
end

function [X, strips] = i_catalogue_run(S, N, o)
    P = i_catalogue();  if ~isempty(o.Patterns), P = P(ismember([P.name], o.Patterns)); end
    X = table();  strips = struct();
    for h = o.Hemis
        C = i_hemi(S, h, o);
        cIdx = randperm(C.nV, o.Centres);  bAx = randn(o.Centres, 3);
        for r = 1:o.Centres
            for p = P'
                [fs, dec, t] = i_time(p, N.fs, o);
                [X0, tgt, path] = i_plant(p, C, cIdx(r), bAx(r, :), t);
                reg = i_region(C, tgt.r, p.E*1e-3);
                Q0 = i_readout(X0, C, fs, p, reg, tgt, path, []);
                Y0 = C.Gn * X0;
                arms = {"direct" "none" X0; "meg" "inf" C.Kn*Y0; ...
                        "meg" "rest" C.Kn*i_addnoise(Y0, N.rest, dec, fs, p, o.SNRdB); ...
                        "meg" "emptyroom" C.Kn*i_addnoise(Y0, N.er, dec, fs, p, o.SNRdB)};
                if p.name == "planar_macro", arms(end+1, :) = {"meglocal" "inf" C.Knl*Y0}; end %#ok<AGROW>
                for a = 1:size(arms, 1)
                    Xa = arms{a, 3};
                    if isempty(Xa), continue, end                     % no empty room
                    if a == 1, Q = Q0; else, Q = i_readout(Xa, C, fs, p, reg, tgt, path, Q0); end
                    key = {h, r, p, arms{a, 1}, arms{a, 2}};
                    X = [X; i_rows(key, "framework", NaN, i_frameworkq(p, Q, Q0))]; %#ok<AGROW>
                    if a <= 3                                          % comparators: direct, meg/inf, meg/rest
                        X = [X; i_comparators(key, Xa, C, fs, p, reg, tgt, Q0, o)]; %#ok<AGROW>
                    end
                    if r == 1 && h == o.Hemis(1) && a <= 3
                        strips.(p.name).(arms{a, 1} + "_" + arms{a, 2}) = i_strip(Xa, fs, p, Q);
                    end
                end
            end
            fprintf('[catalogue] %s centre %d/%d: %d patterns\n', h, r, o.Centres, numel(P));
        end
    end
end

function [fs, dec, t] = i_time(p, fs0, o)
    need = o.SPC * p.f;  nCyc = 4;
    if p.class == "nested", need = 10 * p.ff;  nCyc = 2; end
    if p.class == "drift", nCyc = 5; end
    dec = max(1, floor(fs0 / need));  fs = fs0 / dec;
    t = (0:round(nCyc * fs / p.f) - 1) / fs;
end

function [X, tgt, path] = i_plant(p, C, ic, bax, t)
% The validation's i_plant on the cortex chart (r geodesic on the folded surface, theta in the smoothed tangent plane)
    E = p.E*1e-3;  k = 2*pi/(p.lam*1e-3);  w2 = 2*pi*p.f;  path = [];
    [r, w] = i_chart(C, ic, bax);
    A = exp(-r.^2/(2*E^2));  core = @(rr) tanh(rr/(0.2*E));
    tgt = struct('r', r, 'pos', C.V(ic,:), 'ic', ic, 'centre', C.V(ic,:));
    switch p.class
        case "planar",   phi = k*real(w);
        case "source",   phi = k*r;
        case "sink",     phi = -k*r;
        case "rotor",    phi = angle(w);  A = A.*core(r);
        case "spiral",   phi = angle(w) + k*r;  A = A.*core(r);
        case "pair"
            d = p.sep*1e-3/2;  phi = angle((w-d).*conj(w+d));  A = A.*core(abs(w-d)).*core(abs(w+d));
            [~, i1] = min(abs(w-d));  [~, i2] = min(abs(w+d));  tgt.pos = C.V([i1 i2],:);
        case "saddle",   phi = pi/(p.lam*1e-3*E) * (real(w).^2 - imag(w).^2);
        case "standing", X = (A.*cos(k*real(w))) * cos(w2*t);  return
        case "nested"
            slow = cos(w2*t - k*real(w));  X = A.*(slow + p.mod*(1+slow)/2 .* cos(2*pi*p.ff*t));  return
        case "drift"                 % a rotor whose core moves along the geodesic towards e1 at p.drift m/s
            e1 = i_tangent(C, ic, bax);  L = p.drift * t(end);
            [~, vt] = min(vecnorm(C.Vsm - (C.Vsm(ic,:) + 1.2*L*e1), 2, 2));
            pv = shortestpath(C.ge, ic, vt);  arc = [0 cumsum(vecnorm(diff(C.V(pv,:)), 2, 2))'];
            X = zeros(C.nV, numel(t));  path = zeros(numel(t), 3);  cache = containers.Map('KeyType','double','ValueType','any');
            for j = 1:numel(t)
                iv = pv(find(arc <= p.drift*t(j), 1, 'last'));
                if ~isKey(cache, iv), [rj, wj] = i_chart(C, iv, bax);  cache(iv) = {rj, wj}; end
                rw = cache(iv);  rj = rw{1};  wj = rw{2};
                X(:,j) = exp(-rj.^2/(2*E^2)).*core(rj) .* cos(w2*t(j) - angle(wj));  path(j,:) = C.V(iv,:);
            end
            tgt.pos = path(round(numel(t)/2),:);  tgt.centre = tgt.pos;
            return
    end
    X = A .* cos(w2*t - phi);
end

function [r, w] = i_chart(C, ic, bax)
    r = distances(C.ge, ic)';  e1 = i_tangent(C, ic, bax);  e2 = cross(C.Ns(ic,:), e1);
    d = C.Vsm - C.Vsm(ic,:);  w = r .* exp(1i*atan2(d*e2', d*e1'));
end

function e1 = i_tangent(C, ic, bax)
    n = C.Ns(ic,:);  e1 = bax - (bax*n')*n;  e1 = e1/norm(e1);
end

function reg = i_region(C, r, E)
    F = C.F;  rf = mean(r(F), 2);
    reg.v = r <= E;  reg.f = all(reg.v(F), 2);  reg.rf = rf;  reg.r = r;  reg.E = E;
    reg.fOut = reg.f & rf > 0.2*E;  reg.ring = reg.f & rf >= 0.3*E & rf <= 0.8*E;  reg.core = r <= 0.3*E;
    reg.use = reg.v & r > 0.2*E;
end

function z = i_analytic(X, fs, band)                     % band-pass + analytic signal in one FFT mask
    nT = size(X,2);  f = (0:nT-1)*fs/nT;  m = 2*(f >= band(1) & f <= band(2));
    z = ifft(fft(X,[],2) .* m, [], 2);
end

function fr = i_frames(X, fs, p)
    spc = round(fs/p.f);  fr = spc + (1:min(spc+1, size(X,2)-spc));
    if p.class == "drift", fr = spc+1 : size(X,2)-spc; end      % ⚠ the FFT analytic wraps: drop one cycle each end
end

function Q = i_readout(X, C, fs, p, reg, tgt, path, Q0)
% The validation's i_readout (framework), less its optical flow (bst_of replaces it, in i_comparators)
    Sg = C.Sg;  unitv = @(V) V ./ max(vecnorm(V,2,2), eps);  ang = @(U,V) acosd(max(min(sum(unitv(U).*unitv(V),2),1),-1));
    z = i_analytic(X, fs, p.f*[0.7 1.3]);  fr = i_frames(X, fs, p);
    pg = rheome.flow.phasegradient(z(:,fr), Sg, 'Rate', fs);
    vbar = mean(pg.velocity, 3);  amp = mean(abs(z(:,fr)), 2);  ampF = mean(amp(C.F), 2);
    Q.vbar = vbar;  Q.speed = median(pg.speed(reg.fOut,:), 'all', 'omitnan');
    Q.wl = median(pg.wavelength(reg.fOut,:), 'all', 'omitnan');
    Q.period = 1e3*median(2*pi*reg.rf(reg.ring) ./ median(pg.speed(reg.ring,:), 2, 'omitnan'), 'omitnan');
    if isempty(Q0), Q.dirErr = 0; else, Q.dirErr = median(ang(vbar(reg.fOut,:), Q0.vbar(reg.fOut,:)), 'omitnan'); end
    zr = z(reg.v, fr);  a = C.av(reg.v);
    Q.stand = mean(abs(sum(a.*zr.^2, 1)) ./ sum(a.*abs(zr).^2, 1));
    Xr = real(z(:,fr(1:4:end))) .* reg.v;  En = mean(vertexSpectrum(C.gfb, Xr), 2);  ok = isfinite(C.wl(:));
    wl = C.wl(ok);  En = En(ok);  [~, j] = max(En);
    if j > 1 && j < numel(En), y = log(En(j-1:j+1)); d = (y(1)-y(3))/(2*(y(1)-2*y(2)+y(3)));
        Q.wlWav = 1e3*exp(log(wl(j)) + d*(log(wl(j+1))-log(wl(j))));
    else, Q.wlWav = 1e3*wl(j); end
    Jv = (C.Afv*vbar) ./ C.afv;  Jv(~isfinite(Jv)) = 0;
    Ju = (C.Afv*(ampF.*unitv(vbar))) ./ C.afv;  Ju(~isfinite(Ju)) = 0;
    Q.Jv = Jv;  Q.amp = amp;
    [Q.div, Q.curl] = i_core(Jv, C, reg);
    [Q.nCP, Q.cpType, Q.cpErr, Q.cpCharge] = i_cp(rheome.detect.criticalPoints(reshape(Ju.', [], 1), Sg, 'all', C.op), tgt, reg.E);
    Q.cpChance = 1e3*0.5*sqrt(pi*reg.E^2/max(Q.nCP,1));
    Q.coh = i_coherence(Jv, amp, reg.v);  Q.propOK = i_declared(Q.coh, Q.speed, C.o);
    fs4 = fr(1:4:end);  ns = zeros(numel(fs4), 2);  se = nan(numel(fs4), 1);  ch = nan(numel(fs4), 1);
    for j = 1:numel(fs4)
        ps = rheome.detect.phasesingularity(z(:,fs4(j)), Sg);
        in = vecnorm(ps.pos - tgt.centre, 2, 2) <= reg.E;
        ns(j,:) = [nnz(ps.charge(in) > 0) nnz(ps.charge(in) < 0)];
        if any(in)
            e = zeros(size(tgt.pos,1),1);  q = ps.charge(in);
            for ti = 1:size(tgt.pos,1), [e(ti), ii] = min(vecnorm(ps.pos(in,:) - tgt.pos(ti,:), 2, 2)); if ti == 1, ch(j) = q(ii); end, end
            se(j) = 1e3*max(e);
        end
    end
    Q.nSingPos = median(ns(:,1));  Q.nSingNeg = median(ns(:,2));
    Q.singChance = 1e3*0.5*sqrt(pi*reg.E^2/max(Q.nSingPos+Q.nSingNeg,1));  Q.singErr = median(se, 'omitnan');  Q.chir = mode(ch);
    if p.class == "drift"
        trk = nan(numel(fr), 3);  prev = path(1,:);
        for j = 1:numel(fr)
            ps = rheome.detect.phasesingularity(z(:,fr(j)), Sg);
            if isempty(ps.charge), continue; end
            [m, i] = min(vecnorm(ps.pos - prev, 2, 2));
            if m < 0.4*reg.E, trk(j,:) = ps.pos(i,:); prev = trk(j,:); end
        end
        ok = all(isfinite(trk), 2);  tt = (fr(:)-1)/fs;
        Q.pathErr = 1e3*median(vecnorm(trk(ok,:) - path(fr(ok),:), 2, 2));  Q.cover = mean(ok);
        Q.driftSpeed = NaN;  if nnz(ok) > 2, cf = [ones(nnz(ok),1) tt(ok)] \ trk(ok,:);  Q.driftSpeed = norm(cf(2,:)); end
        ct = [ones(numel(tt),1) tt] \ path(fr,:);  Q.driftTrue = norm(ct(2,:));
    end
    H = rheome.differential.helmholtz(reshape(Jv.', [], 1), Sg);
    es = sum(reshape(H.Vsol,3,[]).'.^2, 2);  ei = sum(reshape(H.Virr,3,[]).'.^2, 2);
    Q.solFrac = sum(C.av(reg.v).*es(reg.v)) / sum(C.av(reg.v).*(es(reg.v)+ei(reg.v)));
    pot = H.Phi;  if any(p.class == ["rotor" "spiral" "pair" "drift"]), pot = H.Psi; end
    rim = reg.v & reg.r >= 0.9*reg.E;  pr = pot(reg.v) - mean(pot(rim));  iv = find(reg.v);  [~, m] = max(abs(pr));
    Q.potErr = 1e3*min(vecnorm(C.V(iv(m),:) - tgt.pos, 2, 2));
    if p.class == "nested"
        zs = i_analytic(X, fs, p.f*[0.7 1.3]);  zf = i_analytic(X, fs, p.ff + p.f*[-1.5 1.5]);
        Af = abs(zf);  mvl = mean(Af.*exp(1i*angle(zs)), 2) ./ mean(Af, 2);
        Q.mi = median(abs(mvl(reg.v)));  Q.pacPhase = angle(sum(C.av(reg.v).*mvl(reg.v)));
        ze = i_analytic(Af - mean(Af,2), fs, p.f*[0.7 1.3]);
        pe = rheome.flow.phasegradient(ze(:,fr), Sg, 'Rate', fs);  Q.envSpeed = median(pe.speed(reg.fOut,:), 'all', 'omitnan');
    end
end

function [dv, cv] = i_core(Jv, C, reg)
    J3 = reshape(Jv.', [], 1);  w = C.av .* reg.core;
    dv = sum(w.*rheome.differential.divergence(J3, [], C.fg))/sum(w);  cv = sum(w.*rheome.differential.curl(J3, [], C.fg))/sum(w);
end

function c = i_coherence(Vv, amp, m)
% amplitude-weighted directional coherence |sum a v^| / sum a over the vertices m (0 = no common direction)
    u = Vv(m,:) ./ max(vecnorm(Vv(m,:), 2, 2), eps);  a = amp(m) .* (vecnorm(Vv(m,:), 2, 2) > 0);
    c = norm(sum(a.*u, 1)) / max(sum(a), eps);
end

function [n, ty, err, q] = i_cp(cp, tgt, E)     % critical points within E of the centre; the one nearest each target
    inR = vecnorm(cp.pos - tgt.centre, 2, 2) <= E;  n = nnz(inR);  ty = "none";  err = NaN;  q = NaN;
    for ti = 1:size(tgt.pos,1)
        if ~any(inR), break; end
        dd = vecnorm(cp.pos - tgt.pos(ti,:), 2, 2);  dd(~inR) = Inf;  [m, i] = min(dd);
        if ti == 1, ty = string(cp.type{i}); q = cp.charge(i); err = 1e3*m; else, err = max(err, 1e3*m); end
    end
end

function Wq = i_frameworkq(p, Q, Q0)
% the framework's quantities: {quantity, unit, nominal, truthOnMesh, recovered}
    c = p.lam*1e-3*p.f;  per = 1e3/p.f;  rot = any(p.class == ["rotor" "spiral" "pair" "drift"]);  cpc = ~any(p.class == ["standing" "planar" "nested"]);
    Q.cpTypeOK = double(Q.cpType == Q0.cpType && Q0.cpType ~= "none");  Q0.cpTypeOK = double(Q0.cpType ~= "none");
    Q.chirOK = double(Q.chir == Q0.chir);  Q0.chirOK = 1;
    if any(p.class == ["source" "sink" "saddle"]), Q.signOK = double(sign(Q.div) == sign(Q0.div));
    else, Q.signOK = double(sign(Q.curl) == sign(Q0.curl)); end
    Q0.signOK = 1;
    if p.class == "nested", Q.pacErr = abs(rad2deg(angle(exp(1i*(Q.pacPhase - Q0.pacPhase)))));  Q0.pacErr = 0; end
    if p.class == "drift", Q0.driftSpeed = Q0.driftTrue; end
    want = {
      'phase speed',                       'm/s', c,    'speed',   ~rot && p.class ~= "standing" && p.class ~= "saddle"
      'phase speed',                       'm/s', NaN,  'speed',   p.class == "spiral" || p.class == "saddle"
      'direction error',                   'deg', 0,    'dirErr',  p.class ~= "standing"
      'wavelength (phase gradient)',       'mm',  p.lam,'wl',      any(p.class == ["planar" "source" "sink" "nested"])
      'wavelength (wavelet scale)',        'mm',  NaN,  'wlWav',   any(p.class == ["planar" "source" "sink" "standing" "nested"])
      'rotation period',                   'ms',  per,  'period',  any(p.class == ["rotor" "drift"])
      'standing index',                    '1',   double(p.class=="standing"), 'stand', true
      'critical points in region',         'n',   NaN,  'nCP',     true
      'critical point type matches class', '1',   1,    'cpTypeOK',cpc
      'critical point index',              '1',   NaN,  'cpCharge',cpc
      'critical point location error',     'mm',  0,    'cpErr',   cpc
      'critical point location, chance',   'mm',  NaN,  'cpChance',cpc
      'phase singularities +1 in region',  'n',   NaN,  'nSingPos',true
      'phase singularities -1 in region',  'n',   NaN,  'nSingNeg',true
      'phase singularity location error',  'mm',  0,    'singErr', rot && p.class ~= "drift"
      'phase singularity location, chance','mm',  NaN,  'singChance', rot && p.class ~= "drift"
      'chirality matches plant',           '1',   1,    'chirOK',  rot
      'divergence at core',                '1/s', NaN,  'div',     any(p.class == ["source" "sink" "spiral" "saddle"])
      'curl at core',                      '1/s', NaN,  'curl',    any(p.class == ["rotor" "spiral" "drift"])
      'sign of div/curl matches plant',    '1',   1,    'signOK',  any(p.class == ["source" "sink" "spiral" "rotor" "drift" "saddle"])
      'Helmholtz extremum location error', 'mm',  0,    'potErr',  any(p.class == ["source" "sink" "rotor" "spiral" "drift"])
      'solenoidal fraction',               '1',   NaN,  'solFrac', p.class ~= "standing"
      'propagation declared',              '1',   NaN,  'propOK',  true
      'PAC modulation index',              '1',   NaN,  'mi',      p.class == "nested"
      'PAC preferred phase error',         'deg', 0,    'pacErr',  p.class == "nested"
      'fast-envelope speed',               'm/s', c,    'envSpeed',p.class == "nested"
      'drift speed',                       'm/s', p.drift,'driftSpeed', p.class == "drift"
      'core path error',                   'mm',  0,    'pathErr', p.class == "drift"
      'core tracked fraction',             '1',   1,    'cover',   p.class == "drift" };
    want = want([want{:,5}], :);
    Wq = cell(size(want,1), 5);
    for i = 1:size(want,1)
        f = want{i,4};  Wq(i,:) = {want{i,1}, want{i,2}, want{i,3}, double(Q0.(f)), double(Q.(f))};
    end
end

function y = i_declared(coh, sp, o)
    y = double(coh >= o.PropCoherence && sp >= o.PropSpeed(1) && sp <= o.PropSpeed(2));
end

function T = i_comparators(key, Xa, C, fs, p, reg, tgt, Q0, o)
% bst_opticalflow at each HornSchunck, and the phase regression, on the same band-limited field
    z = i_analytic(Xa, fs, p.f*[0.7 1.3]);  fr = i_frames(Xa, fs, p);  amp = mean(abs(z(:,fr)), 2);
    c = p.lam*1e-3*p.f;  rot = any(p.class == ["rotor" "spiral" "pair" "drift"]);  T = table();
    for hs = o.HornSchunck
        of = rheome.flow.bstopticalflow(real(z(:,fr)), C.Sg, fs, HornSchunck=hs);
        sp = median(vecnorm(of.velocity(reg.use,:,:), 2, 2), 'all', 'omitnan');
        T = [T; i_rows(key, "bst_of", hs, i_compq(p, mean(of.velocity, 3), sp, amp, C, reg, tgt, Q0, c, rot))]; %#ok<AGROW>
    end
    m = find(reg.r <= 1.3*reg.E);
    pr = rheome.flow.phaseregression(z(:,fr), C.Sg, Rate=fs, Centres=m, RadiusMM=o.PatchMM, Neighbours=C.nbr);
    Vp = zeros(C.nV, 3);  Vp(m,:) = pr.velocity;
    sp = median(pr.speed(reg.use(m)), 'omitnan');
    T = [T; i_rows(key, "phasereg", o.PatchMM, i_compq(p, Vp, sp, amp, C, reg, tgt, Q0, c, rot))];
end

function Wq = i_compq(p, Vv, sp, amp, C, reg, tgt, Q0, c, rot)
    unitv = @(V) V ./ max(vecnorm(V,2,2), eps);
    dirE = median(acosd(max(min(sum(unitv(Vv(reg.use,:)).*unitv(Q0.Jv(reg.use,:)),2),1),-1)), 'omitnan');
    [dv, cv] = i_core(Vv, C, reg);
    Ju = amp .* unitv(Vv);
    [~, ty] = i_cp(rheome.detect.criticalPoints(reshape(Ju.', [], 1), C.Sg, 'all', C.op), tgt, reg.E);
    if any(p.class == ["source" "sink" "saddle"]), sOK = double(sign(dv) == sign(Q0.div)); else, sOK = double(sign(cv) == sign(Q0.curl)); end
    cpc = ~any(p.class == ["standing" "planar" "nested"]);
    Wq = {'phase speed', 'm/s', i_tern(~rot && ~any(p.class == ["standing" "saddle"]), c, NaN), i_tern(p.class == "standing", NaN, Q0.speed), sp
          'direction error', 'deg', 0, 0, dirE
          'propagation declared', '1', NaN, Q0.propOK, i_declared(i_coherence(Vv, amp, reg.v), sp, C.o)};
    if cpc, Wq(end+1,:) = {'critical point type matches class', '1', 1, double(Q0.cpType ~= "none"), double(ty == Q0.cpType && Q0.cpType ~= "none")}; end
    if any(p.class == ["source" "sink" "spiral" "rotor" "drift" "saddle"]), Wq(end+1,:) = {'sign of div/curl matches plant', '1', 1, 1, sOK}; end
    if p.class == "standing", Wq = Wq([1 3], :); end
end

function T = i_rows(key, est, param, Wq)
    n = size(Wq, 1);  p = key{3};  rec = cell2mat(Wq(:,5));  nom = cell2mat(Wq(:,3));  tom = cell2mat(Wq(:,4));
    ref = nom;  ref(isnan(ref)) = tom(isnan(ref));  ratio = rec ./ ref;
    ratio(~ismember(string(Wq(:,1)), ["phase speed" "wavelength (phase gradient)" "wavelength (wavelet scale)" ...
          "rotation period" "drift speed" "fast-envelope speed"])) = NaN;
    T = table(repmat(key{1},n,1), repmat(key{2},n,1), repmat(p.name,n,1), repmat(p.class,n,1), repmat(p.scale,n,1), ...
        repmat(key{4},n,1), repmat(key{5},n,1), repmat(est,n,1), repmat(param,n,1), string(Wq(:,1)), string(Wq(:,2)), ...
        nom, tom, rec, ratio, 'VariableNames', {'hemi','centre','pattern','class','scale','arm','noise','estimator', ...
        'param','quantity','unit','nominal','truthOnMesh','recovered','ratio'});
end

function T = i_catalogue_rows(X)
    T = table();  if isempty(X), return, end
    bin = ["critical point type matches class" "chirality matches plant" "sign of div/curl matches plant" "propagation declared"];
    G = X(:, {'pattern','arm','noise','estimator','quantity'});  G.param = compose("%g", X.param);  [g, K] = findgroups(G);   % ⚠ compose: string(NaN) is <missing>, which findgroups drops
    rat = splitapply(@(v) median(v, 'omitnan'), X.ratio, g);  frac = splitapply(@(v) mean(v, 'omitnan'), X.recovered, g);
    for i = 1:numel(rat)
        b = K.pattern(i) + "_" + K.arm(i) + "_" + K.noise(i) + "_" + K.estimator(i) + i_par(str2double(K.param(i))) + "_" + i_slug(K.quantity(i));
        if ismember(K.quantity(i), bin), T = [T; rheome.scale.rows("catalogue", "preserved_frac", frac(i), "fraction", b)]; %#ok<AGROW>
        elseif isfinite(rat(i)), T = [T; rheome.scale.rows("catalogue", "ratio_median", rat(i), "ratio", b)]; %#ok<AGROW>
        end
    end
end

%% ---------- the false-propagation nulls
function X = i_nulls_run(S, N, o)
    X = table();  fs0 = N.fs;  p = struct('name', "null", 'class', "null", 'scale', "", 'f', 10, 'lam', NaN, 'E', NaN);
    dec = max(1, floor(fs0 / (o.SPC*p.f)));  fs = fs0/dec;  t = (0:round(4*fs/p.f)-1)/fs;  E = o.NullExtentMM*1e-3;
    for h = o.Hemis
        C = i_hemi(S, h, o);
        for ty = ["coherent" "lagged"]
            lag = 0;  if ty == "lagged", lag = o.NullLag; end
            for n = 1:o.NullsPerHemi
                v1 = randi(C.nV);  d1 = distances(C.ge, v1)';  cand = find(d1 >= o.NullSepMM(1)*1e-3 & d1 <= o.NullSepMM(2)*1e-3);
                if isempty(cand), [~, k] = min(abs(d1 - mean(o.NullSepMM)*1e-3)); cand = k; end
                v2 = cand(randi(numel(cand)));  d2 = distances(C.ge, v2)';  pv = shortestpath(C.ge, v1, v2);
                mid = pv(ceil(numel(pv)/2));  r = distances(C.ge, mid)';
                X0 = exp(-d1.^2/(2*E^2)) * cos(2*pi*p.f*t) + exp(-d2.^2/(2*E^2)) * cos(2*pi*p.f*t - lag);
                reg = i_region(C, r, d1(v2)/2 + 2*E);  reg.use = reg.v;
                Y0 = C.Gn * X0;
                arms = {"direct" "none" X0; "meg" "inf" C.Kn*Y0; "meg" "rest" C.Kn*i_addnoise(Y0, N.rest, dec, fs, p, o.SNRdB)};
                for a = 1:3
                    z = i_analytic(arms{a,3}, fs, p.f*[0.7 1.3]);  fr = i_frames(arms{a,3}, fs, p);  amp = mean(abs(z(:,fr)), 2);
                    pg = rheome.flow.phasegradient(z(:,fr), C.Sg, 'Rate', fs);
                    Jv = (C.Afv*mean(pg.velocity, 3)) ./ C.afv;  Jv(~isfinite(Jv)) = 0;
                    sp = median(pg.speed(reg.f,:), 'all', 'omitnan');
                    res = {"framework", NaN, i_coherence(Jv, amp, reg.v), sp};
                    for hs = o.HornSchunck
                        of = rheome.flow.bstopticalflow(real(z(:,fr)), C.Sg, fs, HornSchunck=hs);
                        res(end+1,:) = {"bst_of", hs, i_coherence(mean(of.velocity, 3), amp, reg.v), ...
                                        median(vecnorm(of.velocity(reg.v,:,:), 2, 2), 'all', 'omitnan')}; %#ok<AGROW>
                    end
                    m = find(reg.v);
                    pr = rheome.flow.phaseregression(z(:,fr), C.Sg, Rate=fs, Centres=m, RadiusMM=o.PatchMM, Neighbours=C.nbr);
                    Vp = zeros(C.nV, 3);  Vp(m,:) = pr.velocity;
                    res(end+1,:) = {"phasereg", o.PatchMM, i_coherence(Vp, amp, reg.v), median(pr.speed, 'omitnan')}; %#ok<AGROW>
                    for k = 1:size(res, 1)
                        X = [X; table(h, n, ty, 1e3*d1(v2), arms{a,1}, arms{a,2}, res{k,1}, res{k,2}, res{k,3}, res{k,4}, ...
                             i_declared(res{k,3}, res{k,4}, o), 'VariableNames', {'hemi','null','type','sepMM','arm','noise', ...
                             'estimator','param','coherence','speed','declared'})]; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('[catalognulls] %s: %d nulls x 2 types\n', h, o.NullsPerHemi);
    end
end

function T = i_nulls_rows(X)
    T = table();
    G = X(:, {'type','arm','noise','estimator'});  G.param = compose("%g", X.param);  [g, K] = findgroups(G);   % ⚠ compose: string(NaN) is <missing>, which findgroups drops
    r = splitapply(@mean, X.declared, g);
    for i = 1:numel(r)
        b = K.type(i) + "_" + K.arm(i) + "_" + K.noise(i) + "_" + K.estimator(i) + i_par(str2double(K.param(i)));
        T = [T; rheome.scale.rows("catalognulls", "false_prop_rate", r(i), "fraction", b)]; %#ok<AGROW>
    end
end

%% ---------- the instrument, the noise, the hemisphere
function C = i_hemi(S, h, o)
    H = S.B.(char(h));  Sh = H.S;  C.V = double(Sh.Vertices);  C.F = double(Sh.Faces);  C.nV = size(C.V, 1);
    unit = @(X) X ./ max(vecnorm(X,2,2), eps);
    lbo = H.lbo;  lam = lbo.Lambda(:);  Phi = lbo.Phi;  M = lbo.Mass;  hk = exp(-lam*(o.NormalSigmaMM*1e-3)^2/2);
    Nl = unit(Sh.VertNormals);  C.Ns = unit(Phi*(hk .* (Phi'*(M*Nl))));  C.Vsm = Phi*(hk .* (Phi'*(M*C.V)));
    gv = double(H.gv(:));  K = S.Res.ImagingKernel;
    C.Kn  = C.Ns(:,1).*K(3*gv-2,:) + C.Ns(:,2).*K(3*gv-1,:) + C.Ns(:,3).*K(3*gv,:);     % read along smoothed normals
    C.Knl = Nl(:,1).*K(3*gv-2,:) + Nl(:,2).*K(3*gv-1,:) + Nl(:,3).*K(3*gv,:);            % G3 rule 11: local normals
    C.Gn  = S.G(:,3*gv-2).*Nl(:,1)' + S.G(:,3*gv-1).*Nl(:,2)' + S.G(:,3*gv).*Nl(:,3)';     % the plant is a normal current
    C.Sg = struct('Vertices', C.V, 'Faces', C.F, 'VertNormals', Nl, 'nV', C.nV, 'Hemi', {{(1:C.nV)'}}, ...
                  'HemiLabel', {{char("Cortex " + h)}}, 'Comment', '', 'SurfaceFile', '');
    C.fg = rheome.operators.face_gradient(C.V, C.F);  C.op = rheome.detect.operator(C.Sg);  C.ge = rheome.geom.edgegraph(Sh);
    C.av = full(sum(M, 2));  nF = size(C.F, 1);
    C.Afv = sparse(C.F(:), repmat((1:nF)',3,1), repmat(C.fg.FaceArea,3,1), C.nV, nF);  C.afv = full(sum(C.Afv, 2));
    C.gfb = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'VoicesPerOctave',8, 'SizeLimits',[4 80]*1e-3, ...
        'Transform', rheome.graphtransform.eigen(Phi, M, lam));
    C.wl = wavelengths(C.gfb);  C.o = o;
    pr = rheome.flow.phaseregression(zeros(C.nV, 2), C.Sg, Rate=1, RadiusMM=o.PatchMM);   % every patch, once per hemisphere
    C.nbr = pr.Neighbours;
end

function N = i_noise(name, S)
% the participant's own resting record and empty room on S's channels (empty room by name; absent channels 0)
    st = rheome.load.study(name);  N.fs = st.rec.sfreq;  N.rest = double(st.rec.F(S.iSel, :));
    nm = string(st.chan.Name(S.iSel));  clear st
    N.er = [];  f = fullfile(rheome.load.root(), name, 'noise.mat');
    if ~isfile(f), return, end
    nrec = getfield(builtin('load', f, 'nrec'), 'nrec');
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;  jo = find(ok);
    [~, iE, jn] = intersect(nm, nn(ok), 'stable');
    if isempty(iE), return, end
    assert(abs(nrec.sfreq - N.fs) < 1e-6, 'scale:catalogue:rate', 'empty room at %g Hz, record at %g Hz', nrec.sfreq, N.fs);
    N.er = zeros(numel(nm), size(nrec.F, 2));  N.er(iE, :) = double(nrec.F(jo(jn), :));
end

function y = i_addnoise(y0, B, dec, fs, p, snr)
% a random window of B (at the record's rate), anti-aliased down by dec, scaled to snr dB in the readout band
    if isempty(B), y = []; return, end
    nT = size(y0, 2);  L = nT*dec;  s0 = randi([1, size(B, 2) - L]);  nz = B(:, s0 + (0:L-1));
    if dec > 1, nz = resample(nz.', 1, dec).';  end
    nz = nz(:, 1:nT);  band = p.f*[0.7 1.3];
    ps = sum(real(i_analytic(y0, fs, band)).^2, 'all');  pn = sum(real(i_analytic(nz - mean(nz, 2), fs, band)).^2, 'all');
    y = y0 + sqrt(ps / (pn * 10^(snr/10))) * nz;
end

%% ---------- small things
function s = i_strip(X, fs, p, Q)        % six phases of the middle cycle + the recovered field, for Graphics
    spc = round(fs/p.f);  fr = min(spc + round(linspace(1, spc, 6)), size(X, 2));
    s = struct('frames', single(X(:,fr)), 'timesS', (fr-1)/fs, 'vertexVelocity', single(Q.Jv));
end

function s = i_par(v)
    if isnan(v), s = ""; else, s = "_" + string(v); end
end

function s = i_slug(q)
    s = regexprep(lower(q), '[^a-z0-9]+', '_');  s = regexprep(s, '^_|_$', '');
end

function y = i_tern(c, a, b), if c, y = a; else, y = b; end, end

% Author: Diellor Basha, 2026
