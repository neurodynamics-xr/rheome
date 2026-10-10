function [T, X] = measure_plantfloors(name, S, opts)
% SCALE.MEASURE_PLANTFLOORS  Helmholtz-band resolution floors through this participant's MEG (MS1 G1, G9).
%
%   [T, X] = rheome.scale.measure_plantfloors(name)
%   [T, X] = rheome.scale.measure_plantfloors(name, S, LambdaMM=[46 65 92 130 183], PlantsPerCell=10)
%
% The group version of the Helmholtz-noise planting test (MS1 Table 3, Fig. 4C-E), on the participant's own
% cortex, leadfield, plain minimum norm (S.Res) and 8-16 Hz resting background:
%   field    a source (tangential gradient of a band atom) at one vertex plus a vortex (n x gradient) at
%            another 60-120 mm away, each at unit energy, oscillating at F0 Hz with a random phase
%   sensors  its forward field plus a WindowS-second window of the participant's own rest, 8-16 Hz
%   SNR      signal power / background power over the window, summed over channels (dB)
%   readout  lock-in on the carrier -> complex sensor pattern -> MNE -> its real part at the pattern's
%            own phase -> rheome.differential.helmholtzbands; the band map's peak is the located vertex
% ⭐ BALANCED CELLS. Per hemisphere, PlantsPerCell fields per planted band: the source bands are
% repelem(bands, PlantsPerCell) and the vortex bands a random permutation of them, so every
% (type, band, hemisphere) cell holds exactly PlantsPerCell plants (the single-subject test drew bands at random: n 13-19).
% ⭐ SPLIT HALVES. Half the plants of every cell (source and vortex alike) take their background from
% the first half of the recording, the rest from the second; column `half` carries it for the ICC.
%
% Per plant and condition (X, one row per plant x type x condition):
%   errMM        geodesic distance from the planted vertex to the peak of the band map at the KNOWN band
%   blindErrMM   the same at the band of maximum energy of that part (rule 2), blindShift its band offset
%   chanceMM     to a uniformly random vertex of the same hemisphere (the chance baseline)
%   normalFrac   energy share of the normal part of the read field (G9: spurious normal energy)
%   condition    "meg" (through leadfield + MNE, at snrDB) or "direct" (the planted field itself, no
%                instrument; snrDB Inf)
% Returns rheome.scale.rows (analysis "plantfloors"), pooled over hemispheres:
%   err_mm        median errMM per type x band x SNR             band "<type>_<lambda>mm_<snr>"
%   direct_err_mm, chance_mm   median per type x band              band "<type>_<lambda>mm"
%   floor_mm      the finest planted band whose median err_mm <= FloorMM (33 mm = one depth-7 tile;
%                 NaN when none), per type x SNR                   band "<type>_<snr>"
%   blind_shift   median band offset of the blind choice, noise-free, per type
%   mixed_err_mm  noise-free median over all bands, per type (MS1 Fig. 4C: vortex 20, source 80)
%   normal_frac   noise-free median spurious normal share through MEG (single-subject reference: 0.27)
% plus the same floors per half (band "<type>_<snr>_h1|h2") for the split-half reliability.
%
% See also: rheome.differential.helmholtzbands, rheome.scale.run, rheome.scale.measure_movingvortex
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.LambdaMM double = [46 65 92 130 183]
        opts.SNRdB double = [Inf 20 10 5 0 -5 -10]
        opts.PlantsPerCell (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.Hemis string = ["L" "R"]
        opts.F0 (1,1) double = 10
        opts.WindowS (1,1) double = 1
        opts.FloorMM (1,1) double = 33
        opts.Seed (1,1) double = 29
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    Nb = filtfilt(bp, double(st.rec.F(S.iSel, :)).').';  clear st
    nW = round(opts.WindowS * fs);  t = (0:nW-1) / fs;  car = exp(-1i * 2*pi*opts.F0*t);
    nS = size(Nb, 2);  edges = [fs, floor(nS/2); floor(nS/2) + 1, nS - nW - fs];   % background span per half
    rng(opts.Seed);  nL = numel(opts.LambdaMM);  ppc = opts.PlantsPerCell;  nSNR = numel(opts.SNRdB);

    X = table();
    for h = opts.Hemis
        H = S.B.(char(h));  Sh = H.S;  lbo = H.lbo;  nV = size(Sh.Vertices, 1);
        rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
        Kinv = S.Res.ImagingKernel(rows, :);  GL = S.G(:, rows);
        lam = lbo.Lambda(:);  P = lbo.Phi;  Mm = lbo.Mass;
        nrm = Sh.VertNormals ./ max(vecnorm(Sh.VertNormals, 2, 2), eps);
        fg = rheome.operators.face_gradient(Sh.Vertices, double(Sh.Faces));  ge = rheome.geom.edgegraph(Sh);
        gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2);
        G = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0));
        wl = 1e3 * wavelengths(gfb);
        [~, mb] = min(abs(wl(:) - opts.LambdaMM(:)'), [], 1);          % planted band per requested lambda
        atom = @(v, g) P * (g .* (P' * (Mm * full(sparse(v, 1, 1, nV, 1)))));
        bs = repelem(1:nL, ppc);  hfs = 1 + (mod(0:numel(bs)-1, ppc) >= ppc/2);  bv = bs;
        for q = 1:2, k = find(hfs == q);  bv(k) = bs(k(randperm(numel(k)))); end   % vortex cells balanced per half too
        for n = 1:numel(bs)
            v1 = randi(nV);  d1 = distances(ge, v1)';  cand = find(d1 > 0.06 & d1 < 0.12);
            if isempty(cand), [~, o] = sort(abs(d1 - 0.09));  cand = o(1:10); end
            v2 = cand(randi(numel(cand)));  d2 = distances(ge, v2)';
            ms = mb(bs(n));  mv = mb(bv(n));
            gs = i_tgrad(atom(v1, G(:, ms)), fg, nrm);  gw = i_tgrad(atom(v2, G(:, mv)), fg, nrm);
            J0 = reshape(gs', [], 1) / norm(gs(:)) + reshape(cross(nrm, gw, 2)', [], 1) / norm(gw(:));
            sig = (GL * J0) * cos(2*pi*opts.F0*t + 2*pi*rand);
            hf = hfs(n);  s0 = randi(edges(hf, :));
            noise = Nb(:, s0 + (0:nW-1));  ps = sum(sig(:).^2) / nW;  pn = sum(noise(:).^2) / nW;
            Jr = zeros(numel(J0), nSNR + 1);
            for k = 1:nSNR
                if isinf(opts.SNRdB(k)), y = sig; else, y = sqrt(pn / ps * 10^(opts.SNRdB(k)/10)) * sig + noise; end
                Jc = Kinv * (y * car.' * 2 / nW);  Jr(:, k) = real(Jc * exp(-1i * 0.5 * angle(sum(Jc.^2))));
            end
            Jr(:, end) = J0;                                                 % no instrument
            Bd = rheome.differential.helmholtzbands(Jr, Sh, lbo, Maps=true);
            assert(isequal(size(Bd.G), size(G)) && max(abs(Bd.G(:) - G(:))) < 1e-9, ...
                'scale:plantfloors:bank', 'helmholtzbands built a different bank from the planting bank');
            ch = [d1(randi(nV)) d2(randi(nV))];
            for k = 1:nSNR + 1
                [~, mi] = max(Bd.Eirr(:, k));  [~, mo] = max(Bd.Esol(:, k));
                [~, pS] = max(abs(Bd.PhiBand(:, ms, k)));  [~, pV] = max(abs(Bd.PsiBand(:, mv, k)));
                [~, pSb] = max(abs(Bd.PhiBand(:, mi, k)));  [~, pVb] = max(abs(Bd.PsiBand(:, mo, k)));
                if k <= nSNR, cond = "meg"; snr = opts.SNRdB(k); else, cond = "direct"; snr = Inf; end
                X = [X; table([h; h], [n; n], [hf; hf], ["source"; "vortex"], opts.LambdaMM([bs(n); bv(n)])', ...
                    reshape(wl([ms mv]), [], 1), [cond; cond], [snr; snr], 1e3*[d1(pS); d2(pV)], 1e3*[d1(pSb); d2(pVb)], ...
                    [mi - ms; mo - mv], 1e3*ch(:), Bd.normalFrac([k; k])', ...
                    'VariableNames', {'hemi','plant','half','type','lambdaMM','bandMM','condition','snrDB', ...
                    'errMM','blindErrMM','blindShift','chanceMM','normalFrac'})]; %#ok<AGROW>
            end
        end
        fprintf('[plantfloors %s] %s: %d mixed fields x %d SNR\n', name, h, numel(bs), nSNR);
    end
    T = i_metrics(X, opts);
end

function T = i_metrics(X, opts)
    T = table();  r = @(m, v, u, b) rheome.scale.rows("plantfloors", m, v, u, b);
    meg = X(X.condition == "meg", :);  dir = X(X.condition == "direct", :);
    for ty = ["source" "vortex"]
        for l = opts.LambdaMM
            k = dir.type == ty & dir.lambdaMM == l;  b = ty + "_" + l + "mm";
            T = [T; r(["direct_err_mm" "chance_mm"], [median(dir.errMM(k)) median(dir.chanceMM(k))], "mm", b)]; %#ok<AGROW>
            for s = opts.SNRdB
                T = [T; r("err_mm", median(meg.errMM(meg.type == ty & meg.lambdaMM == l & meg.snrDB == s)), "mm", ...
                          b + "_" + i_snr(s))]; %#ok<AGROW>
            end
        end
        for s = opts.SNRdB
            for hf = 0:2                                               % 0 = all plants, 1/2 = the halves
                k = meg.type == ty & meg.snrDB == s & (hf == 0 | meg.half == hf);
                fl = i_floor(meg(k, :), opts);  b = ty + "_" + i_snr(s);  if hf, b = b + "_h" + hf; end
                T = [T; r("floor_mm", fl, "mm", b)]; %#ok<AGROW>
            end
        end
        k = meg.type == ty & isinf(meg.snrDB);
        T = [T; r(["blind_shift" "mixed_err_mm"], [median(meg.blindShift(k)) median(meg.errMM(k))], ["bands" "mm"], ty)]; %#ok<AGROW>
    end
    T = [T; r("normal_frac", median(meg.normalFrac(isinf(meg.snrDB))), "fraction", "")];
end

function fl = i_floor(Y, opts)
% the finest planted band whose median error is within one depth-7 tile (pre-declared, Fig. 4 dashed line)
    fl = NaN;
    for l = sort(opts.LambdaMM)
        if median(Y.errMM(Y.lambdaMM == l)) <= opts.FloorMM, fl = l; return, end
    end
end

function s = i_snr(s)
    if isinf(s), s = "inf"; else, s = string(s) + "dB"; end
end

function g = i_tgrad(psi, fg, nrm)
    g = [fg.W*(fg.Gx*psi) fg.W*(fg.Gy*psi) fg.W*(fg.Gz*psi)];  g = g - sum(g.*nrm, 2) .* nrm;
end

% Author: Diellor Basha, 2026
