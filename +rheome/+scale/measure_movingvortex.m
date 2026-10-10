function [T, X] = measure_movingvortex(name, S, opts)
% SCALE.MEASURE_MOVINGVORTEX  A moving vortex planted through this participant's MEG and read back (MS1 G11).
%
%   [T, X] = rheome.scale.measure_movingvortex(name)
%   [T, X] = rheome.scale.measure_movingvortex(name, S, Placements=10, SigmaMM=[44 65 92 130], SNRdB=[Inf 10 0])
%
% The group version of MS1 Fig. 11 / section 7.2 (one placement on one participant):
%   plant    rheome.flow.movingvortex, SpeedMS 0.10 along a geodesic (126 mm at the default 140 mm atom),
%            alpha carrier in quadrature, 4 s; Placements start vertices per hemisphere, drawn at random
%            among vertices >= 20 mm from the gauge's singular faces (the default start is deterministic)
%   sensors  forwarded through the participant's leadfield, plus a 4 s stretch of the participant's own
%            rest at SNRdB (8-16 Hz power over the atom's temporal support); Inf = the planted field alone
%   readout  8-16 Hz Butterworth + Hilbert at the sensors -> the fused stream-function kernel
%            (rheome.differential.helmholtz applied to the hemisphere's MNE kernel, i.e. Psi of the
%            minimum-norm current) -> amplitude in one mexhat graph-wavelet member of width sigma
%            (bank 12-160 mm, 1 voice, exact eigenbasis); the core is the peak of that amplitude summed
%            over +-0.1 s around the frame
%   frames   every fs/20-th sample inside the support (12 at 600 Hz, as the single-participant figure)
% ⭐ The no-instrument arm reads the TRUE current's stream function with the same band and frames, so
% what the instrument costs is coreErr - directErr at the same sigma.
% ⚠ Core distances are Euclidean, as the single-participant figure; the path and speed are geodesic.
% ⚠ The direct arm's analytic signal takes the Hilbert transform of the quadrature time factors and
% holds the waypoint weights fixed, which is exact for a stationary atom and close for a slow one.
%
% Per placement x sigma x SNR (X): median/IQR/max core error (mm), the net displacement of the true and
% the recovered core and their ratio, the median chance error (a random vertex of the atom's support
% per frame), the direct (no-instrument) median error, the path length and realised speed.
% Placements 1..P/2 take their background from the first half of the recording (`half`).
% Returns rheome.scale.rows (analysis "movingvortex"), medians over placements and hemispheres:
%   core_err_mm, net_ratio per sigma x SNR (band "s<sigma>_<snr>"), chance_mm and direct_err_mm per
%   sigma (band "s<sigma>"), and core_err_mm per half (band "s<sigma>_<snr>_h1|h2").
%
% See also: rheome.flow.movingvortex, rheome.differential.helmholtz, rheome.scale.measure_plantfloors
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Placements (1,1) double {mustBeInteger, mustBePositive} = 10
        opts.SigmaMM double = [44 65 92 130]
        opts.SNRdB double = [Inf 10 0]
        opts.SpeedMS (1,1) double = 0.10
        opts.DurationS (1,1) double = 4
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 31
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    [bb, aa] = butter(3, [8 16] / (fs/2), 'bandpass');
    nS = size(F, 2);  rng(opts.Seed);  X = table();
    for h = opts.Hemis
        Hm = S.B.(char(h));  Sh = Hm.S;  lbo = Hm.lbo;  P = Sh.Vertices;  nV = size(P, 1);
        rows = reshape((double(Hm.gv(:))' - 1) * 3 + (1:3)', [], 1);
        GL = S.G(:, rows);
        Hk = rheome.differential.helmholtz(S.Res.ImagingKernel(rows, :), Sh);  Kpsi = Hk.Psi;  clear Hk   % [nV x nCh]
        gfb = rheome.graphfilterbank(lbo.Lambda(:), 'Wavelet', 'mexhat', 'VoicesPerOctave', 1, ...
              'SizeLimits', [12e-3 160e-3], 'Transform', rheome.graphtransform.eigen(lbo.Phi, lbo.Mass, lbo.Lambda(:)));
        sg = 1e3 * widths(gfb);  [~, mS] = min(abs(sg(:) - opts.SigmaMM(:)'), [], 1);
        band = @(m, Z) lbo.Phi * (feval(gain(gfb, m), lbo.Lambda(:)) .* (lbo.Phi' * (lbo.Mass * Z)));
        Kb = arrayfun(@(m) band(m, Kpsi), mS, 'uni', 0);  clear Kpsi              % one fused kernel per sigma
        g = rheome.operators.gauge(P, double(Sh.Faces), Method="diffusion");
        ok = true(nV, 1);
        if ~isempty(g.singular)
            Fc = double(Sh.Faces(g.singular, :));
            fc = (P(Fc(:,1),:) + P(Fc(:,2),:) + P(Fc(:,3),:)) / 3;  ok = min(pdist2(P, fc), [], 2) >= 0.02;
        end
        okv = find(ok);
        for p = 1:opts.Placements
            v0 = okv(randi(numel(okv)));
            mv = rheome.flow.movingvortex(name, SpeedMS=opts.SpeedMS, Duration=opts.DurationS, SampleRate=fs, ...
                                          Vertex=v0, Hemi=h, Bases=S.B, Gauge=g);
            nT = size(mv.W, 2);
            Bc = (GL * mv.JA) * (mv.W .* mv.a(:)') + (GL * mv.JB) * (mv.W .* mv.b(:)');
            envv = abs(mv.a(:) + 1i*mv.b(:));  sup = envv > 0.3 * max(envv);
            hf = 1 + (p > opts.Placements/2);
            lo = [fs, floor(nS/2)];  hi = [floor(nS/2) - nT, nS - nT - fs];  s0 = randi([lo(hf) hi(hf)]);
            Nz = F(:, s0 + (1:nT));  Nz = Nz - mean(Nz, 2);
            Pc = mean(mean(filtfilt(bb, aa, Bc.').'.^2, 1)' .* sup) / mean(sup);
            Pn = mean(mean(filtfilt(bb, aa, Nz.').'.^2, 1)' .* sup) / mean(sup);
            fr = find(sup);  fr = fr(1:round(fs/20):end);  fr = fr(:)';  nF = numel(fr);
            wp = P(mv.path(:), :);  tru = zeros(nF, 1);
            for i = 1:nF
                w = mv.W(:, fr(i));  [~, tru(i)] = min(vecnorm(P - (w(:)' * wp) / sum(w), 2, 2));
            end
            Jpk = max(abs(mv.JA), [], 2);  Jpk = max(reshape(Jpk, 3, []), [], 1)';
            supV = find(Jpk > 0.1 * max(Jpk));
            chance = median(vecnorm(P(supV(randi(numel(supV), nF, 1)), :) - P(tru, :), 2, 2)) * 1e3;
            ttrue = 1e3 * norm(P(tru(end), :) - P(tru(1), :));
            % no instrument: the true stream function, analytic through its quadrature time factors
            HA = rheome.differential.helmholtz(mv.JA, Sh);  HB = rheome.differential.helmholtz(mv.JB, Sh);  psiA = HA.Psi;  psiB = HB.Psi;
            ha = hilbert(mv.a(:)).';  hb = hilbert(mv.b(:)).';
            for si = 1:numel(mS)
                dirV = zeros(nF, 1);
                for i = 1:nF
                    u = max(1, fr(i)-60):10:min(nT, fr(i)+60);
                    Zt = band(mS(si), psiA * (mv.W(:, u) .* ha(u)) + psiB * (mv.W(:, u) .* hb(u)));
                    [~, dirV(i)] = max(sum(abs(Zt).^2, 2));
                end
                dirErr = 1e3 * vecnorm(P(dirV, :) - P(tru, :), 2, 2);
                for snr = opts.SNRdB
                    if isinf(snr), B = Bc; else, B = sqrt(Pn/Pc * 10^(snr/10)) * Bc + Nz; end
                    Zb = hilbert(filtfilt(bb, aa, B.')).';
                    rec = zeros(nF, 1);
                    for i = 1:nF
                        u = max(1, fr(i)-60):10:min(nT, fr(i)+60);
                        [~, rec(i)] = max(sum(abs(Kb{si} * Zb(:, u)).^2, 2));
                    end
                    err = 1e3 * vecnorm(P(rec, :) - P(tru, :), 2, 2);
                    netR = 1e3 * norm(P(rec(end), :) - P(rec(1), :));
                    X = [X; table(h, p, hf, v0, sg(mS(si)), opts.SigmaMM(si), snr, nF, median(err), ...
                         prctile(err, 25), prctile(err, 75), max(err), ttrue, netR, netR / max(ttrue, eps), ...
                         chance, median(dirErr), mv.pathLengthMM, mv.check.speedGeodesicMS, ...
                         'VariableNames', {'hemi','placement','half','startVertex','sigmaMM','sigmaReqMM','snrDB', ...
                         'nFrames','coreErrMM','coreErrQ1','coreErrQ3','coreErrMax','netTrueMM','netRecMM', ...
                         'netRatio','chanceMM','directErrMM','pathMM','speedMS'})]; %#ok<AGROW>
                end
            end
        end
        fprintf('[movingvortex %s] %s: %d placements x %d sigma x %d SNR\n', name, h, opts.Placements, ...
                numel(mS), numel(opts.SNRdB));
    end
    T = table();  r = @(m, v, u, b) rheome.scale.rows("movingvortex", m, v, u, b);
    for s = opts.SigmaMM
        k = X.sigmaReqMM == s;  b = "s" + s;
        T = [T; r(["chance_mm" "direct_err_mm"], [median(X.chanceMM(k)) median(X.directErrMM(k))], "mm", b)]; %#ok<AGROW>
        for snr = opts.SNRdB
            ks = k & X.snrDB == snr;  bs = b + "_" + i_snr(snr);
            T = [T; r(["core_err_mm" "net_ratio"], [median(X.coreErrMM(ks)) median(X.netRatio(ks))], ["mm" "ratio"], bs)]; %#ok<AGROW>
            for hf = 1:2
                T = [T; r("core_err_mm", median(X.coreErrMM(ks & X.half == hf)), "mm", bs + "_h" + hf)]; %#ok<AGROW>
            end
        end
    end
end

function s = i_snr(s)
    if isinf(s), s = "inf"; else, s = string(s) + "dB"; end
end

% Author: Diellor Basha, 2026
