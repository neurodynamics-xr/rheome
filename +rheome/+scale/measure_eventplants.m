function [T, X, Z] = measure_eventplants(name, S, opts)
% SCALE.MEASURE_EVENTPLANTS  Planted alpha episodes of known timing in the participant's own rest, boundaries read back (MS1 SPEC1 §8 A).
%
%   [T, X] = rheome.scale.measure_eventplants(name)
%   [T, X] = rheome.scale.measure_eventplants(name, S, SNRs=[5 10 20], NPer=20, Depth=8)
%
% The through-MEG test of rheome.detect.eventtable's boundaries, in the injection design of
% rheome.scale.measure_inject (same source, carrier and readout):
%   source    a Laplace-Beltrami heat-kernel blob (SigmaMM 15) along the normal at a random vertex, a 10 Hz
%             carrier of random phase, amplitude a(t) = S((t - tOn)/tauR) S(-(t - tOff)/tauF); duration
%             U(DurRange), tauR, tauF U(TauRange); forwarded through the participant's leadfield
%   background  a random 8 s segment of the participant's own clean rest (the episode centred in it)
%   amplitude  SNR (dB) = 20 log10 of the planted episode's full-height envelope on its tile (read alone,
%             through the same readout) over the median of the REAL envelope on that tile
%   readout   8-16 Hz band-pass and Hilbert on the segment, minimum norm, |J| per vertex, LB modes, the
%             depth-Depth tile mean (rheome.flow.envelopemodes' chain), at 300 frames/s
%   boundary  rheome.detect.eventtable on the planted tile, candidate [tOn, tOff] each moved by
%             U(-CandJitterS, CandJitterS) (a detector's boundary error), at TauMax and TauMaxAlt, with each of
%             risefall's Models (linear, quadrature, posterior) -- the evidence measure_events' Model is chosen by
% A plant PASSES if its fit counts (r2 >= R2Min), |tOn error| and |tOff error| <= TimeTileS (one time tile,
% 165 ms for alpha) and tauR, tauF within TauRatio of planted. The PARTICIPANT passes at an SNR if the
% medians over its plants meet the same rule (unfitted plants counted as failures of the timing medians:
% their errors are Inf). SPEC1 §8 A: the boundaries are accepted if >= 80% of participants pass at 10 dB
% (rheome.scale.groupevents).
%
% X (eventplants.csv): one row per hemisphere x plant x SNR x bound: the truth, the fit, errors, pass.
% T (rheome.scale.rows, analysis "eventplants"), band "<snr>dB/tau<bound>/<model>": frac_fit, on_err_med_s,
%   off_err_med_s, tauR_ratio_med, tauF_ratio_med, frac_plants_pass, participant_pass.
%
% ⚠ The candidate comes from the truth plus jitter, not from the tracker: this tests the boundary fit
%   through the instrument, not detection (which measure_inject tests).
% ⚠ The planted tile only: the envelope is the tile under the blob's peak vertex.
%
% See also: rheome.detect.eventtable, rheome.detect.risefall, rheome.scale.measure_inject, rheome.scale.measure_events
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depth (1,1) double = 8
        opts.SNRs double = [5 10 20]
        opts.NPer (1,1) double {mustBeInteger, mustBePositive} = 20
        opts.DurRange (1,2) double = [0.5 2]
        opts.TauRange (1,2) double = [0.1 0.5]
        opts.CandJitterS (1,1) double = 0.25
        opts.TauMax (1,1) double = 1
        opts.TauMaxAlt (1,1) double = 2
        opts.R2Min (1,1) double = 0.7
        opts.TimeTileS (1,1) double = 0.165
        opts.TauRatio (1,2) double = [0.8 1.25]
        opts.SigmaMM (1,1) double = 15
        opts.Hemis string = ["L" "R"]
        opts.Seed (1,1) double = 41
        opts.EdgeS (1,1) double = 10
        opts.Models string = ["linear" "quadrature" "posterior"]
        opts.Footprint (1,1) double = 0
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    rng(opts.Seed);
    FR = 300;  F0 = 10;  SEG = 8;  PAD = 1;                      % 8 s segment, readout over its central 6 s
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', 8, 'HalfPowerFrequency2', 16, 'SampleRate', fs);
    E = rheome.flow.envelopemodes(name, Band=[8 16], Rate=FR);
    nSeg = SEG*fs;  tSeg = (0:nSeg-1)/fs;  cIdx = PAD*fs + 1:round(fs/FR):(SEG - PAD)*fs;  tc = (cIdx - 1)/fs;
    Sg = @(x) 1 ./ (1 + exp(-x));  X = table();  Z = table();  keepZ = nargout > 2;
    for hh = opts.Hemis
        h = char(hh);  H = S.B.(h);  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  lam = H.lbo.Lambda(:);
        av = full(sum(Mm, 2));  gv = double(H.gv(:));  nV = numel(gv);  nrm = double(H.S.VertNormals);
        rows = reshape((gv' - 1)*3 + (1:3)', [], 1);  K = S.Res.ImagingKernel(rows, :);  GL = S.G(:, rows);
        hk = exp(-lam*(opts.SigmaMM*1e-3)^2/2);
        Tr = rheome.geom.tree(H.S, L=H.lbo.L, M=Mm, MaxDepth=opts.Depth);
        G = rheome.geom.tiles(Tr, H.S, opts.Depth);  P = double(G.P);  Wd = (P'*(av.*Phi)) ./ (P'*av);
        Yt = Wd*double(E.(h).C);  med = median(Yt(:, ceil(CS.first*FR/fs):floor(CS.last*FR/fs)), 2);  clear Yt
        A = struct(h, P'*av);
        for q = 1:opts.NPer
            v0 = randi(nV);  k0 = G.tileOf(v0);
            d = opts.DurRange(1) + rand*diff(opts.DurRange);  tau = opts.TauRange(1) + rand(1, 2)*diff(opts.TauRange);
            on = SEG/2 - d/2;  off = SEG/2 + d/2;
            env = Sg((tSeg - on)/tau(1)) .* Sg(-(tSeg - off)/tau(2));
            g0 = GL*i_blob(v0, Phi, Mm, hk, nrm);
            bsig = g0 * (env .* cos(2*pi*F0*tSeg + 2*pi*rand));
            s0 = randi([CS.first, CS.last - nSeg + 1]) + (0:nSeg-1);
            U = Wd * i_envmodes(bsig, bp, cIdx, K, Phi, Mm);  u = U(k0, :);  pk = max(U, [], 2) / max(u);
            fp = find(pk >= opts.Footprint);  if opts.Footprint <= 0 || ~ismember(k0, fp), fp = k0; end
            zk = [k0; setdiff(find(pk >= 0.25), k0)];
            j = opts.CandJitterS*(2*rand(1, 2) - 1);
            Ev = table(hh, {fp}, on + j(1), off + j(2), 'VariableNames', {'hemi','tiles','t0','t1'});
            if keepZ
                n0 = Wd(zk, :) * i_envmodes(F(:, s0), bp, cIdx, K, Phi, Mm);  yZ = cell(1, numel(opts.SNRs));
            end
            for snr = opts.SNRs
                amp = 10^(snr/20) * med(k0) / max(u);
                Y = struct(h, zeros(size(Wd, 1), numel(cIdx)));
                Yk = Wd(union(fp, zk), :) * i_envmodes(F(:, s0) + amp*bsig, bp, cIdx, K, Phi, Mm);
                Y.(h)(union(fp, zk), :) = Yk;
                if keepZ, yZ{snr == opts.SNRs} = Y.(h)(zk, :) / (amp * max(u)); end   % in units of the plant's height
                for tm = [opts.TauMax opts.TauMaxAlt], for md = opts.Models
                    R = rheome.detect.eventtable(Ev, Y, A, FR, T0=tc(1), TauMax=tm, R2Min=opts.R2Min, Model=md);
                    eOn = R.tOn - on;  eOff = R.tOff - off;  rR = R.tauR / tau(1);  rF = R.tauF / tau(2);
                    ps = R.fit && abs(eOn) <= opts.TimeTileS && abs(eOff) <= opts.TimeTileS && ...
                         all([rR rF] >= opts.TauRatio(1) & [rR rF] <= opts.TauRatio(2));
                    X = [X; table(hh, q, snr, tm, md, on, off, tau(1), tau(2), R.tOn, R.tOff, R.tauR, R.tauF, R.r2, R.fit, ...
                        R.atBound, eOn, eOff, rR, rF, ps, 'VariableNames', {'hemi','plant','snr_db','tau_max','model','tOn_true', ...
                        'tOff_true','tauR_true','tauF_true','tOn','tOff','tauR','tauF','r2','fit','atBound','on_err_s', ...
                        'off_err_s','tauR_ratio','tauF_ratio','pass'})]; %#ok<AGROW>
                end, end
            end
            if keepZ
                Z = [Z; table(hh, q, on, off, tau(1), tau(2), on + j(1), off + j(2), {U(zk, :) / max(u)}, {n0 / max(u)}, {yZ}, ...
                    {10.^(opts.SNRs/20) * med(k0) / max(u)}, {zk}, {pk(zk)}, {A.(h)(zk)}, 'VariableNames', {'hemi','plant', ...
                    'tOn_true','tOff_true','tauR_true','tauF_true','t0','t1','u','n','y','amp','tiles','peakfrac','area'})]; %#ok<AGROW>
            end
        end
        fprintf('[eventplants %s] %s: %d plants x %d SNRs\n', name, h, opts.NPer, numel(opts.SNRs));
    end
    if keepZ, Z.Properties.UserData = struct('t', tc, 'fs', FR, 'snr_db', opts.SNRs); end
    T = table();
    for snr = opts.SNRs
        for tm = [opts.TauMax opts.TauMaxAlt], for md = opts.Models
            q = X(X.snr_db == snr & X.tau_max == tm & X.model == md, :);
            eOn = abs(q.on_err_s);  eOff = abs(q.off_err_s);  eOn(~q.fit) = Inf;  eOff(~q.fit) = Inf;
            rR = q.tauR_ratio;  rF = q.tauF_ratio;  rR(~q.fit) = NaN;  rF(~q.fit) = NaN;
            m = [mean(q.fit) median(eOn) median(eOff) median(rR, 'omitnan') median(rF, 'omitnan') mean(q.pass)];
            pp = m(2) <= opts.TimeTileS && m(3) <= opts.TimeTileS && all(m(4:5) >= opts.TauRatio(1) & m(4:5) <= opts.TauRatio(2));
            T = [T; rheome.scale.rows("eventplants", ["frac_fit" "on_err_med_s" "off_err_med_s" "tauR_ratio_med" ...
                 "tauF_ratio_med" "frac_plants_pass" "participant_pass"], [m pp], ...
                 ["fraction" "s" "s" "ratio" "ratio" "fraction" "bool"], sprintf('%gdB/tau%g/%s', snr, tm, md))]; %#ok<AGROW>
        end, end
    end
end

% Heat-kernel blob at vertex v along the normal, interleaved [x1;y1;z1;x2;...] (as measure_inject).
function j = i_blob(v, Phi, Mm, hk, nrm)
    d = zeros(size(Phi, 1), 1);  d(v) = 1;
    a = Phi*(hk.*(Phi'*(Mm*d)));  a = a/max(a);
    j = reshape((a.*nrm)', [], 1);
end

% The envelope of a sensor segment as LB mode coefficients at the readout frames (as measure_inject).
function Cw = i_envmodes(Fs, bp, cIdx, K, Phi, Mm)
    Fa = hilbert(filtfilt(bp, Fs.')).';
    Cw = Phi'*(Mm*rheome.flow.activation(K*Fa(:, cIdx)));
end

% Author: Diellor Basha, 2026
