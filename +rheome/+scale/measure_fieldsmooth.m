function [T, X, maps] = measure_fieldsmooth(name, S, opts)
% SCALE.MEASURE_FIELDSMOOTH  How smooth the cortical current, div and curl are, per vertex and band-limited.
%
%   [T, X] = rheome.scale.measure_fieldsmooth(name)
%   [T, X, maps] = rheome.scale.measure_fieldsmooth(name, S, CutoffsMM=[65 92 130], NumFrames=20)
%
% The diagnosis behind MS1 Fig. 12 (C-E): the arrows and the div/curl maps look noisy although the
% minimum-norm current is smooth at the scale MEG resolves. At the NumFrames strongest alpha moments of
% the clean span (8-12 Hz, 3rd-order Butterworth zero-phase, Hilbert at the sensors, real part; peaks of
% the channel-mean envelope, >= 1 s apart), left hemisphere, plain minimum norm (S.Res, free
% orientation, amplitude), it measures with rheome.flow.smoothness:
%   band "raw"         J (ambient 3-vector), Jt (its tangential reading = what the arrows draw),
%                      |J|, pointwise div and curl (rheome.differential.divergence / curl), and the
%                      Helmholtz potentials phi, psi
%   band "bl<cut>"     the same after rheome.flow.bandlimit at CutoffMM = cut: J, div, curl, phi, psi
% and, per frame, the shares of the current energy that are NORMAL to the cortex and harmonic.
%
% Returns rheome.scale.rows (analysis "fieldsmooth"): for each map <m> in J Jt absJ div curl phi psi,
%   <m>_wavelength_mm   median over frames of the Dirichlet wavelength                mm
%   <m>_coherence       (vectors) median magnitude-weighted neighbour cosine           cos
% and normal_frac, harmonic_frac (band ""), finest_sigma_mm per cutoff. X is the per-frame table.
% maps (third output, for the figure): the strongest frame's raw and band-limited fields at every cutoff,
% the hemisphere surface, the frame time.
%
% See also: rheome.flow.bandlimit, rheome.flow.smoothness, rheome.differential.helmholtzbands, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.CutoffsMM double = [65 92 130]
        opts.NumFrames (1,1) double {mustBeInteger, mustBePositive} = 20
        opts.Band (1,2) double = [8 12]
        opts.Hemi (1,1) string = "L"
        opts.EdgeS (1,1) double = 10
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    H = S.B.(char(opts.Hemi));  Sh = H.S;  lbo = H.lbo;
    rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
    K = S.Res.ImagingKernel(rows, :);

    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);
    [b, a] = butter(3, opts.Band / (fs/2), 'bandpass');
    Fa = hilbert(filtfilt(b, a, F.')).';
    env = mean(abs(Fa), 1);  env(CS.mask) = 0;  env([1:CS.first-1, CS.last+1:end]) = 0;
    [~, fr] = findpeaks(env, 'MinPeakDistance', round(fs), 'SortStr', 'descend', 'NPeaks', opts.NumFrames);
    J = K * real(Fa(:, fr));  clear Fa
    nF = numel(fr);

    n = Sh.VertNormals ./ max(vecnorm(Sh.VertNormals, 2, 2), eps);  Jt = J;
    jn = J(1:3:end,:).*n(:,1) + J(2:3:end,:).*n(:,2) + J(3:3:end,:).*n(:,3);
    for c = 1:3, Jt(c:3:end, :) = J(c:3:end, :) - jn .* n(:, c); end
    Hb = rheome.differential.helmholtzbands(J, Sh, lbo, Voices=2);
    raw = struct('J', J, 'Jt', Jt, 'absJ', sqrt(J(1:3:end,:).^2 + J(2:3:end,:).^2 + J(3:3:end,:).^2), ...
                 'div', rheome.differential.divergence(J, Sh), 'curl', rheome.differential.curl(J, Sh), ...
                 'phi', Hb.H.Phi, 'psi', Hb.H.Psi);
    X = table();  X = i_measure(X, raw, "raw", fr/fs, Sh, lbo);
    maps = struct('t0', fr(1)/fs, 'surface', Sh, 'gv', H.gv, 'hemi', opts.Hemi, 'raw', i_frame(raw, 1), 'bl', struct([]));
    finest = zeros(size(opts.CutoffsMM));
    for i = 1:numel(opts.CutoffsMM)
        Bl = rheome.flow.bandlimit(J, Sh, lbo, CutoffMM=opts.CutoffsMM(i));
        bl = struct('J', Bl.J, 'div', Bl.Div, 'curl', Bl.Curl, 'phi', Bl.Phi, 'psi', Bl.Psi);
        tag = "bl" + opts.CutoffsMM(i);
        X = i_measure(X, bl, tag, fr/fs, Sh, lbo);
        finest(i) = Bl.finestSigmaMM;
        m = i_frame(bl, 1);  m.cutoffMM = opts.CutoffsMM(i);  m.finestSigmaMM = Bl.finestSigmaMM;
        maps.bl = [maps.bl, m];
    end

    T = rheome.scale.rows("fieldsmooth", ["n_frames" "normal_frac" "harmonic_frac"], ...
        [nF median(Hb.normalFrac) median(Hb.harmFrac)], ["frames" "fraction" "fraction"]);
    T = [T; rheome.scale.rows("fieldsmooth", repmat("finest_sigma_mm", numel(finest), 1), finest(:), "mm", ...
         "bl" + string(opts.CutoffsMM(:)))];
    for g = unique(X.band, 'stable')'
        for m = unique(X.map(X.band == g), 'stable')'
            k = X.band == g & X.map == m;
            T = [T; rheome.scale.rows("fieldsmooth", m + "_wavelength_mm", median(X.wavelengthMM(k)), "mm", g)]; %#ok<AGROW>
            if any(~isnan(X.coherence(k)))
                T = [T; rheome.scale.rows("fieldsmooth", m + "_coherence", median(X.coherence(k)), "cos", g)]; %#ok<AGROW>
            end
        end
    end
    fprintf('[fieldsmooth %s] %d frames: Jt %.0f mm raw -> %.0f mm at %g mm cutoff; div %.0f -> %.0f mm; normal share %.2f\n', ...
        name, nF, median(X.wavelengthMM(X.band=="raw" & X.map=="Jt")), ...
        median(X.wavelengthMM(X.band=="bl"+opts.CutoffsMM(end) & X.map=="J")), opts.CutoffsMM(end), ...
        median(X.wavelengthMM(X.band=="raw" & X.map=="div")), ...
        median(X.wavelengthMM(X.band=="bl"+opts.CutoffsMM(end) & X.map=="div")), median(Hb.normalFrac));
end

function X = i_measure(X, maps, band, tS, Sh, lbo)
    for m = string(fieldnames(maps))'
        s = rheome.flow.smoothness(maps.(m), Sh, lbo);  n = numel(tS);
        X = [X; table(repmat(band, n, 1), repmat(m, n, 1), tS(:), s.wavelengthMM(:), s.coherence(:), ...
             'VariableNames', {'band', 'map', 'tS', 'wavelengthMM', 'coherence'})]; %#ok<AGROW>
    end
end

function m = i_frame(maps, k)
    m = structfun(@(x) x(:, k), maps, 'UniformOutput', false);
end

% Author: Diellor Basha, 2026
