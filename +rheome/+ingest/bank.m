function [fb, bands, frame] = bank(nT, fs, cfg)
% INGEST.BANK  The full-range constant-Q bank of a record, and its scales grouped into octave bands.
%
%   [fb, bands, frame] = rheome.ingest.bank(nT, fs, cfg)
%
% TWO BANKS. cfg.Bank = "frame" (default): a timefilterbank -- a designed tight constant-Q
% frame on the log-frequency axis anchored at cfg.Anchor, evaluated per member at its own
% rate; bands are the integer octaves of the anchor plus a "below" (low-pass) and an
% "above" (high-pass) band, so the bands partition sum x^2 exactly (A = B = 1). The rest
% of this header describes cfg.Bank = "morse": MATLAB's cwtfilterbank at the full rate,
% the reference bank, 7x slower. Both return the same bands / frame structures.
%
% fb is MATLAB's cwtfilterbank at the record's own length. With cfg.FrequencyLimits empty
% (the default) the bottom is the wavelet's own limit (cwtfreqbounds: where the support
% still fits in the record), raised to support/cfg.MaxSupport if that is set, and the
% top is THE HIGHEST SCALE WHOSE COEFFICIENTS ARE REPRODUCIBLE, see below. Nothing else
% is truncated (design §3.2). Each scale is a full-rate pass, so the floor is the build
% time: 151 scales to 0.0055 Hz, 76 to 1 Hz, at 600 Hz.
%
% ⚠ THE TOP SCALE IS CAPPED WHERE THE WAVELET'S GAIN AT NYQUIST FALLS TO 1e-6, about 0.7
% of cwtfreqbounds' default. Measured (docs/2026-09-22-ingest-notes.md): a scale whose
% response is still appreciable at Nyquist is not band-limited on the discrete grid, and
% its coefficients then DEPEND ON THE RECORD LENGTH -- the same filter (gain curves equal
% to 1e-4) applied to the same samples in a 4000-sample record and in a 256-sample page
% differs by 127 % at the default top (gain 0.5 at Nyquist), 1.4 % at gain 2.6e-3, 5e-6 at
% gain 2e-6. Below the cap, paged and whole-record coefficients agree to single precision.
% The cap is not a choice about what is interesting: the energy above it stays in sumX2
% and in the residual; only its band decomposition is withheld, because no two builds
% would agree on it. At 600 Hz the top is ~182 Hz (was 260), at 2400 Hz ~730 (was 1042).
%
% THE MASTER GRID IS ANALYTIC: fc(k) = fmax * 2^(-k/V). It equals centerFrequencies of
% the bank to 6e-15, and a SUB-BANK whose upper limit is a master scale reproduces the
% master's scales at ANY SignalLength and ANY fs (4e-16) -- which is what lets a page or
% a decimated copy carry the same band axis (rheome.ingest.reducepaged).
%
% BANDS are groups of VoicesPerOctave consecutive scales counted DOWN from the top scale.
% Each scale owns half a voice on either side, so a full band spans exactly one octave and
% neighbouring bands are contiguous in log-frequency:
%
%     fLo = fc(lowest scale) * 2^(-1/(2V)),   fHi = fc(highest scale) * 2^(1/(2V))
%
% The leftover scales at the bottom form one PARTIAL band (.partial = true, extent < 1).
%
% MEASURED EXTENTS COME FROM A PROBE and are scale-covariant. powerbw, waveletsupport and
% freqz on a 1.44 M-sample bank would materialise gigabytes, so they are read from a probe
% bank of min(nT, 2^16) samples at the same fs, and the per-scale constants are fitted on
% its untruncated scales: support*fc = 6.46 (spread 0.9 %), lo/fc = 0.893, hi/fc = 1.107
% (5e-5). Every master scale gets its extents from those constants, with the time support
% CAPPED at the record duration (waveletsupport truncates at the boundary anyway).
%
% ⚠ frame.A and frame.B are the frame bounds of S(f)/2, S(f) = sum_m |psi_m(f)|^2, on the
% INTERIOR of the probe's range. THE HALF IS THE ANALYTIC CONVENTION: freqz is one-sided
% with gain 2, so for a REAL signal a tone at f0 gives sum_m |W_m|^2 = S(f0)/4 per sample
% against x^2 = 1/2, i.e. S(f0)/2 times sum x^2. Measured: 6.60 against S/2 = 6.61. S is
% flat on the interior (B/A = 1.0003 at 10 voices), so the probe's A, B are the record's.
%
% ⚠ frame.coiConst: the cone of influence is coi(t) = min(fmax, coiConst / dist(t)) with
% dist the seconds to the nearer record edge. Fitted on the probe (1.6475 for Morse at 10
% voices, spread 3e-4). Both reduce paths count nCoi from it.
%
% INPUTS:
%   nT, fs   record length (samples) and rate (Hz)
%   cfg      rheome.ingest.config
% OUTPUT:
%   fb     cwtfilterbank at SignalLength = nT (lazy; used by the whole-record path)
%   bands  table, top band first: j, scales {master indices, highest first}, fLo, fHi,
%          fCenter, fExtent (octaves), fLoMeasured, fHiMeasured, tSupport (s), partial,
%          naturalLevel (first level whose tile >= tSupport; the store's diagonal)
%   frame  .fc [1 x nF] master scales  .nF  .fmin .fmax  .f .S (probe, raw sum |psi|^2)
%          .A .B (bounds of S/2)  .interior  .Q  .coiConst  .support (per-scale constants)
%          .probeLength
%
% See also: rheome.ingest.config, rheome.ingest.reduce, rheome.ingest.reducepaged, cwtfilterbank
%
% Author: Diellor Basha, 2026

    arguments
        nT  (1,1) double {mustBeInteger, mustBePositive}
        fs  (1,1) double {mustBePositive}
        cfg (1,1) struct
    end
    if isfield(cfg, 'Bank') && strcmp(cfg.Bank, 'frame')
        [fb, bands, frame] = i_framebank(nT, fs, cfg);
        return
    end
    NYQUIST_GAIN = 1e-6;

    V = cfg.VoicesPerOctave;
    common = {'SamplingFrequency', fs, 'Wavelet', cfg.Wavelet, 'VoicesPerOctave', V};

    % the top scale: where the wavelet's gain at Nyquist is negligible (NyquistGain)
    [fminD, fmaxD] = cwtfreqbounds(nT, fs, 'Wavelet', cfg.Wavelet);
    Np = min(nT, 2^16);
    p0 = cwtfilterbank('SignalLength', Np, common{:});          % probe at the default limits
    [H0, f0] = freqz(p0);
    fcp0 = centerFrequencies(p0);
    % r* = f/fc where a filter's gain has fallen to NYQUIST_GAIN on its high side. Read off
    % a MID filter, whose whole response lies inside [0, fs/2]; scale covariance makes the
    % ratio the same for every scale (the top filter itself never gets there before Nyquist).
    m0 = max(1, round(numel(fcp0) / 2));
    gm = H0(m0, :) / max(H0(m0, :));
    above = f0(:)' > fcp0(m0);
    rStar = min(f0(above & gm < NYQUIST_GAIN)) / fcp0(m0);
    if isempty(rStar), rStar = 2; end
    fmaxCap = (fs / 2) / rStar;
    % the floor: the record (cwtfreqbounds) or cfg.MaxSupport, whichever is higher
    ws0 = waveletsupport(p0);
    ok0 = ~isnan(ws0.TimeSupport) & ws0.TimeSupport < 0.25 * Np / fs;
    if nnz(ok0) < 3, ok0 = ~isnan(ws0.TimeSupport); end
    cS0 = median(ws0.TimeSupport(ok0) .* fcp0(ok0));
    if isempty(cfg.FrequencyLimits)
        fmin = fminD;  fmax = min(fmaxD, fmaxCap);
        if isfinite(cfg.MaxSupport), fmin = max(fmin, cS0 / cfg.MaxSupport); end
    else
        fmin = cfg.FrequencyLimits(1);  fmax = cfg.FrequencyLimits(2);
        if fmax > fmaxCap * (1 + 1e-9)
            warning('ingest:bank:nyquist', ...
                    'FrequencyLimits(2) = %.4g Hz exceeds the reproducible top %.4g Hz (gain %.0e at Nyquist); coefficients there depend on record length.', ...
                    fmax, fmaxCap, NYQUIST_GAIN);
        end
    end
    fb = cwtfilterbank('SignalLength', nT, common{:}, 'FrequencyLimits', [fmin fmax]);
    fc = centerFrequencies(fb);  fc = fc(:)';          % descending; equals the analytic grid
    nF = numel(fc);
    fmax = fc(1);  fmin = min(fmin, fc(end));

    % probe: same fs, wavelet and top, short enough to measure on
    if Np == nT
        p = fb;
    else
        pmin = max(fmin, cwtfreqbounds(Np, fs, 'Wavelet', cfg.Wavelet));
        p = cwtfilterbank('SignalLength', Np, common{:}, 'FrequencyLimits', [pmin fmax]);
    end
    fcp = centerFrequencies(p);  fcp = fcp(:)';
    pb  = powerbw(p);
    ws  = waveletsupport(p);
    ok  = ~isnan(ws.TimeSupport) & ws.TimeSupport < 0.25 * Np / fs;   % untruncated scales
    if nnz(ok) < 3, ok = ~isnan(ws.TimeSupport); end
    cS  = median(ws.TimeSupport(ok) .* fcp(ok)');
    cLo = median(pb.LowFrequencyBorder(ok)  ./ fcp(ok)');
    cHi = median(pb.HighFrequencyBorder(ok) ./ fcp(ok)');
    support = min(cS ./ fc, nT / fs);
    loM = cLo * fc;  hiM = cHi * fc;

    % bands
    nB = ceil(nF / V);
    half = 2^(1/(2*V));
    j = (1:nB)';
    scales = cell(nB, 1);
    [fLo, fHi, fLoM, fHiM, tSup] = deal(zeros(nB, 1));
    partial = false(nB, 1);
    for b = 1:nB
        s = (b-1)*V + (1:V);  s = s(s <= nF);
        scales{b} = s;
        fLo(b)  = fc(s(end)) / half;
        fHi(b)  = fc(s(1)) * half;
        fLoM(b) = min(loM(s));
        fHiM(b) = max(hiM(s));
        tSup(b) = max(support(s));
        partial(b) = numel(s) < V;
    end
    fCenter = sqrt(fLo .* fHi);
    fExtent = log2(fHi ./ fLo);
    % natural level: the first grid level whose tile is at least the band's support
    % (design §3.13). Below it a tile's coefficient energy is not that tile's sample
    % energy (the residual is negative there), so the store carries the band from here up.
    F = cfg.FrameFloor;
    naturalLevel = max(0, ceil(log2(tSup / F)));
    bands = table(j, scales, fLo, fHi, fCenter, fExtent, fLoM, fHiM, tSup, partial, naturalLevel, ...
                  'VariableNames', {'j','scales','fLo','fHi','fCenter','fExtent', ...
                                    'fLoMeasured','fHiMeasured','tSupport','partial','naturalLevel'});

    % frame constants from the probe
    [H, f] = freqz(p);
    f = f(:)';
    S = sum(abs(H).^2, 1);
    interior = [min(fcp) * sqrt(2), max(fcp) / sqrt(2)];
    in = f >= interior(1) & f <= interior(2);
    [~, ~, coi] = wt(p, zeros(Np, 1));
    d = (0:Np-1) / fs;  d = min(d, (Np-1)/fs - d);
    mid = d > 0.05 * Np / fs & d < 0.4 * Np / fs;
    coiConst = median(coi(mid)' .* d(mid));
    frame = struct('fc', fc, 'nF', nF, 'fmin', fmin, 'fmax', fmax, 'f', f, 'S', S, ...
                   'A', min(S(in))/2, 'B', max(S(in))/2, 'interior', interior, ...
                   'Q', qfactor(p), 'coiConst', coiConst, 'wavelet', cfg.Wavelet, ...
                   'support', struct('timeTimesFc', cS, 'loOverFc', cLo, 'hiOverFc', cHi), ...
                   'probeLength', Np);
end

function [fb, bands, frame] = i_framebank(nT, fs, cfg)
% the timefilterbank path: bands are integer octaves of the anchor, plus the two edge bands
    V = cfg.VoicesPerOctave;  a = cfg.Anchor;
    lim = cfg.FrequencyLimits;
    fb = rheome.timefilterbank(nT, 'SamplingFrequency', fs, 'VoicesPerOctave', V, 'Anchor', a, ...
                        'FrequencyLimits', lim, 'Oversample', cfg.Oversample);
    sup = support(fb);
    if isfinite(cfg.MaxSupport)
        fLo = sup.timeTimesFc / cfg.MaxSupport;
        if isempty(lim), lim = [fLo, fs/2]; else, lim(1) = max(lim(1), fLo); end
        fb = rheome.timefilterbank(nT, 'SamplingFrequency', fs, 'VoicesPerOctave', V, 'Anchor', a, ...
                            'FrequencyLimits', lim, 'Oversample', cfg.Oversample);
    end
    fc = centerFrequencies(fb);                              % ascending, edges included
    kd = kinds(fb);
    M = fb.NumMembers;
    isBand = strcmp(kd, 'band');
    k = round(V * log2(fc(isBand) / a));                    % voice index of each band member
    oct = floor(k / V);                                      % integer octave of the anchor
    octs = unique(oct);
    nB = numel(octs) + 2;
    scales = cell(nB, 1);  kind = cell(nB, 1);
    [fLo, fHi, fLoM, fHiM, tSup] = deal(zeros(nB, 1));  partial = false(nB, 1);
    % row 1: the "above" band (high-pass member), then octaves descending, then "below"
    scales{1} = M;  kind{1} = 'above';
    fLo(1) = fc(M-1) ;  fHi(1) = fs/2;  fLoM(1) = fc(M-1);  fHiM(1) = fs/2;  tSup(1) = sup.timeTimesFc / fc(M-1);
    bandIdx = find(isBand);
    for i = 1:numel(octs)
        o = octs(end - i + 1);  r = i + 1;
        mem = bandIdx(oct == o);                             % ascending
        scales{r} = fliplr(mem(:)');                         % top first, so scales{b}(end) is the lowest
        kind{r} = 'octave';
        fLo(r) = a * 2^o;  fHi(r) = a * 2^(o+1);
        kk = k(oct == o);
        % the outermost octaves are partial: they hand over to the edge bands at the
        % outermost member centres, so the nominal extents stay contiguous
        if o == octs(end), fHi(r) = min(fHi(r), fc(M-1)); end
        if o == octs(1),   fLo(r) = max(fLo(r), fc(2)); end
        fLoM(r) = a * 2^((min(kk) - 1) / V);  fHiM(r) = a * 2^((max(kk) + 1) / V);   % one voice of overhang each side
        tSup(r) = sup.timeTimesFc / (a * 2^(min(kk) / V));
        partial(r) = numel(mem) < V || fHi(r) < a * 2^(o+1) || fLo(r) > a * 2^o;
    end
    scales{nB} = 1;  kind{nB} = 'below';
    fLo(nB) = fs / nT;  fHi(nB) = fc(2);  fLoM(nB) = fs / nT;  fHiM(nB) = fc(2);  tSup(nB) = nT / fs;
    fCenter = sqrt(fLo .* fHi);
    fExtent = log2(fHi ./ fLo);
    naturalLevel = max(0, ceil(log2(tSup / cfg.FrameFloor)));
    j = (1:nB)';
    bands = table(j, scales, fLo, fHi, fCenter, fExtent, fLoM, fHiM, tSup, partial, naturalLevel, kind, ...
                  'VariableNames', {'j','scales','fLo','fHi','fCenter','fExtent', ...
                                    'fLoMeasured','fHiMeasured','tSupport','partial','naturalLevel','kind'});
    [H, f] = freqz(fb);
    S = sum(H.^2, 1);
    b = framebounds(fb);
    frame = struct('fc', fc, 'nF', M, 'fmin', fc(2), 'fmax', fs/2, 'f', f, 'S', S, ...
                   'A', b.A, 'B', b.B, 'interior', [fc(2), fc(M-1)], 'Q', sup.Q, ...
                   'coiConst', NaN, 'wavelet', 'frame', 'bank', 'frame', ...
                   'support', struct('timeTimesFc', sup.timeTimesFc, 'loOverFc', sup.loOverFc, 'hiOverFc', sup.hiOverFc), ...
                   'probeLength', 2^16);
end
% Author: Diellor Basha, 2026
