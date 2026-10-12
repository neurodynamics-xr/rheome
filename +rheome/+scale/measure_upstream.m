function [T, X] = measure_upstream(name, S, opts)
% SCALE.MEASURE_UPSTREAM  The AGE-1 participant readout: U, U_surrogate, U_total, U_of and their covariates per recording.
%
%   [T, X] = rheome.scale.measure_upstream(name)
%   [T, X] = rheome.scale.measure_upstream(name, rheome.scale.sensors(name), Gradient='/path/g1_fsaverage.mat')
%   [T, X] = rheome.scale.measure_upstream(name, S, Passes=["periodic" "surrogate" "total"])
%   [T, X] = rheome.scale.measure_upstream(name, S, Data=Fplanted)    % a planted copy of the S.iSel rows
%   [T, X] = rheome.scale.measure_upstream(name, S, Readout="v2")     % tile values by fwdLS (see READOUT)
%
% Why: whether large-scale slow-band activity at rest runs up the sensory -> association gradient, and whether
% that survives the instrument, is read per participant by ONE frozen definition: rheome.flow.upstreamshare,
% "participant readout (version 1)" (rheome 2a2fe62; TP2 gate S). This function is that readout from sensors to
% one row per recording, steps a-f of upstreamshare's help, and nothing in it is chosen per recording.
%
%   a. SPAN       rheome.scale.cleanspan; at least MinSpanS (180 s), else the recording is not read. Frames on
%                 an interior artefact block (cleanspan's .mask) cast no vote. Halves (gate G): the first and
%                 second half of the span, each through b-f on its own (its own aperiodic fit).
%   b. BAND       per channel, Welch PSD of the span (8 s Hann, 50 %; rheome.scale.cleanwelch), aperiodic fit with
%                 a knee on 1-45 Hz without 60/120 +- 2 Hz (rheome.spectral.decompose, fitPower = that PSD in the
%                 FFT's units, x fs N / 2), gain h on the span's FFT; the ANALYTIC slow band is the inverse FFT
%                 of 2 h C on the positive bins in [1 8] Hz, taken at Rate (50) frames per second.
%   c-d. FIELD    the unconstrained minimum-norm current (S.Res), all three components, as the area-weighted
%                 mean over each tile of rheome.geom.tree at DEPTH 3 per hemisphere (one [tiles x channels]
%                 operator per component: K then the tile mean, exactly linear); Normals = tile mean normal.
%   e. GRADIENT   the Margulies et al. 2016 principal gradient (G1, normative) on the FreeSurfer sphere,
%                 carried to the recording's cortex through its registered sphere (nearest vertex); signed so
%                 pericalcarine, precentral and postcentral (Desikan-Killiany) are low and inferiorparietal and
%                 medialorbitofrontal high; g = tile mean.
%                 Readout="v2" ("fwdLS") replaces this step only, see READOUT below.
%   f. RESOLUTION own-tile fraction (rheome.flow.crosstalk on the kernel's three components) >= 0.25 votes; a
%                 recording whose keptWeight < 0.5 in either hemisphere is not read. U pools the hemispheres
%                 by their vote weight. The mask and the exclusion depend on the kernel only, so every pass
%                 reads the same tiles of the same participants.
%
% READOUT (Readout=; "v1" by default, frozen):
%   "v1"  steps c-d as above: the MNE tile means of the three components.
%   "v2"  "fwdLS": the tile values are the regularised least-squares amplitudes of the tiles' OWN forward fields,
%         q = (Gp'Gp + mu I)^-1 Gp' Za, Gp = [L | R] the sensor field of uniform normal current over each depth-3
%         tile (the v1 null's generators), mu = 1e-3 trace(Gp'Gp)/nTiles; upstreamshare on the scalar q per
%         hemisphere. Everything else (a, b, e, f -- the mask and exclusion stay the MNE kernel's) is v1's.
%         Why: a folded depth-3 tile carrying uniform current is near-silent at the sensors, so the MNE tile-mean
%         kernel has almost no zero order and each tile's phase follows a kernel centroid hundreds of mm away.
%         Gate S on the Margulies G1 axis (planted plane waves, 0 and 10 dB, MEG): v1 reads the planted direction
%         reversed at lambda 200-500 mm and right only at 1000 mm; v2 reads it right at every lambda 200-1000 mm.
%         Every pass that reads tiles (periodic, surrogate, total) uses the same readout; "of" is unchanged.
%         ⚠ v2 IS A GATE S DEVIATION awaiting approval: it is not the frozen readout. Plants and readout share one
%         gain (the inverse crime), so its gate S margin is an upper bound; it aliases at lambda 100 mm.
%
% PASSES (Passes=; all four by default; each is one extra reading of the same span, tiles, mask and gradient):
%   "periodic"   U: h = hper, the periodic gain (the primary readout)
%   "surrogate"  U_surrogate: gate N. The periodic slow band of step b (real, at Rate) replaced by its
%                lagged-source null (rheome.forward.laggedrest: in 4 s windows, the frontal and posterior tile
%                generators of each hemisphere refitted with no travel, through the participant's own gain, at
%                the window's lag and ratio, + the phase-randomised residual), made analytic on the same bins,
%                and read by c-f. U* = U - U_surrogate. Read on the whole span and on each half.
%                Null="v1" (default, frozen, NOT VALIDATED: gate N failed 2026-10-11) is the above. Null="v2" (the
%                redesign): candidates are DEPTH-6 tiles (64 per hemisphere), the two
%                fitted generators keep their own time courses (Keep), and the residual keeps its zero-lag
%                spatial structure and loses only its lags (Residual "zerolag"); see rheome.forward.laggedpair.
%                ⚠ v2 IS NOT VALIDATED EITHER (2026-10-12, measure_upstreamnull on a held-out OMEGA MEG head):
%                lagged 20/40 mm discs leave the 95th-pct |U*| at 0.07-0.13 (readout v1) and 0.10-0.20
%                (v2 fwdLS), 13-44 % of the lag cancelled. With the TRUE disc generators as candidates the same
%                null cancels 87-97 % (95th pct 0.007-0.018): the limit is that one tile does not fit a compact
%                source's shape, not the residual. Do not read U* as travel.
%   "total"      U_total: h = 1, the periodic gain switched off (WAVE-OR-LEAKAGE H3b, Delta_ap = U_total - U).
%                The aperiodic fit still runs and gives the covariates. Whole span only.
%   "of"         U_of: the comparator (P2). Horn-Schunck optical flow (rheome.flow.apparent, Alpha 1) of the
%                periodic slow-band envelope |J| per vertex and hemisphere; the share of sum |v| |grad g| a over
%                vertices and frame pairs with v . grad g > 0. Whole span only. Never pooled with U.
%                ⚠ solved in chunks of OfChunk frames (each chunk normalised on its own by opticalflow_scalar).
%
% COVARIATES (from the aperiodic fit of the span, whatever the pass): aperiodic_exponent and aperiodic_offset
% (channel means, offset in log10 PSD units), snr_slow_dB (periodic / aperiodic Welch power in 1-8 Hz, per
% channel in dB, channel mean), aperiodic_ok (gate P: every channel's fit finite). r50 is the resolution
% analysis's (rheome.scale.measure_resolution), not repeated here.
%
% OUTPUT
%   T  metrics rows (analysis "upstream"): read (1/0), span_s, U, U_surrogate, U_total, U_of, U_star, and the
%      halves U_h1, U_h2, U_surrogate_h1, U_surrogate_h2, U_star_h1, U_star_h2; aperiodic_exponent,
%      aperiodic_offset, snr_slow_dB, aperiodic_ok (whole span; band "h1"/"h2" for the halves);
%      kept_weight_L/_R, g1_tilecheck_L/_R, n_frames
%   X  one row per (span, pass, hemisphere|"both"): span, pass, hemi, U, weight, kept_weight, n_edges,
%      n_frames, reason ("" when read; why the recording was not read otherwise)
%
% ⚠ WHAT IT CANNOT RESOLVE (state with any result): travel slower than 2 d f between tile centres reads
%   backwards (0.24 m/s at 1 Hz, 1.9 m/s at 8 Hz for d = 120 mm); 8 tiles per hemisphere make U a lobe-scale
%   share; zero-lag point spread pulls U to 0.5, lagged generators that do not travel can bias it (hence U*).
% ⚠ GATE G READS ONLY THE ICC AND SD OF U* (rheome.scale.groupupstream), never a mean, before the plans' approval.
%
% See also: rheome.flow.upstreamshare, rheome.forward.laggedrest, rheome.flow.crosstalk, rheome.scale.run,
%           rheome.scale.groupupstream
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Passes string {mustBeMember(opts.Passes, ["periodic" "surrogate" "total" "of"])} = ["periodic" "surrogate" "total" "of"]
        opts.Halves (1,1) logical = true
        opts.Gradient char = ''              % default: RHEOME_G1
        opts.SignCheck (1,1) logical = true
        opts.MinSpanS (1,1) double = 180
        opts.Seed (1,1) double = 1
        opts.OfChunk (1,1) double {mustBeInteger, mustBePositive} = 500
        opts.OwnMin (1,1) double = 0.25     % v1: 0.25; lower it only on a toy head (tests)
        opts.Data double = []               % [nSel x nT] read instead of the recording's S.iSel rows (planted copies)
        opts.Null (1,1) string {mustBeMember(opts.Null, ["v1" "v2"])} = "v1"   % the surrogate pass's null
        opts.Readout (1,1) string {mustBeMember(opts.Readout, ["v1" "v2"])} = "v1"   % the tile values (READOUT)
    end
    V1 = struct('Depth', 3, 'Band', [1 8], 'Rate', 50, 'FitRange', [1 45], 'LineHz', 60, 'SegS', 8, ...
                'OwnMin', opts.OwnMin, 'KeptMin', 0.5, 'EdgeS', 10, 'WindowS', 4, 'MaxLagS', 0.5, 'OfAlpha', 1);   % readout v1, frozen
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    if ~isempty(opts.Data)
        assert(isequal(size(opts.Data), size(F)), 'scale:upstream:data', 'Data must be [%d x %d], the recording''s selected rows', size(F, 1), size(F, 2));
        F = opts.Data;
    end
    CS = rheome.scale.cleanspan(F, fs, EdgeS=V1.EdgeS);
    spanS = (CS.last - CS.first + 1) / fs;
    H = i_geometry(name, S, V1, opts);
    T = rheome.scale.rows("upstream", ["span_s" "kept_weight_L" "kept_weight_R" "g1_tilecheck_L" "g1_tilecheck_R" "readout_version"], ...
        [spanS H.L.keptWeight H.R.keptWeight H.L.tilecheck H.R.tilecheck str2double(extractAfter(opts.Readout, 1))], ...
        ["s" "fraction" "fraction" "flag" "flag" "version"]);
    X = table();
    why = "";
    if spanS < opts.MinSpanS, why = sprintf("clean span %.1f s < %g s", spanS, opts.MinSpanS); end
    for h = ["L" "R"]
        if why == "" && H.(h).keptWeight < V1.KeptMin, why = sprintf("kept weight %.3f < %g in %s", H.(h).keptWeight, V1.KeptMin, h); end
    end
    if why ~= ""                                                  % not read: recorded before any U exists
        T = [T; rheome.scale.rows("upstream", "read", 0, "flag")];
        X = i_row("whole", "none", "both", NaN, NaN, NaN, 0, 0, why);
        fprintf('[upstream %s] not read: %s\n', name, why);  return
    end
    T = [T; rheome.scale.rows("upstream", "read", 1, "flag")];

    mid = floor((CS.first + CS.last) / 2);
    spans = {"whole", CS.first:CS.last};
    if opts.Halves, spans = [spans; {"h1", CS.first:mid}; {"h2", mid+1:CS.last}]; end
    U = struct();
    for i = 1:size(spans, 1)
        sn = spans{i, 1};  sp = spans{i, 2};  bandName = "";  sfx = "";
        if sn ~= "whole", bandName = sn;  sfx = "_" + sn; end
        try
            B = i_band(F, fs, sp, CS.mask, V1);
        catch ME                                                      % no clean Welch segment: this span is not read
            if ~startsWith(ME.identifier, "scale:cleanwelch:"), rethrow(ME); end
            X = [X; i_row(sn, "none", "both", NaN, NaN, NaN, 0, 0, string(ME.message))]; %#ok<AGROW>
            if sn == "whole", T.value(T.metric == "read") = 0; end
            fprintf('[upstream %s] %s not read: %s\n', name, sn, ME.message);
            if sn == "whole", return, else, continue, end
        end
        T = [T; rheome.scale.rows("upstream", ["aperiodic_exponent" "aperiodic_offset" "snr_slow_dB" "aperiodic_ok"], ...
             [B.exponent B.offset B.snr B.ok], ["chi" "log10 PSD" "dB" "flag"], bandName)]; %#ok<AGROW>
        for p = opts.Passes(:)'
            if sn ~= "whole" && ismember(p, ["total" "of"]), continue; end
            switch p
                case "periodic",  Za = i_analytic(B, B.hper, V1);  ok = B.valid;
                case "total",     Za = i_analytic(B, ones(size(B.hper)), V1);  ok = B.valid;
                case "surrogate"
                    Zp = i_analytic(B, B.hper, V1);
                    rs = rng;  rng(opts.Seed + i, 'twister');            % reproducible null, global stream restored
                    Nl = H.null;
                    [Ys, li] = rheome.forward.laggedrest(real(Zp), Nl.Gp, Nl.front, B.rate, WindowS=V1.WindowS, ...
                        MaxLagS=V1.MaxLagS, Group=Nl.group, Mask=~B.valid, Keep=Nl.keep, Residual=Nl.residual);
                    rng(rs);
                    Za = i_reanalytic(Ys, B.rate, V1.Band);  ok = B.valid & li.read;
                case "of"
                    Za = i_analytic(B, B.hper, V1);
                    uo = i_of(S, H, Za, B.valid, B.rate, V1, opts.OfChunk);
                    U.("U_of") = uo;
                    X = [X; i_row(sn, p, "both", uo, NaN, NaN, NaN, nnz(B.valid), "")]; %#ok<AGROW>
                    continue
            end
            [u, Xr] = i_read(H, Za(:, ok), sn, p);  X = [X; Xr]; %#ok<AGROW>
            key = struct('periodic', "U", 'surrogate', "U_surrogate", 'total', "U_total");
            U.(key.(p) + sfx) = u;
            if p == "periodic", T = [T; rheome.scale.rows("upstream", "n_frames", nnz(ok), "frames", bandName)]; end %#ok<AGROW>
        end
        if isfield(U, "U" + sfx) && isfield(U, "U_surrogate" + sfx), U.("U_star" + sfx) = U.("U" + sfx) - U.("U_surrogate" + sfx); end
        fprintf('[upstream %s] %s read (%.1f s, chi %.2f, SNR %.1f dB)\n', name, sn, numel(sp)/fs, B.exponent, B.snr);
    end
    f = string(fieldnames(U));
    T = [T; rheome.scale.rows("upstream", f, cellfun(@(k) U.(k), cellstr(f)), "share")];
end

function H = i_geometry(name, S, V1, opts)
% per hemisphere: depth-3 tiles, the tile-mean kernel per component, normals, gradient, own-tile mask; and the
% tile generators of the lagged null (uniform normal current over each tile through the participant's gain)
    K = S.Res.ImagingKernel;  gfile = opts.Gradient;  if isempty(gfile), gfile = getenv('RHEOME_G1'); end
    if isempty(gfile) || ~isfile(gfile)
        error('scale:upstream:gradient', 'no G1 map: pass Gradient= or set RHEOME_G1 (a .mat with L/R .sphere, .g), got "%s"', gfile);
    end
    G1 = load(gfile);
    if ~isfield(S.Sf, 'Sphere') || isempty(S.Sf.Sphere), error('scale:upstream:sphere', '%s: the cortex has no registered sphere', name); end
    gAll = nan(S.nV, 1);  H = struct('Gp', [], 'front', false(1, 0), 'group', zeros(1, 0));
    for h = ["L" "R"]
        Hh = S.B.(h);  gv = double(Hh.gv(:));  Sh = Hh.S;  lbo = Hh.lbo;  a = full(sum(lbo.Mass, 2));
        Tr = rheome.geom.tree(Sh, L=lbo.L, M=lbo.M, MaxDepth=V1.Depth);  G = rheome.geom.tiles(Tr, Sh, V1.Depth);
        Pa = double(G.P) .* a;  Pa = Pa ./ sum(Pa, 1);                     % area-weighted tile mean [nVh x nTile]
        Wc = cell(1, 3);  for c = 1:3, Wc{c} = Pa' * K((gv-1)*3 + c, :); end
        nrm = i_normals(S.Sf, gv, Sh);  N = Pa' * nrm;
        ref = double(G1.(h).sphere);  ref = ref ./ vecnorm(ref, 2, 2);
        us = double(S.Sf.Sphere(gv, :));  us = us ./ vecnorm(us, 2, 2);
        gvx = double(G1.(h).g(:));  gvx = gvx(knnsearch(ref, us));  gAll(gv) = gvx;
        Mem = sparse(size(G.P, 2), S.nV);  Mem(:, gv) = G.P';             %#ok<SPRIX> tile membership on the whole cortex
        own = rheome.flow.crosstalk(K, S.G, Mem);
        Gh = nrm(:,1)' .* S.G(:, (gv-1)*3 + 1) + nrm(:,2)' .* S.G(:, (gv-1)*3 + 2) + nrm(:,3)' .* S.G(:, (gv-1)*3 + 3);
        Gp = Gh * (double(G.P) .* a);  xc = Pa' * double(Sh.Vertices(:, 1));   % SCS: +x anterior
        H.(h) = struct('gv', gv, 'S', Sh, 'P', G.P, 'Pa', Pa, 'W', [Wc{1}; Wc{2}; Wc{3}], 'N', N, 'A', G.A, ...
                       'mask', own >= V1.OwnMin, 'own', own, 'a', a, 'gv_g', gvx, 'tilecheck', NaN);
        H.Gp = [H.Gp Gp];  H.front = [H.front (xc > median(xc))'];  H.group = [H.group repmat(double(h == "R") + 1, 1, numel(xc))];
    end
    H.Pinv = [];
    if opts.Readout == "v2"                                               % fwdLS: the tiles' own forward fields
        LL = H.Gp' * H.Gp;  H.Pinv = (LL + 1e-3 * trace(LL) / size(LL, 1) * eye(size(LL, 1))) \ H.Gp';
    end
    H.null = struct('Gp', H.Gp, 'front', H.front, 'group', H.group, 'keep', false, 'residual', "perchan");
    if opts.Null == "v2"                                                  % depth-6 tile generators, fitted the same way
        H.null = struct('Gp', [], 'front', false(1, 0), 'group', zeros(1, 0), 'keep', true, 'residual', "zerolag");
        for h = ["L" "R"]
            Hh = H.(h);  Tr = rheome.geom.tree(Hh.S, L=S.B.(h).lbo.L, M=S.B.(h).lbo.M, MaxDepth=6);
            G6 = rheome.geom.tiles(Tr, Hh.S, 6);  P6 = double(G6.P) .* Hh.a;  nrm = i_normals(S.Sf, Hh.gv, Hh.S);  gv = Hh.gv;
            Gh = nrm(:,1)' .* S.G(:, (gv-1)*3 + 1) + nrm(:,2)' .* S.G(:, (gv-1)*3 + 2) + nrm(:,3)' .* S.G(:, (gv-1)*3 + 3);
            xc = ((P6 ./ sum(P6, 1))' * double(Hh.S.Vertices(:, 1)))';
            H.null.Gp = [H.null.Gp Gh * P6];  H.null.front = [H.null.front xc > median(xc)];
            H.null.group = [H.null.group repmat(double(h == "R") + 1, 1, numel(xc))];
        end
    end
    sgn = 1;
    if opts.SignCheck                                                     % G1's sign is arbitrary in the normative map
        A = rheome.load.atlas(name, 'Desikan-Killiany');  lab = lower(string(A.Label(:)));
        lo = any(double(A.Membership(contains(lab, ["pericalcarine" "precentral" "postcentral"]), :)), 1)';
        hi = any(double(A.Membership(contains(lab, ["inferiorparietal" "medialorbitofrontal"]), :)), 1)';
        if ~any(lo) || ~any(hi), error('scale:upstream:sign', '%s: no Desikan-Killiany labels for the G1 sign check', name); end
        sgn = sign(mean(gAll(hi), 'omitnan') - mean(gAll(lo), 'omitnan'));
    end
    for h = ["L" "R"]
        Hh = H.(h);  has = isfinite(Hh.gv_g);                            % G1 has no value on the medial wall: the
        g = sgn * (Hh.Pa' * (has .* fillmissing(Hh.gv_g, 'constant', 0)) ./ (Hh.Pa' * has));  H.(h).g = g;   % tile mean of the vertices that have one
        if any(~isfinite(g)), error('scale:upstream:gradient', '%s: a %s tile has no G1 value', name, h); end
        if opts.SignCheck                                                 % the tile holding most of each region
            [~, tl] = max(Hh.P' * (lo(Hh.gv) .* Hh.a));  [~, th] = max(Hh.P' * (hi(Hh.gv) .* Hh.a));
            H.(h).tilecheck = double(g(tl) < g(th));
        end
        H.(h).gradg = sgn * i_grad(Hh.S, Hh.gv_g);
        o = rheome.flow.upstreamshare(ones(numel(g), 1, 2), g, Hh.A, Mask=Hh.mask);   % keptWeight: kernel and g only
        H.(h).keptWeight = o.keptWeight;
    end
end

function B = i_band(F, fs, sp, mask, V1)
% step b on one span: the FFT, the aperiodic fit, the periodic gain on the [1 8] Hz bins, the covariates
    nf = 2^nextpow2(V1.SegS * fs);  N = numel(sp);
    [Pw, fw] = rheome.scale.cleanwelch(F, fs, sp, mask, nf);
    mf = fw >= V1.FitRange(1) & fw <= V1.FitRange(2);
    for l = V1.LineHz * [1 2], mf = mf & ~(fw > l-2 & fw < l+2); end
    C = fft(F(:, sp), [], 2);  f = (0:N-1)' * fs / N;
    bins = find(f >= V1.Band(1) & f <= V1.Band(2));  sc = fs * N / 2;     % |C|^2 = PSD x fs N / 2
    D = rheome.spectral.decompose(C(:, bins), f(bins), struct('fitPower', (Pw(mf, :) * sc)', 'fitF', fw(mf), 'knee', true));
    Pap = D.ap.Pap / sc;  Pk = Pw(mf, :);  inb = fw(mf) >= V1.Band(1) & fw(mf) <= V1.Band(2);
    snr = 10 * log10(sum(max(Pk(inb, :) - Pap(inb, :), 0), 1) ./ sum(Pap(inb, :), 1));
    m = max(1, round(fs / V1.Rate));  idx = 1:m:N;
    B = struct('C', C(:, bins), 'bins', bins, 'N', N, 'hper', D.hper, 'idx', idx, 'rate', fs / m, ...
               'valid', ~mask(sp(idx)), 'exponent', mean(D.ap.exponent), 'offset', mean(D.ap.offset) - log10(sc), ...
               'snr', mean(snr(isfinite(snr))), 'ok', double(all(isfinite(D.ap.exponent) & isfinite(D.ap.offset))));
end

function Za = i_analytic(B, h, ~)
% the analytic slow band: inverse FFT of 2 h C on the positive [1 8] Hz bins, taken at the readout's frames
    Cf = zeros(size(B.C, 1), B.N);  Cf(:, B.bins) = 2 * h .* B.C;
    Za = ifft(Cf, [], 2);  Za = Za(:, B.idx);
end

function Za = i_reanalytic(Y, rate, band)
% a real band-limited signal at the readout's rate made analytic on the same bins
    N = size(Y, 2);  f = (0:N-1)' * rate / N;  C = fft(Y, [], 2);
    k = f >= band(1) & f <= band(2);  Cf = zeros(size(C));  Cf(:, k) = 2 * C(:, k);
    Za = ifft(Cf, [], 2);
end

function [U, X] = i_read(H, Za, sn, p)
% steps c-f: tile values (v1: MNE tile means of the three components; v2: fwdLS), upstreamshare per hemisphere,
% hemispheres pooled by weight
    X = table();  num = 0;  den = 0;  off = 0;
    if ~isempty(H.Pinv), q = H.Pinv * Za; end
    for h = ["L" "R"]
        Hh = H.(h);  n = numel(Hh.g);
        if isempty(H.Pinv)
            o = rheome.flow.upstreamshare(reshape(Hh.W * Za, n, 3, []), Hh.g, Hh.A, Normals=Hh.N, Mask=Hh.mask);
        else                                                              % H.Gp's columns: L tiles, then R tiles
            o = rheome.flow.upstreamshare(reshape(q(off + (1:n), :), n, 1, []), Hh.g, Hh.A, Mask=Hh.mask);
        end
        off = off + n;
        num = num + o.U * o.weight;  den = den + o.weight;
        X = [X; i_row(sn, p, h, o.U, o.weight, o.keptWeight, size(o.edges, 1), size(Za, 2), "")]; %#ok<AGROW>
    end
    U = num / den;
    X = [X; i_row(sn, p, "both", U, den, NaN, NaN, size(Za, 2), "")];
end

function U = i_of(S, H, Za, valid, rate, V1, chunk)
% the comparator: optical flow of the slow-band envelope |J| per hemisphere, share of |v||grad g| a up g
    K = S.Res.ImagingKernel;  num = 0;  den = 0;  nT = size(Za, 2);
    for h = ["L" "R"]
        Hh = H.(h);  gv = Hh.gv;  r = reshape(3*gv(:)' - [2; 1; 0], [], 1);  Kh = K(r, :);  nV = numel(gv);
        gn = vecnorm(Hh.gradg, 2, 2);
        for s = 1:chunk:nT-1
            fr = s:min(s + chunk, nT);  if numel(fr) < 2, continue; end
            J = reshape(Kh * Za(:, fr), 3, nV, []);  A = squeeze(sqrt(sum(abs(J).^2, 1)));
            o = rheome.flow.apparent(A, Hh.S, Rate=rate, Band=V1.Band, Alpha=V1.OfAlpha);
            v = permute(reshape(o.velocity, 3, nV, []), [2 1 3]);              % [nV x 3 x nPairs]
            ok = reshape(valid(fr(1:end-1)) & valid(fr(2:end)), 1, 1, []);
            up = sum(v .* Hh.gradg, 2) > 0;  w = vecnorm(v, 2, 2) .* gn .* Hh.a .* ok;
            num = num + sum(w .* up, 'all');  den = den + sum(w, 'all');
        end
    end
    U = num / den;
end

function gg = i_grad(Sh, g)
% per-vertex surface gradient of g: per-face gradient, area-averaged to the vertices
    fg = rheome.operators.face_gradient(Sh.Vertices, Sh.Faces);
    gg = fg.W * [fg.Gx * g, fg.Gy * g, fg.Gz * g];
    gg(~isfinite(gg)) = 0;                                                % next to the medial wall (no G1): no weight
end

function N = i_normals(Sf, gv, Sh)
    if isfield(Sf, 'VertNormals') && ~isempty(Sf.VertNormals), N = double(Sf.VertNormals(gv, :));
    else, N = double(Sh.VertNormals); end
    N = N ./ vecnorm(N, 2, 2);
end

function R = i_row(sn, p, h, U, w, kw, ne, nf, why)
    R = table(string(sn), string(p), string(h), U, w, kw, ne, nf, string(why), 'VariableNames', ...
              {'span','pass','hemi','U','weight','kept_weight','n_edges','n_frames','reason'});
end

% Author: Diellor Basha, 2026
