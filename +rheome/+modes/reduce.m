function [C, Ts] = reduce(fitFiles, centroidsFile, outDir, opts)
% MODES.REDUCE  The keep rule per mode x band over a cohort's back-fits: Level 1, phase-layout labels, Level 2.
%
%   [C, Ts] = rheome.modes.reduce(fitFiles, centroidsFile, outDir)
%   [C, Ts] = rheome.modes.reduce("fits.txt", "centroids.mat", outDir, Cohort="omega-H", Alpha=0.05)
%
% fitFiles: the participants' modesfit.csv (rheome.modes.fit via rheome.scale.run), as paths or a text file
% of paths. Protocol section 5, fixed before the data: every test a two-sided one-sample Wilcoxon (signrank)
% on the participant-wise statistic minus its null, Holm across modes x bands within this cohort, Alpha.
%   Level 1 (real spatial structure), all of:
%     (i)   r_halves - r_n1_sources > 0 (the participant's own mode matches across halves better than it
%           matches its leakage null) AND the centroid's |r| to every N1 centroid of the discovery < N1R (0.7)
%     (ii)  N2 spin: r_spin_obs - r_spin_mean > 0
%     (iii) N4: occupancy on the periodic part - occupancy on the aperiodic part > 0
%     (iv)  split-half ICC(2,1) of the mode's gev >= 0.5 with the lower bound of its 95 % CI > 0.3 (percentile
%           bootstrap over participants, NBoot)
%   Labels: standing (median pg < StandingRad); travelling or rotating (r_n3_sim - r_n3 > N3Margin under Holm --
%     the mode is further from its best two-phase source than that source reads when simulated in the
%     participant's own data; see rheome.modes.fit on why not pg - pg_n3 -- AND the sign
%     of div_index -- travelling -- or of winding_net -- rotating -- shared by >= SignShare of participants);
%     lagged-pair (neither). The label is read only for Level-1 modes and says nothing about speed.
%   Level 2 (dynamics beyond second order): dwell_median_ms (the mode's) and trans_entropy (overall, under
%     mode 0) minus the N5 surrogate mean; plus the share of participants with N5 p significant at FDR q
%     (Benjamini-Hochberg over participants) against 5 % (binomial test).
% Characterisation (protocol 6): participant medians (IQR) of extent, Dirichlet wavelength, occupancy,
% occurrence, dwell, gev; the mean direction vector, its anatomical axis label and concentration; speed always
% "unresolved" at depths below 7; the discovery's r_dmd, r_cpca, robust.
%
% Writes outDir/modes_catalogue.csv (one row per band x mode), modes_tests.csv (one row per test),
% reduce_provenance.json.
%
% See also: rheome.modes.fit, rheome.modes.discover
%
% Author: Diellor Basha, 2026

    arguments
        fitFiles string
        centroidsFile (1,:) char
        outDir (1,:) char
        opts.Cohort (1,1) string = ""
        opts.Alpha (1,1) double = 0.05
        opts.N1R (1,1) double = 0.7
        opts.StandingRad (1,1) double = pi/8
        opts.SignShare (1,1) double = 0.8
        opts.Fdr (1,1) double = 0.05
        opts.NBoot (1,1) double = 2000
        opts.Seed (1,1) double = 1
        opts.N3Margin (1,1) double = 0      % 0 = the protocol as fixed; a floor on the median N3 excess
    end
    if isscalar(fitFiles) && endsWith(fitFiles, ".txt"), fitFiles = strip(splitlines(strtrim(fileread(fitFiles)))); fitFiles(fitFiles == "") = []; end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    L = load(centroidsFile, 'M');  M = L.M;
    if ~strcmp(rheome.modes.hash(M.bands.C), M.sha256), error('modes:reduce:hash', 'centroids.mat does not match its own hash'); end
    D = table();
    for f = fitFiles(:)'
        t = readtable(f, 'TextType', 'string', 'Delimiter', ',');
        if ~ismember("subject", t.Properties.VariableNames), t.subject = repmat(string(f), height(t), 1); end
        D = [D; t(:, {'subject','band','source','half','method','mode','metric','value'})]; %#ok<AGROW>
    end
    subs = unique(D.subject);  nS = numel(subs);  rs = RandStream('mt19937ar', 'Seed', opts.Seed);
    get = @(b, src, half, meth, k, metric) i_vec(D, subs, b, src, half, meth, k, metric);
    Ts = table();  C = table();
    for b = 1:numel(M.bands)
        bn = M.bands(b).band;  K = M.bands(b).K;
        for k = 1:K
            dd = @(m) get(bn, "data", "all", "backfit", k, m);
            row = table(opts.Cohort, bn, k, K, nS, 'VariableNames', {'cohort','band','mode','K','n'});
            % Level 1
            s1 = dd("r_halves") - dd("r_n1_sources");
            s2 = dd("r_spin_obs") - dd("r_spin_mean");
            s3 = dd("occupancy") - get(bn, "n4_aperiodic", "all", "backfit", k, "occupancy");
            sN3 = dd("r_n3_sim") - dd("r_n3");
            sL2 = dd("dwell_median_ms") - get(bn, "n5_phaserand", "all", "backfit", k, "dwell_median_ms");
            for tst = {"i_n1", s1; "ii_spin", s2; "iii_n4", s3; "n3_lagged", sN3; "l2_dwell", sL2}'
                [p, md] = i_wilcox(tst{2});
                Ts = [Ts; table(opts.Cohort, bn, k, string(tst{1}), md, p, nnz(isfinite(tst{2})), 'VariableNames', ...
                      {'cohort','band','mode','test','median','p','n'})]; %#ok<AGROW>
            end
            r1 = max([0; i_nullr(M.bands(b), k, "n1_sources"); i_nullr(M.bands(b), k, "n1_emptyroom")]);
            g1 = get(bn, "data", "h1", "backfit", k, "gev");  g2 = get(bn, "data", "h2", "backfit", k, "gev");
            [icc, lo, hi] = i_icc(g1, g2, opts.NBoot, rs);
            pgm = median(dd("pg"), 'omitnan');
            dv = dd("div_index");  wn = dd("winding_net");
            shDiv = max(mean(dv > 0, 'omitnan'), mean(dv < 0, 'omitnan'));  shRot = max(mean(wn > 0, 'omitnan'), mean(wn < 0, 'omitnan'));
            pL2 = get(bn, "n5_phaserand", "all", "backfit", k, "dwell_median_ms_p");
            [shareL2, pBin] = i_share(pL2, opts.Fdr);
            v = [median(dd("dir_x"), 'omitnan') median(dd("dir_y"), 'omitnan') median(dd("dir_z"), 'omitnan')];
            row = [row, table(r1, icc, lo, hi, pgm, shDiv, shRot, shareL2, pBin, i_axis(v), norm(v), ...
                   median(dd("extent_mm2"), 'omitnan'), median(dd("extent_radius_mm"), 'omitnan'), ...
                   median(dd("dirichlet_wavelength_mm"), 'omitnan'), median(dd("tile_spacing_mm"), 'omitnan'), ...
                   median(dd("occupancy"), 'omitnan'), median(dd("occurrence_per_s"), 'omitnan'), ...
                   median(dd("dwell_median_ms"), 'omitnan'), median(dd("gev"), 'omitnan'), median(dd("r_centroid"), 'omitnan'), ...
                   median(dd("div_index"), 'omitnan'), median(dd("winding_net"), 'omitnan'), ...
                   "unresolved", M.bands(b).r_dmd(k), M.bands(b).r_cpca(k), ...
                   'VariableNames', {'r_n1_centroid_max','icc_gev','icc_lo','icc_hi','pg_median','div_sign_share', ...
                   'winding_sign_share','l2_share_sig','l2_share_p','direction','direction_concentration', ...
                   'extent_mm2','extent_radius_mm','dirichlet_wavelength_mm','tile_spacing_mm','occupancy', ...
                   'occurrence_per_s','dwell_median_ms','gev','r_centroid','div_index','winding_net','speed', ...
                   'r_dmd','r_cpca'})]; %#ok<AGROW>
            C = [C; row]; %#ok<AGROW>
        end
        % overall transition entropy against N5 (mode 0)
        sH = get(bn, "data", "all", "backfit", 0, "trans_entropy") - get(bn, "n5_phaserand", "all", "backfit", 0, "trans_entropy");
        [p, md] = i_wilcox(sH);
        Ts = [Ts; table(opts.Cohort, bn, 0, "l2_entropy", md, p, nnz(isfinite(sH)), 'VariableNames', {'cohort','band','mode','test','median','p','n'})]; %#ok<AGROW>
    end
    % Holm within each test family across modes x bands, then the rule
    Ts.p_holm = nan(height(Ts), 1);
    for tst = unique(Ts.test)', j = Ts.test == tst;  Ts.p_holm(j) = i_holm(Ts.p(j)); end
    pass = @(tst, b, k) i_pass(Ts, tst, b, k, opts.Alpha);
    for r = 1:height(C)
        b = C.band(r);  k = C.mode(r);
        C.pass_i(r) = pass("i_n1", b, k) && C.r_n1_centroid_max(r) < opts.N1R;
        C.pass_ii(r) = pass("ii_spin", b, k);  C.pass_iii(r) = pass("iii_n4", b, k);
        C.pass_iv(r) = C.icc_gev(r) >= 0.5 && C.icc_lo(r) > 0.3;
        C.level1(r) = C.pass_i(r) && C.pass_ii(r) && C.pass_iii(r) && C.pass_iv(r);
        n3 = pass("n3_lagged", b, k) && Ts.median(Ts.test == "n3_lagged" & Ts.band == b & Ts.mode == k) > opts.N3Margin;
        if C.pg_median(r) < opts.StandingRad, lab = "standing";
        elseif n3 && C.winding_sign_share(r) >= opts.SignShare, lab = "rotating";
        elseif n3 && C.div_sign_share(r) >= opts.SignShare, lab = "travelling";
        else, lab = "lagged-pair"; end
        C.label(r) = lab;  C.level2(r) = pass("l2_dwell", b, k) || pass("l2_entropy", b, 0);
    end
    rb = arrayfun(@(b, k) M.bands(string(b) == [M.bands.band]).r_dmd(k) >= 0.7 || M.bands(string(b) == [M.bands.band]).r_cpca(k) >= 0.7, C.band, C.mode);
    C.robust = rb;
    C.mixture_share = arrayfun(@(b, k) M.bands(string(b) == [M.bands.band]).mixture_share(k), C.band, C.mode);
    C.centroids_sha256 = repmat(string(M.sha256), height(C), 1);
    writetable(C, fullfile(outDir, 'modes_catalogue.csv'));  writetable(Ts, fullfile(outDir, 'modes_tests.csv'));
    pv = struct('cohort', opts.Cohort, 'participants', nS, 'sha256', M.sha256, 'options', opts, 'matlab', version);
    fid = fopen(fullfile(outDir, 'reduce_provenance.json'), 'w');  fprintf(fid, '%s', jsonencode(pv, PrettyPrint=true));  fclose(fid);
end

function v = i_vec(D, subs, b, src, half, meth, k, metric)
    j = D.band == b & D.source == src & D.half == half & D.method == meth & D.mode == k & D.metric == metric;
    v = nan(numel(subs), 1);  [ok, w] = ismember(D.subject(j), subs);  x = D.value(j);  v(w(ok)) = x(ok);
end

function [p, md] = i_wilcox(x)
    x = x(isfinite(x));  md = median(x);  p = NaN;
    if numel(x) >= 2 && any(x ~= 0), p = signrank(x); end
end

function ok = i_pass(Ts, tst, b, k, a)
    j = Ts.test == tst & Ts.band == b & Ts.mode == k;
    ok = any(j) && Ts.p_holm(find(j, 1)) < a && Ts.median(find(j, 1)) > 0;
end

function q = i_holm(p)
    q = nan(size(p));  ok = isfinite(p);  [ps, o] = sort(p(ok));  m = numel(ps);
    adj = min(1, cummax((m - (1:m)' + 1) .* ps));  t = find(ok);  q(t(o)) = adj;
end

function [icc, lo, hi] = i_icc(a, b, nb, rs)
% ICC(2,1) (two-way random, absolute agreement, single measure; Shrout & Fleiss) and a percentile bootstrap CI
    ok = isfinite(a) & isfinite(b);  Y = [a(ok) b(ok)];  icc = NaN;  lo = NaN;  hi = NaN;
    if size(Y, 1) < 3, return; end
    f = @(Y) i_icc21(Y);  icc = f(Y);  n = size(Y, 1);  bs = zeros(nb, 1);
    for s = 1:nb, bs(s) = f(Y(randi(rs, n, n, 1), :)); end
    lo = prctile(bs, 2.5);  hi = prctile(bs, 97.5);
end

function r = i_icc21(Y)
    [n, k] = size(Y);  g = mean(Y(:));
    msr = k * sum((mean(Y, 2) - g).^2) / (n - 1);  msc = n * sum((mean(Y, 1) - g).^2) / (k - 1);
    sse = sum((Y - mean(Y, 2) - mean(Y, 1) + g).^2, 'all');  mse = sse / ((n - 1) * (k - 1));
    r = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n);
end

function [sh, p] = i_share(pv, q)
% share of participants significant at FDR q (Benjamini-Hochberg over participants), binomial test vs 5 %
    pv = pv(isfinite(pv));  m = numel(pv);  sh = NaN;  p = NaN;
    if m == 0, return; end
    ps = sort(pv);  k = find(ps <= (1:m)' / m * q, 1, 'last');  if isempty(k), k = 0; end
    sh = k / m;  p = 1 - binocdf(k - 1, m, 0.05);
end

function r = i_nullr(B, k, src)
    r = [];  if isfield(B.nullC, src) && ~isempty(B.nullC.(src)), r = max(abs(B.C(:, k)' * B.nullC.(src))); end
end

function s = i_axis(v)
% the dominant anatomical axis of a direction in template (SCS) coordinates: x anterior, y left, z up
    if any(~isfinite(v)) || norm(v) == 0, s = "none";  return; end
    lab = ["A->P" "P->A"; "L->R" "R->L"; "D->V" "V->D"];  [~, a] = max(abs(v));  s = lab(a, 1 + (v(a) > 0));
end

% Author: Diellor Basha, 2026
