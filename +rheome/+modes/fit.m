function [T, X] = fit(name, S, M, opts)
% MODES.FIT  Back-fit the frozen centroids to one participant: states, physical features and the per-person nulls.
%
%   [T, X] = rheome.modes.fit(name, S, M)                     % M = the struct in centroids.mat
%   [T, X] = rheome.modes.fit(name, S, M, Sha256="...", NSpin=1000, NSurr=99)
%
% Protocol sections 4.2, 5 and 6, per participant. The fields are recomputed by rheome.modes.fields with the
% variant and seed of the discovery (deterministic: the same commit gives the same series; the series hash
% is in the output) and the centroids are refused unless rheome.modes.hash(M.bands.C) equals M.sha256 (and
% Sha256 when given).
%
% Per band, on every valid sample (snapshots normalised; label = argmax_k |<z, c_k>|^2):
%   backfit (data, h1, h2, n1_sources, n1_emptyroom, n4_aperiodic): per mode occupancy, occurrence_per_s,
%     dwell_median_ms (runs cut by an invalid sample or a recording end are not dwell times), gev (share of
%     the participant's power the mode explains, sum |<z, c_k>|^2 over its samples / sum |z|^2), r_centroid
%     (|<m_k, c_k>|, m_k the participant's own mode: principal eigenvector of its samples' scatter); overall
%     gev_total and trans_entropy (bits; entropy of the off-diagonal transitions, weighted by visits)
%   hmm (data): a K-state HMM over the SAME centroids, complex Watson emissions exp(kappa |<z, c_k>|^2),
%     kappa = (d-1)/(1 - mean max_k |<z, c_k>|^2) (its large-kappa moment), transitions by Baum-Welch
%     (10 iterations), Viterbi labels -- the same temporal features (sensitivity, protocol 4.2)
%   Level-1 per person: r_halves = |<m_k(h1), m_k(h2)>|, r_n1 = |<m_k, m_k(n1_sources)>| (and empty room);
%     occupancy minus N4 occupancy is the reduction's
%   N2 spin: r_obs = |<m_k, c_k>| against NSpin rotations of m_k on the template sphere (one random rotation
%     per spin, mirrored for the right hemisphere; each tile takes the value of the tile nearest its rotated
%     centre): r_spin_mean, p_spin = (1 + #{r_spin >= r_obs}) / (1 + NSpin)
%   N3 lagged pair: rheome.modes.lagged -- the best source with exactly TWO phases (two disjoint, extended,
%     nonnegative generators, any ratio and lag) seen through the participant's gain and the pipeline (R, the
%     tile resolution): r_n3 = |r| of m_k to it. Then the null SIMULATED in the participant's own data: in every
%     sample of state k the mode's component is replaced by that fit (noise and the rest kept), the mode is
%     re-estimated and refitted: r_n3_sim. The reduction tests r_n3_sim - r_n3 > 0 (the mode is further from
%     any lagged pair than a lagged pair in this person's data would read). pg_n3 and the lag are reported.
%     ⚠ Not the protocol's "phase gradient exceeds N3": two lagged generators make a phase STEP whose mean
%     |phase difference| over edges exceeds a smooth gradient's (a planted A->P wave read pg 1.57 against its
%     pair's 2.01). ⚠ The separation is small at these scales: through one MEG participant's head at depth 3 a 200 mm wave
%     reads r_n3 0.96, planted pairs 0.993-0.997, a rotation 0.994
%   N3 calibration (protocol section 13, method "n3cal", mode 0): on this head, a planted 200 mm A->P plane wave
%     (source planted_wave200_ap) and two planted lagged pairs (planted_pair_fo, planted_pair_lateral) through
%     W * gain replace the participant's leading component in every valid sample (noise and the rest kept); the
%     re-estimated mode goes through the same N3: r_n3, r_n3_sim, n3_excess = r_n3_sim - r_n3, r_plant. The
%     frozen N3Margin is the midpoint between the largest pair excess and the smallest wave excess (Science)
%   N5 (Level 2): NSurr multivariate phase-randomised surrogates of the periodic tile series through the same
%     band and back-fit: dwell_median_ms and trans_entropy as surrogate mean, sd and two-sided p
%   physical features of m_k (protocol 6), on the participant's tile geometry:
%     extent_mm2 / extent_radius_mm  area of the tiles at >= half the maximum amplitude, equivalent radius
%     dirichlet_wavelength_mm        2 pi / sqrt(lambda), lambda the tile-graph Rayleigh quotient
%                                    2 mean_E |m_i - m_j|^2/d_ij^2 / mean_E (|m_i|^2+|m_j|^2)/2 (k^2 for a plane wave)
%     tile_spacing_mm                median centre distance of adjacent tiles: below ~2x this nothing is resolved
%     pg                             phase-gradient statistic: amplitude-weighted mean |wrap(phi_i - phi_j)| over
%                                    edges between tiles >= half max; NaN with fewer than 3 such tiles
%     dir_x dir_y dir_z, dir_concentration, dir_azimuth_deg   direction of travel (-grad phi; phase lead to
%                                    lag) as the amplitude-weighted mean unit vector in template coordinates (SCS:
%                                    x anterior, y left, z up) and its sphere-frame azimuth (0 = north, 90 = u x north)
%     phase_wavelength_mm, speed_m_s, speed_resolved    f_centre x 2 pi/|grad phi|; resolved only at depth >=
%                                    ResolvableDepth (MS1: speed is recovered at depths 7-8 only), so at depths
%                                    2-4 it is ALWAYS 0: "direction only"
%     div_index (rad/mm^2, + = diverging from a phase lead), winding_net, winding_count (signed and unsigned
%                                    phase windings around tile triangles, outward orientation). ⚠ Not s^-1: a
%                                    flow rate needs the speed, which these depths do not resolve
%
% JOINT ATOMS (M.kind "joint": centroids of rheome.modes.discover on rheome.modes.jointfields; protocol section 14):
% the joint fields are recomputed instead (same seed); F (f_fieldspread) takes N1's place (r_f_fieldspread);
% N2 spins within each depth (the same rotation for every depth); physical features and N3 are read on the
% atom's dominant depth (the depth holding most of its energy; dominant_depth), N3 reported and never gated;
% per atom also scale_mm (energy-weighted geometric mean of the depths' scales), t_bar (energy-weighted T(s) of
% this participant, rheome.modes.scaletransfer) and scale_edge_mm; N5 surrogates through rheome.modes.octave;
% no N3 calibration plants.
%
% Returns T (rheome.scale.rows "modesfit": per band gev_total, n_valid, kappa, ...) and X with tables
%   modesfit  (band, source, half, method, mode, metric, value, unit; mode 0 = overall)
%   modesmaps (band, source, half, mode, hemi, tile, amp, phase_rad)
%   modestrans (band, source, method, from, to, p)
%
% See also: rheome.modes.discover, rheome.modes.fields, rheome.modes.reduce
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S
        M (1,1) struct
        opts.Sha256 (1,:) char = ''
        opts.NSpin (1,1) double = 1000
        opts.NSurr (1,1) double = 99
        opts.Seed (1,1) double = 1
        opts.ResolvableDepth (1,1) double = 7
        opts.HmmIter (1,1) double = 10
    end
    h0 = rheome.modes.hash(M.bands.C);
    if ~strcmp(h0, M.sha256) || (~isempty(opts.Sha256) && ~strcmpi(opts.Sha256, h0))
        error('modes:fit:hash', 'centroids hash %s does not match the frozen %s %s', h0, M.sha256, opts.Sha256);
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    vs = M.variantSpec;  jt = isfield(M, 'kind') && M.kind == "joint";
    if jt
        o = M.fieldOptions;
        [Z, Y] = rheome.modes.jointfields(name, S, Seed=opts.Seed, Depths=o.Depths, Octaves=o.Octaves, GenDepth=o.GenDepth, ...
                                          WorkRate=o.WorkRate, FitRange=o.FitRange);  Pn = strings(0);  Pm = [];
    else
        [Z, Y, Op] = rheome.modes.fields(name, S, CutoffMM=vs.cutoff_mm, Depth=vs.depth, SensCutoffMM=[], SensDepths=[], Seed=opts.Seed);
        [Pn, Pm] = i_planted(S, Op{1}.W);  clear Op
    end
    tl = Z.tiles{1};
    if height(tl) ~= height(M.tiles) || any(tl.tile ~= M.tiles.tile) || any(tl.hemi ~= M.tiles.hemi)
        error('modes:fit:tiles', '%s: tiles differ from the frozen centroids''', name);
    end
    G = i_geom(tl, Z.A);  Gd = {};  dep = 0;
    if jt
        dep = unique(tl.depth);  Gd = arrayfun(@(d) i_geom(tl(tl.depth == d, :), Z.A(tl.depth == d, tl.depth == d)), dep, 'uni', 0);
        Tq = Z.transfer.table.T(arrayfun(@(d) find(Z.transfer.table.depth == d), tl.depth));
    end
    rs = RandStream('mt19937ar', 'Seed', opts.Seed + 11);
    F = cell(0, 8);  Mp = table();  Tr = table();  T = rheome.scale.rows("modesfit", strings(0,1), [], "");

    function put(b, src, half, meth, k, metric, v, unit)
        F(end+1, :) = {b, src, half, meth, k, metric, v, unit};
    end
    function maps(b, src, half, Mk)
        for jj = 1:size(Mk, 2)
            Mp = [Mp; table(repmat(b, height(tl), 1), repmat(src, height(tl), 1), repmat(half, height(tl), 1), ...
                  repmat(jj, height(tl), 1), tl.hemi, tl.tile, abs(Mk(:, jj)), angle(Mk(:, jj)), 'VariableNames', ...
                  {'band','source','half','mode','hemi','tile','amp','phase_rad'})];
        end
    end
    function feats(b, src, half, meth, R, K)
        for jj = 1:K
            put(b, src, half, meth, jj, "occupancy", R.occ(jj), "fraction");
            put(b, src, half, meth, jj, "occurrence_per_s", R.rate(jj), "1/s");
            put(b, src, half, meth, jj, "dwell_median_ms", R.dwell(jj), "ms");
            put(b, src, half, meth, jj, "gev", R.gev(jj), "fraction");
            put(b, src, half, meth, jj, "r_centroid", R.rc(jj), "r");
        end
        put(b, src, half, meth, 0, "gev_total", sum(R.gev), "fraction");
        put(b, src, half, meth, 0, "trans_entropy", R.H, "bits");
        put(b, src, half, meth, 0, "dwell_median_ms", R.dwellAll, "ms");
        [ii, jj] = ndgrid(1:K, 1:K);
        Tr = [Tr; table(repmat(b, K*K, 1), repmat(src, K*K, 1), repmat(meth + "_" + half, K*K, 1), ii(:), jj(:), R.P(:), ...
              'VariableNames', {'band','source','method','from','to','p'})];
    end

    for b = 1:numel(M.bands)
        bn = M.bands(b).band;  C = M.bands(b).C;  K = size(C, 2);
        bi = find(Z.bands.band == bn);  rate = Z.bands.rate(bi);  edges = [Z.bands.lo(bi) Z.bands.hi(bi)];
        valid = Y.(sprintf('valid_%s', bn));  z = double(Y.(sprintf('z_data_%s_v1', bn)));
        fc = mean(edges);
        % data, its halves, and the HMM
        Rd = i_backfit(z, valid, C, rate);  feats(bn, "data", "all", "backfit", Rd, K);  maps(bn, "data", "all", Rd.m);
        iv = find(valid);  hv = iv(ceil(numel(iv)/2));
        for hh = 1:2
            vh = valid;  if hh == 1, vh(hv+1:end) = false; else, vh(1:hv) = false; end
            Rh{hh} = i_backfit(z, vh, C, rate);  feats(bn, "data", "h" + hh, "backfit", Rh{hh}, K);  maps(bn, "data", "h" + hh, Rh{hh}.m); %#ok<AGROW>
        end
        [lab, kap] = i_hmm(z, valid, C, opts.HmmIter);
        feats(bn, "data", "all", "hmm", i_backfit(z, valid, C, rate, lab), K);
        if isfield(Z, 'periodic'), put(bn, "data", "all", "backfit", 0, "periodic_fraction_median", median(Z.periodic.periodic_fraction(Z.periodic.band == bn)), "fraction"); end
        T = [T; rheome.scale.rows("modesfit", ["gev_total" "n_valid" "kappa" "K" "band_lo_hz" "band_hi_hz" "rate_hz"], ...
             [sum(Rd.gev) nnz(valid) kap K edges rate], ["fraction" "samples" "1" "modes" "Hz" "Hz" "Hz"], bn)]; %#ok<AGROW>
        for kk = 1:K
            put(bn, "data", "all", "backfit", kk, "r_halves", abs(Rh{1}.m(:, kk)' * Rh{2}.m(:, kk)), "r");
        end
        % the instrument and background nulls through the same back-fit
        for src = Z.sources(ismember(Z.sources, ["n1_sources" "n1_emptyroom" "f_fieldspread" "n4_aperiodic"]))
            Rn = i_backfit(double(Y.(sprintf('z_%s_%s_v1', src, bn))), valid, C, rate);
            feats(bn, src, "all", "backfit", Rn, K);  maps(bn, src, "all", Rn.m);
            if startsWith(src, "n1") || src == "f_fieldspread"
                for kk = 1:K, put(bn, "data", "all", "backfit", kk, "r_" + src, abs(Rd.m(:, kk)' * Rn.m(:, kk)), "r"); end
            end
        end
        % per-mode: physical features, N2 spin, N3 lagged pair
        for kk = 1:K
            m = Rd.m(:, kk);
            if jt                                                       % the dominant depth carries geometry and N3
                e = arrayfun(@(d) sum(abs(m(tl.depth == d)).^2), dep);  [~, id] = max(e);  in = tl.depth == dep(id);
                Ph = i_phys(m(in), Gd{id}, fc, dep(id), opts.ResolvableDepth);
                Ph.dominant_depth = dep(id);  Ph.scale_mm = exp(sum(abs(m).^2 .* log(tl.scale_mm)) / sum(abs(m).^2));
                Ph.t_bar = sum(abs(m).^2 .* Tq) / sum(abs(m).^2);  Ph.scale_edge_mm = Z.transfer.edge_mm;
            else
                Ph = i_phys(m, G, fc, vs.depth, opts.ResolvableDepth);  in = true(size(m));  id = 1;
            end
            for f = string(fieldnames(Ph))', put(bn, "data", "all", "backfit", kk, f, Ph.(f), i_unit(f)); end
            [ro, rsm, ps] = i_spin(m, C(:, kk), G, opts.NSpin, rs);
            put(bn, "data", "all", "backfit", kk, "r_spin_obs", ro, "r");  put(bn, "data", "all", "backfit", kk, "r_spin_mean", rsm, "r");
            put(bn, "data", "all", "backfit", kk, "p_spin", ps, "p");
            if jt, L3 = i_lagged(m(in), Z.R{id}, Gd{id}, z(in, Rd.lab == kk)); else, L3 = i_lagged(m, Z.R, G, z(:, Rd.lab == kk)); end
            for f = string(fieldnames(L3))', put(bn, "data", "all", "backfit", kk, f, L3.(f), i_unit(f)); end
        end
        % N3 calibration (protocol section 13): planted patterns in this participant's own band data
        for pp = 1:numel(Pn)
            Lc = i_n3plant(Pm(:, pp), Z.R, G, z(:, valid));
            for f = string(fieldnames(Lc))', put(bn, Pn(pp), "all", "n3cal", 0, f, Lc.(f), i_unit(f)); end
        end
        % N5: dwell and transition entropy against phase-randomised surrogates (Level 2)
        sd = nan(opts.NSurr, K + 2);
        for s = 1:opts.NSurr
            if jt, zs = rheome.modes.octave(rheome.modes.phaserand(double(Y.yper), rs), Z.workRate, edges, Y.mask);
            else, zs = rheome.modes.analytic(rheome.modes.phaserand(double(Y.yper), rs), Z.workRate, edges, Y.mask); end
            Rs = i_backfit(double(zs), valid, C, rate);  sd(s, :) = [Rs.dwell(:)' Rs.dwellAll Rs.H];
        end
        dv = [Rd.dwell(:)' Rd.dwellAll Rd.H];  mu = mean(sd, 1, 'omitnan');  sg = std(sd, 0, 1, 'omitnan');
        p = (1 + sum(abs(sd - mu) >= abs(dv - mu), 1)) ./ (1 + sum(isfinite(sd), 1));
        for kk = 1:K + 2
            if kk <= K, k = kk;  mt = "dwell_median_ms";  u = "ms";
            elseif kk == K + 1, k = 0;  mt = "dwell_median_ms";  u = "ms";
            else, k = 0;  mt = "trans_entropy";  u = "bits"; end
            put(bn, "n5_phaserand", "all", "backfit", k, mt, mu(kk), u);
            put(bn, "n5_phaserand", "all", "backfit", k, mt + "_sd", sg(kk), u);
            put(bn, "n5_phaserand", "all", "backfit", k, mt + "_p", p(kk), "p");
        end
        fprintf('[modes.fit %s] %s: K=%d gev %.2f, dwell %.0f ms (N5 %.0f), entropy %.2f (N5 %.2f)\n', name, bn, K, ...
                sum(Rd.gev), Rd.dwellAll, mu(K+1), Rd.H, mu(K+2));
    end
    X.modesfit = cell2table(F, 'VariableNames', {'band','source','half','method','mode','metric','value','unit'});
    X.modesmaps = Mp;  X.modestrans = Tr;
    X.modesseries = Z.series(:, {'source','band','variant','n_valid','sha256'});
    T = [T; rheome.scale.rows("modesfit", "centroids_sha256_ok", 1, "flag", string(h0))];
end

%% ---------- back-fit
function R = i_backfit(z, valid, C, rate, lab)
    K = size(C, 2);  n = vecnorm(z, 2, 1);  Pw = abs(C' * z).^2;
    if nargin < 5, [~, lab] = max(Pw, [], 1); end
    lab(~valid) = 0;  nv = nnz(valid);  tot = sum(n(valid).^2);  R.lab = lab;
    R.occ = accumarray(lab(valid)', 1, [K 1])' / nv;
    R.gev = arrayfun(@(k) sum(Pw(k, lab == k)), 1:K) / tot;
    % runs: a run touching an invalid sample or an end is not a dwell time
    d = diff([-1 lab -1]);  st = find(d(1:end-1) ~= 0);  en = [st(2:end) - 1, numel(lab)];
    rl = lab(st);  len = en - st + 1;
    cut = [true, lab(max(st(2:end) - 1, 1)) == 0] | [lab(min(en(1:end-1) + 1, numel(lab))) == 0, true];
    ok = rl > 0;
    R.rate = accumarray(rl(ok)', 1, [K 1])' / (nv / rate);
    full_ = ok & ~cut;
    R.dwell = arrayfun(@(k) median(len(full_ & rl == k)) / rate * 1e3, 1:K);
    R.dwellAll = median(len(full_)) / rate * 1e3;
    a = lab(1:end-1);  bb = lab(2:end);  tt = a > 0 & bb > 0 & a ~= bb;
    Cn = accumarray([a(tt)' bb(tt)'], 1, [K K]);  R.P = Cn ./ max(sum(Cn, 2), 1);
    pi_ = sum(Cn, 2) / max(sum(Cn, 'all'), 1);  lp = R.P .* log2(R.P);  lp(R.P == 0) = 0;
    R.H = -sum(pi_ .* sum(lp, 2));
    R.m = zeros(size(C));  R.rc = nan(1, K);
    for k = 1:K
        zk = z(:, lab == k);
        if size(zk, 2) < 2, continue; end
        [V, L] = eig(zk * zk', 'vector');  [~, j] = max(real(L));  m = V(:, j);
        m = m * exp(-1i * angle(C(:, k)' * m));  R.m(:, k) = m;  R.rc(k) = abs(C(:, k)' * m);
    end
end

function [lab, kap] = i_hmm(z, valid, C, nit)
% K-state HMM over fixed centroids, complex Watson emissions, Baum-Welch on the transitions, Viterbi
    K = size(C, 2);  d = size(C, 1);  x = z(:, valid) ./ vecnorm(z(:, valid));  Pw = abs(C' * x).^2;
    kap = (d - 1) / max(1e-6, 1 - mean(max(Pw, [], 1)));
    E = exp(kap * (Pw - max(Pw, [], 1)));  Tn = size(E, 2);
    A = 0.9 * eye(K) + 0.1 / max(K - 1, 1) * ~eye(K);  if K == 1, A = 1; end
    p0 = ones(K, 1) / K;
    for it = 1:nit
        al = zeros(K, Tn);  c = zeros(1, Tn);  a = p0 .* E(:, 1);  c(1) = sum(a);  al(:, 1) = a / c(1);
        for t = 2:Tn, a = (A' * al(:, t-1)) .* E(:, t);  c(t) = sum(a);  al(:, t) = a / c(t); end
        be = ones(K, Tn);
        for t = Tn-1:-1:1, be(:, t) = A * (E(:, t+1) .* be(:, t+1)) / c(t+1); end
        Xi = (al(:, 1:end-1) * (E(:, 2:end) .* be(:, 2:end) ./ c(2:end))') .* A;
        A = Xi ./ sum(Xi, 2);  g = al(:, 1) .* be(:, 1);  p0 = g / sum(g);
    end
    lA = log(A);  dl = log(p0) + kap * Pw(:, 1);  bp = zeros(K, Tn, 'uint8');
    for t = 2:Tn, [dl, bp(:, t)] = max(dl + lA, [], 1);  dl = dl' + kap * Pw(:, t); end
    q = zeros(1, Tn);  [~, q(Tn)] = max(dl);
    for t = Tn:-1:2, q(t-1) = bp(q(t), t); end
    lab = zeros(1, numel(valid));  lab(valid) = q;
end

%% ---------- geometry, physical features, nulls
function G = i_geom(tl, A)
    G.A = A;  G.area = tl.area_mm2;  G.x = [tl.x_mm tl.y_mm tl.z_mm];  G.t = [tl.tx tl.ty tl.tz];
    G.u = [tl.u1 tl.u2 tl.u3];  G.u = G.u ./ vecnorm(G.u, 2, 2);  G.hemi = tl.hemi;
    gk = tl.hemi;  if ismember('depth', tl.Properties.VariableNames), gk = gk + "_" + tl.depth; end
    [~, ~, gi] = unique(gk);  G.groups = accumarray(gi, (1:numel(gi))', [], @(x) {x});   % spin within hemisphere x depth
    [i, j] = find(triu(A));  G.e = [i j];  G.d = vecnorm(G.x(i, :) - G.x(j, :), 2, 2);
    ax = [0 0 1];  e1 = ax - (G.u * ax') .* G.u;  G.e1 = e1 ./ vecnorm(e1, 2, 2);  G.e2 = cross(G.u, G.e1, 2);
    tri = zeros(0, 3);
    for a = 1:size(A, 1)
        nb = find(A(a, :));  nb = nb(nb > a);
        for p = nb, q = nb(nb > p);  q = q(A(p, q));  tri = [tri; repmat([a p], numel(q), 1) q(:)]; end %#ok<AGROW>
    end
    G.tri = tri;
end

function P = i_phys(m, G, fc, depth, rdepth)
    amp = abs(m) / max(abs(m));  ph = angle(m);  on = amp >= 0.5;  i = G.e(:, 1);  j = G.e(:, 2);
    P.extent_mm2 = sum(G.area(on));  P.extent_radius_mm = sqrt(P.extent_mm2 / pi);
    mn = m / max(abs(m));
    lam = 2 * mean(abs(mn(i) - mn(j)).^2 ./ G.d.^2) / mean((abs(mn(i)).^2 + abs(mn(j)).^2) / 2);
    P.dirichlet_wavelength_mm = 2*pi / sqrt(lam);  P.tile_spacing_mm = median(G.d);
    P.pg = i_pg(m, G);
    w = @(x) angle(exp(1i * x));
    g3 = nan(numel(m), 3);  gx = nan(numel(m), 3);  az = nan(numel(m), 1);
    for t = find(on)'
        nb = find(G.A(t, :) & amp' >= 0.25);
        if numel(nb) < 2, continue; end
        dp = w(ph(nb) - ph(t));
        g3(t, :) = (pinv(G.t(nb, :) - G.t(t, :)) * dp)';            % rad / mm, template coordinates
        gx(t, :) = (pinv(G.x(nb, :) - G.x(t, :)) * dp)';            % rad / mm, the participant's cortex
        D = G.u(nb, :) - G.u(t, :);  gs = [D * G.e1(t, :)' D * G.e2(t, :)'] \ dp;
        az(t) = atan2d(-gs(2), -gs(1));
    end
    ok = all(isfinite(g3), 2) & vecnorm(g3, 2, 2) > 0;
    v = -g3(ok, :) ./ vecnorm(g3(ok, :), 2, 2);  wa = amp(ok);
    mv = sum(v .* wa, 1) / max(sum(wa), eps);  if ~any(ok), mv = nan(1, 3); end
    P.dir_x = mv(1);  P.dir_y = mv(2);  P.dir_z = mv(3);  P.dir_concentration = norm(mv);
    P.dir_azimuth_deg = angle(sum(wa .* exp(1i * deg2rad(az(ok))))) * 180/pi;
    gm = median(vecnorm(gx(ok, :), 2, 2));
    P.phase_wavelength_mm = 2*pi / gm;  P.speed_m_s = fc * P.phase_wavelength_mm / 1e3;
    P.speed_resolved = double(depth >= rdepth);
    Lp = zeros(numel(m), 1);
    for e = 1:size(G.e, 1)
        a = i(e);  b = j(e);  dd = w(ph(b) - ph(a)) / G.d(e)^2;  Lp(a) = Lp(a) + dd;  Lp(b) = Lp(b) - dd;
    end
    P.div_index = -sum(amp(on) .* Lp(on)) / max(sum(amp(on)), eps);
    wn = 0;  wc = 0;
    for r = 1:size(G.tri, 1)
        t = G.tri(r, :);  if any(amp(t) < 0.25), continue; end
        k = round((w(ph(t(2)) - ph(t(1))) + w(ph(t(3)) - ph(t(2))) + w(ph(t(1)) - ph(t(3)))) / (2*pi));
        k = k * sign(det(G.u(t, :)));  wn = wn + k;  wc = wc + abs(k);
    end
    P.winding_net = wn;  P.winding_count = wc;
end

function s = i_pg(m, G)
% amplitude-weighted mean |phase difference| over edges between tiles at >= half max; NaN with < 3 such tiles
    amp = abs(m) / max(abs(m));  on = amp >= 0.5;  s = NaN;
    if nnz(on) < 3, return; end
    i = G.e(:, 1);  j = G.e(:, 2);  k = on(i) & on(j);
    if ~any(k), return; end
    wt = amp(i(k)) .* amp(j(k));  s = sum(wt .* abs(angle(m(i(k)) .* conj(m(j(k)))))) / sum(wt);
end

function [ro, rm, p] = i_spin(m, c, G, n, rs)
% Alexander-Bloch spin on the template sphere, the right hemisphere mirrored (x -> -x); within each hemisphere
% (x depth, for joint atoms) group, one rotation per spin for all groups
    nr = @(a, b) abs(a' * b) / (norm(a) * norm(b));  ro = nr(c, m);  r = zeros(n, 1);  Fx = diag([-1 1 1]);
    for s = 1:n
        [Q, Rq] = qr(randn(rs, 3));  Q = Q * diag(sign(diag(Rq)));  if det(Q) < 0, Q(:, 1) = -Q(:, 1); end
        ms = m;
        for g = 1:numel(G.groups)
            ix = G.groups{g};  Qg = Q;  if G.hemi(ix(1)) == "R", Qg = Fx * Q * Fx; end
            ms(ix) = m(ix(knnsearch(G.u(ix, :), G.u(ix, :) * Qg)));
        end
        r(s) = nr(c, ms);
    end
    rm = mean(r);  p = (1 + nnz(r >= ro)) / (1 + n);
end

function L = i_lagged(m, R, G, zk)
% N3: the best two-phase source through the instrument (rheome.modes.lagged), and the same model SIMULATED in
% the participant's own data -- the mode's component of every sample in its state replaced by the fit, the rest
% (noise, other structure) kept -- re-estimated and refitted, so the estimator's own misfit is the null's.
% A state the back-fit gave < 2 samples has no mode (all zeros): no N3, every field NaN (it read r_n3 0 and
% r_n3_sim 1 with 0 samples, and crashed rheome.modes.lagged on a NaN mode with 1 sample)
    if ~any(m)
        L = struct('r_n3', NaN, 'r_n3_sim', NaN, 'pg_n3', NaN, 'n3_lag_rad', NaN, 'n3_tiles_a', NaN, 'n3_tiles_b', NaN);  return
    end
    F = rheome.modes.lagged(m, R);  f = F.fit / norm(F.fit);  f = f * exp(-1i * angle(m' * f));  mu = m / norm(m);
    a = mu' * zk;  zs = zk - mu * a + f * a;
    [V, Lm] = eig(zs * zs', 'vector');  [~, j] = max(real(Lm));  Fs = rheome.modes.lagged(V(:, j), R);
    L.r_n3 = F.r_n3;  L.r_n3_sim = Fs.r_n3;  L.pg_n3 = i_pg(F.fit, G);
    L.n3_lag_rad = F.lag_rad;  L.n3_tiles_a = nnz(F.x1 > 0);  L.n3_tiles_b = nnz(F.x2 > 0);
end

function [names, P] = i_planted(S, W)
% the N3 calibration patterns (protocol section 13) on this head, as in the methods benchmark:
% a 200 mm plane wave travelling A->P along SCS x over the whole cortex, and two lagged pairs of 30 mm blobs
% (frontal-occipital, frontal leading by 90 deg; lateral, 60 deg) -- the normal current seen through W * gain
    V = S.Sf.Vertices;  nr = S.Sf.VertNormals ./ vecnorm(S.Sf.VertNormals, 2, 2);  G = S.G;
    Gn = nr(:,1)' .* G(:, 1:3:end) + nr(:,2)' .* G(:, 2:3:end) + nr(:,3)' .* G(:, 3:3:end);
    c0 = mean(V(V(:,1) <= prctile(V(:,1), 10), :), 1);  c1 = mean(V(V(:,1) >= prctile(V(:,1), 90), :), 1);
    g = @(c) exp(-vecnorm(V - c, 2, 2).^2 / (2 * 0.030^2));
    P = W * Gn * [exp(1i * 2*pi/0.200 * V(:,1)), g(c1) + exp(-1i*pi/2) * g(c0), ...
                  g(c1 + [0 0.04 0]) + exp(-1i*pi/3) * g(c0 - [0 0.04 0])];
    names = ["planted_wave200_ap" "planted_pair_fo" "planted_pair_lateral"];
end

function L = i_n3plant(p, R, G, z)
% the participant's leading component in every valid sample replaced by the planted pattern (its coefficients,
% the noise and the rest kept), the mode re-estimated and read by the same N3 as a real mode: n3_excess =
% r_n3_sim - r_n3, the quantity the travelling label tests against N3Margin
    [V, Lm] = eig(z * z', 'vector');  [~, j] = max(real(Lm));  mu = V(:, j);
    p = p / norm(p);  a = mu' * z;  zp = z - mu * a + p * a;
    [V, Lm] = eig(zp * zp', 'vector');  [~, j] = max(real(Lm));  m = V(:, j);
    L = i_lagged(m, R, G, zp);  L.r_plant = abs(p' * m);  L.n3_excess = L.r_n3_sim - L.r_n3;
end

function u = i_unit(f)
    persistent U
    if isempty(U)
        U = dictionary(["extent_mm2" "extent_radius_mm" "dirichlet_wavelength_mm" "tile_spacing_mm" "pg" ...
            "dir_x" "dir_y" "dir_z" "dir_concentration" "dir_azimuth_deg" "phase_wavelength_mm" "speed_m_s" ...
            "speed_resolved" "div_index" "winding_net" "winding_count" "pg_n3" "r_n3" "r_n3_sim" "n3_lag_rad" ...
            "n3_tiles_a" "n3_tiles_b" "r_plant" "n3_excess" "dominant_depth" "scale_mm" "t_bar" "scale_edge_mm"], ["mm2" "mm" "mm" "mm" "rad" "unit" "unit" "unit" "r" "deg" "mm" "m/s" ...
            "flag" "rad/mm2" "windings" "windings" "rad" "r" "r" "rad" "tiles" "tiles" "r" "r" "depth" "mm" "T" "mm"]);
    end
    u = U(f);
end

% Author: Diellor Basha, 2026
