function [T, X] = measure_atlas(name, S, item, opts)
% SCALE.MEASURE_ATLAS  MS1 group P8 (G12) on one participant: correspondence, atlas events, the shared-gauge
% alpha tensor and the connectome roll-up -- the sub-0002 atlas, gauge and multimodal checks, per participant.
%
%   [T, X] = rheome.scale.measure_atlas(name, S, "correspondence")
%   [T, X] = rheome.scale.measure_atlas(name, S, "atlasevents", NSurr=20)
%   [T, X] = rheome.scale.measure_atlas(name, S, "gaugetensor")
%   [T, X] = rheome.scale.measure_atlas(name, S, "connectome", Gammas=[0 0.01 0.1 1 10])
%
% Spec: desk/evidence/3377657c.../group-analysis-plan.md section 2 G12 (request d447727e). The geometry parts
% of G12 (roll-up exactness, atom-tile overlap, gauge singularity count and smoothness) are P2's
% rheome.scale.measure_geometry; this file is the rest. Every item reads the paper's kernel, S.Res (the
% unconstrained, noise-whitened minimum norm of rheome.scale.sensors, SnrFixed 3), never dSPM.
%
% ⭐ THE SHARED SPACE IS A TEMPLATE CORTEX REACHED THROUGH THE REGISTERED SPHERE. A Brainstorm template
% folder (default: RHEOME_TEMPLATE, else RHEOME_BRAINSTORM/defaults/anatomy/ICBM152) gives
% tess_cortex_pial_low.mat (Reg.Sphere, Desikan-Killiany) and subjectimage_T1.mat (SCS, NCS). The group
% tiles are rheome.geom.tree on the template's hemispheres (depth <= max(Depths)); a participant's vertex
% belongs to the tile of its nearest template vertex on the sphere, so tile k is the same anatomical
% place in every participant (the atlas's group tree, VR1 2.9, with Rheome's tree instead of nsp's ladder).
% The group gauge is rheome.geom.sphereframe: the sphere's meridian pushed to the cortex (21cc8d9b).
%
% ITEMS (rows of T; X is the item's per-unit table, written as <item>.csv by rheome.scale.run):
%   correspondence  per hemisphere: Desikan-Killiany agreement between the participant's own labels and the
%                   template's at the nearest vertex THROUGH THE SPHERE (dk_agree_sphere, by vertex and by
%                   area) and through RAW COORDINATES (dk_agree_raw: nearest template vertex in MNI mm, each
%                   cortex placed by its own MRI's SCS -> NCS, Brainstorm's linear MNI transform). Frame
%                   deviation from the shared meridian: the trivial connection solved on the cortex with +1 at
%                   the faces of the sphere's two poles (rheome.operators.gauge) against sphereframe's e1, after
%                   the one global rotation the connection leaves free: median, p90, p99 in degrees.
%                   (VR1 2.9 on 5 PREVENT-AD stores: 94 % through the sphere; 1.3-1.7 deg median.)
%   atlasevents     Fig. 10 (F12) per hemisphere x depth in Depths: alpha (Band, 6.1-12.1 Hz as F12) power of
%                   the tile-summed current per joint cell (group tile x TimeTileS = 0.1648 s, tower level
%                   -19), in dB against the tile's own median; an EVENT is a connected set of cells above
%                   ThresholdDB (5 dB), linked through touching tiles (a shared mesh edge) in one time tile
%                   and through the same or a touching tile in consecutive time tiles. Rows: events_per_min,
%                   duration_median_s, tiles_visited_median, frac_on, the same averaged over NSurr
%                   phase-randomised surrogates (one random phase per frequency bin shared by every tile
%                   series, so every cross-spectrum is kept), and p_rate = (1 + #{surrogate rate >= data})
%                   / (1 + NSurr).
%   gaugetensor     section 12.3: the alpha (TensorBand) orientation tensor of the tangential current in the
%                   shared gauge, per depth-TensorDepth group tile: area-weighted sum over the tile's vertices
%                   of E' K C K' E, C the sensors' band cross-spectrum, E = [north west]. Per tile: the
%                   trace-normalised components t11 (north-north), t22, t12, anisotropy (l1-l2)/(l1+l2), axis
%                   (deg, 0 = north, 90 = west) and the normal share of the band power.
%   connectome      ea515f28 per participant. (a) Destrieux -> DK against the exact dyadic roll-up (depth 6 ->
%                   5, the participant's own rheome.geom.tree), as F11 C: area outside the parent, share of
%                   finer units > 5 % outside, edge-weight loss ||Ag' C_fine Ag - C_coarse||_F / ||C_coarse||_F
%                   and the CV of the finer units' areas, on the alpha co-power connectome of the kernel and,
%                   with fibres, on the fibre-count connectome. (b) With fibres (rheome.connectome.resolve: the
%                   participant's OWN tractography first, else the HCP-1065 template): the connectome wavelet
%                   (mexhat on the LB-Connectome pencil) per member at Gammas x the auto gamma, over NSeeds
%                   seeds: area-weighted RMS width (mm, Euclidean from the seed), energy in the other
%                   hemisphere, energy beyond 3 sigma (VR1 2.1d), and the coupling scale 2 pi / sqrt(lambda_2).
%
% ⚠ RAW COORDINATES NEED THE MRI. The subject's SCS/NCS come from the protocol's subjectimage (staged by nsp
% cf-atlas, cached as mri.mat by rheome.scale.importsubject); without it dk_agree_raw is NaN. NCS is
% Brainstorm's linear MNI fit; ICBM152 is MNI by construction.
% ⚠ THE POLAR TILE MIXES DIRECTIONS. The +z pole is precentral; a depth-3 tile around it sees the meridian
% turn through 360 deg. colat_median_deg is in the gaugetensor table so the group step can flag it.
% ⚠ Group statistics (paired Wilcoxon of data vs surrogate rates, sphere vs raw; bootstrap CI per tile) are
% the reduction's, across participants; this returns one participant's numbers.
%
% See also: rheome.scale.measure_geometry, rheome.geom.sphereframe, rheome.geom.tree, rheome.operators.gauge,
%           rheome.connectome.resolve, rheome.operators.lb_connectome, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        item (1,1) string {mustBeMember(item, ["correspondence" "atlasevents" "gaugetensor" "connectome"])} = "correspondence"
        opts.Template char = ''
        opts.Depths double = 1:3
        opts.Band (1,2) double = [6.1 12.1]
        opts.TimeTileS (1,1) double = 0.1648
        opts.ThresholdDB (1,1) double = 5
        opts.NSurr (1,1) double = 20
        opts.Seed (1,1) double = 1
        opts.TensorBand (1,2) double = [8 12]
        opts.TensorDepth (1,1) double = 3
        opts.Gammas double = [0 0.01 0.1 1 10]
        opts.Kc (1,1) double = 400
        opts.NSeeds (1,1) double = 60
        opts.Hemis string = ["L" "R"]
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    dmax = max([opts.Depths opts.TensorDepth]);
    switch item
        case "correspondence", Tp = i_template(opts, dmax);  [T, X] = i_correspondence(name, S, Tp, opts);
        case "atlasevents",    Tp = i_template(opts, dmax);  [T, X] = i_events(name, S, Tp, opts);
        case "gaugetensor",    Tp = i_template(opts, dmax);  [T, X] = i_tensor(name, S, Tp, opts);
        case "connectome",     [T, X] = i_connectome(name, S, opts);
    end
end

%% ---------- correspondence: DK through the sphere vs raw coordinates; frame vs the shared meridian
function [T, X] = i_correspondence(name, S, Tp, o)
    labS = i_labels(rheome.load.atlas(name, 'Desikan-Killiany'), S.Sf.nV);
    f = fullfile(rheome.load.root(), name, 'mri.mat');  mri = [];  if isfile(f), mri = load(f); end
    T = rheome.scale.rows("correspondence", strings(0,1), [], "");  X = table();
    for h = o.Hemis(:)'
        Hh = S.B.(h);  gv = Hh.gv(:);  th = Tp.(h);  a = full(sum(Hh.lbo.Mass, 2));
        ls = labS(gv);
        [is, ds] = knnsearch(th.u, i_unit(S.Sf.Sphere(gv, :)));
        [agS, agSa] = i_agree(ls, Tp.dk(th.gv(is)), a);
        agR = NaN;  agRa = NaN;  dR = NaN;
        ps = i_mni(S.Sf.Vertices(gv, :), mri);  pt = i_mni(Tp.S.Vertices(th.gv, :), Tp.mri);
        if ~isempty(ps) && ~isempty(pt)
            [ir, dr] = knnsearch(pt, ps);  [agR, agRa] = i_agree(ls, Tp.dk(th.gv(ir)), a);  dR = median(dr);
        end

        % the frame: trivial connection with +1 at the sphere's poles vs the pushed-forward meridian
        V = double(Hh.S.Vertices);  Fc = double(Hh.S.Faces);  Sph = S.Sf.Sphere(gv, :);
        fr0 = rheome.geom.sphereframe(V, Fc, i_normals(S.Sf, gv, V, Fc), Sph);
        fp = [find(any(Fc == fr0.north, 2), 1) find(any(Fc == fr0.south, 2), 1)];
        g = rheome.operators.gauge(V, Fc, Method="trivial", Singular=fp);
        fr = rheome.geom.sphereframe(V, Fc, g.normal, Sph);
        ring = false(size(V, 1), 1);  ring(Fc(g.singular, :)) = true;
        ok = ~fr.singular & ~ring & all(isfinite(g.e1), 2);
        th1 = atan2(sum(g.e1(ok,:) .* fr.e2(ok,:), 2), sum(g.e1(ok,:) .* fr.e1(ok,:), 2));
        dev = abs(angle(exp(1i * (th1 - angle(mean(exp(1i * th1))))))) * 180/pi;

        T = [T; rheome.scale.rows("correspondence", ...
             ["dk_agree_sphere" "dk_agree_sphere_area" "dk_agree_raw" "dk_agree_raw_area" "sphere_nn_rad_median" ...
              "raw_nn_mm_median" "frame_dev_median_deg" "frame_dev_p90_deg" "frame_dev_p99_deg" "frame_singular_vertices" "gauge_singular_faces"], ...
             [agS agSa agR agRa median(ds) dR median(dev) prctile(dev, 90) prctile(dev, 99) nnz(fr.singular) numel(g.singular)], ...
             ["fraction" "fraction" "fraction" "fraction" "rad" "mm" "deg" "deg" "deg" "vertices" "faces"], h)]; %#ok<AGROW>
        X = [X; table(h, agS, agSa, agR, agRa, median(ds), dR, median(dev), prctile(dev, 90), prctile(dev, 99), nnz(fr.singular), numel(g.singular), ...
             'VariableNames', {'hemi','dk_agree_sphere','dk_agree_sphere_area','dk_agree_raw','dk_agree_raw_area','sphere_nn_rad_median', ...
                               'raw_nn_mm_median','frame_dev_median_deg','frame_dev_p90_deg','frame_dev_p99_deg','frame_singular_vertices','gauge_singular_faces'})]; %#ok<AGROW>
        fprintf('[atlas %s] %s: DK sphere %.3f raw %.3f; frame %.2f deg median\n', name, h, agS, agR, median(dev));
    end
end

function [f, fa] = i_agree(ls, lt, a)
    ok = ls ~= "" & lt ~= "";  m = ok & ls == lt;
    f = nnz(m) / nnz(ok);  fa = sum(a(m)) / sum(a(ok));
end

%% ---------- atlas events (Fig. 10) on the group tiles, against phase-randomised surrogates
function [T, X] = i_events(name, S, Tp, o)
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    N = size(F, 2);  K = S.Res.ImagingKernel;  Y = [];  blk = struct('h', {}, 'l', {}, 'rows', {}, 'A', {});
    for h = o.Hemis(:)'
        gv = S.B.(h).gv(:);
        for l = o.Depths
            [P, A] = i_tiles(S, Tp, h, l);  r0 = size(Y, 1);
            for k = 1:3, Y = [Y; (P' * K((gv - 1)*3 + k, :)) * F]; end %#ok<AGROW>
            blk(end+1) = struct('h', h, 'l', l, 'rows', r0 + (1:3*size(P, 2)), 'A', A); %#ok<AGROW>
        end
    end
    clear F
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', o.Band(1), ...
                    'HalfPowerFrequency2', o.Band(2), 'SampleRate', fs);
    nb = round(o.TimeTileS * fs);  mins = N / fs / 60;
    D = i_detect(Y, blk, bp, nb, o.ThresholdDB, mins);
    Su = cell(1, o.NSurr);
    for s = 1:o.NSurr, rng(o.Seed + s);  Su{s} = i_detect(i_phaserand(Y), blk, bp, nb, o.ThresholdDB, mins); end
    T = rheome.scale.rows("atlasevents", strings(0,1), [], "");  X = table();
    for b = 1:numel(blk)
        d = D(b, :);  sr = cell2mat(cellfun(@(z) z(b, :), Su, 'UniformOutput', false)');
        lab = blk(b).h + ".d" + blk(b).l;
        T = [T; rheome.scale.rows("atlasevents", ...
             ["events_per_min" "duration_median_s" "tiles_visited_median" "frac_on" ...
              "surr_events_per_min" "surr_duration_median_s" "surr_tiles_visited_median" "surr_frac_on" "p_rate"], ...
             [d(2:5) mean(sr(:, 2:5), 1, 'omitnan') (1 + nnz(sr(:, 2) >= d(2))) / (1 + o.NSurr)], ...
             ["1/min" "s" "tiles" "fraction" "1/min" "s" "tiles" "fraction" "p"], lab)]; %#ok<AGROW>
        src = ["data"; repmat("surrogate", o.NSurr, 1)];  M = [d; sr];
        X = [X; table(repmat(blk(b).h, o.NSurr + 1, 1), repmat(blk(b).l, o.NSurr + 1, 1), src, (0:o.NSurr)', ...
             M(:,1), M(:,2), M(:,3), M(:,4), M(:,5), M(:,6), 'VariableNames', ...
             {'hemi','depth','source','surrogate','n_events','events_per_min','duration_median_s','tiles_visited_median','frac_on','duration_p90_s'})]; %#ok<AGROW>
        fprintf('[atlas %s] events %s: %.1f /min (surrogates %.1f), %.2f s, %.1f tiles\n', name, lab, d(2), mean(sr(:, 2)), d(3), d(4));
    end
end

function D = i_detect(Y, blk, bp, nb, thr, mins)
% per block: [n_events events_per_min duration_median_s tiles_visited_median frac_on duration_p90_s]
    E = abs(hilbert(filtfilt(bp, Y.')).').^2;  nbin = floor(size(E, 2) / nb);  dt = nb / bp.SampleRate;
    D = zeros(numel(blk), 6);
    for b = 1:numel(blk)
        nT = numel(blk(b).rows) / 3;  pw = 0;
        for k = 1:3, pw = pw + E(blk(b).rows((k-1)*nT + (1:nT)), :); end
        pw = squeeze(mean(reshape(pw(:, 1:nbin*nb), nT, nb, nbin), 2));  pw = reshape(pw, nT, nbin);
        on = 10*log10(pw ./ median(pw, 2)) > thr;
        [ev, dur, nt] = i_events_of(on, blk(b).A, dt);
        D(b, :) = [ev, ev / mins, median(dur), median(nt), mean(on, 'all'), prctile(dur, 90)];
        if ev == 0, D(b, [3 4 6]) = NaN; end
    end
end

function [n, dur, ntiles] = i_events_of(on, A, dt)
% connected sets of on-cells: touching tiles within a time tile, same or touching tile across consecutive ones
    [nT, nb] = size(on);  id = reshape(1:nT*nb, nT, nb);
    % ⚠ e(:), t(:): with one edge (depth 1, two tiles) on(i1,:) is a row and find returns rows (MATLAB:catenate)
    [i1, j1] = find(triu(A));  [e, t] = find(on(i1, :) & on(j1, :));  e = e(:);  t = t(:);
    s = id(sub2ind([nT nb], i1(e), t));  d = id(sub2ind([nT nb], j1(e), t));
    [i2, j2] = find(A | speye(nT));  [e, t] = find(on(i2, 1:end-1) & on(j2, 2:end));  e = e(:);  t = t(:);
    s = [s; id(sub2ind([nT nb], i2(e), t))];  d = [d; id(sub2ind([nT nb], j2(e), t + 1))];
    c = conncomp(graph(sparse(s, d, 1, nT*nb, nT*nb) + sparse(d, s, 1, nT*nb, nT*nb)));
    v = find(on(:));  [~, ~, g] = unique(c(v));  n = max([g; 0]);
    if n == 0, dur = NaN;  ntiles = NaN;  return; end
    tile = mod(v - 1, nT) + 1;  bin = floor((v - 1) / nT) + 1;
    dur = (accumarray(g, bin, [], @max) - accumarray(g, bin, [], @min) + 1) * dt;
    ntiles = accumarray(g, tile, [], @(x) numel(unique(x)));
end

function Y = i_phaserand(Y)
% one random phase per frequency bin, shared by every row: every auto- and cross-spectrum is kept
    N = size(Y, 2);  Z = fft(Y, [], 2);  h = floor((N - 1) / 2);
    Z(:, 2:h+1) = Z(:, 2:h+1) .* exp(2i*pi*rand(1, h));  Z(:, N-h+1:N) = conj(fliplr(Z(:, 2:h+1)));
    Y = real(ifft(Z, [], 2));
end

%% ---------- the alpha orientation tensor per group tile, in the shared gauge
function [T, X] = i_tensor(name, S, Tp, o)
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', o.TensorBand(1), ...
                    'HalfPowerFrequency2', o.TensorBand(2), 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';  C = real(Fa * Fa') / size(Fa, 2);  clear F Fa
    K = S.Res.ImagingKernel;  T = rheome.scale.rows("gaugetensor", strings(0,1), [], "");  X = table();
    for h = o.Hemis(:)'
        Hh = S.B.(h);  gv = Hh.gv(:);  V = double(Hh.S.Vertices);  Fc = double(Hh.S.Faces);
        Nv = i_normals(S.Sf, gv, V, Fc);  fr = rheome.geom.sphereframe(V, Fc, Nv, S.Sf.Sphere(gv, :));
        ok = ~fr.singular;  a = full(sum(Hh.lbo.Mass, 2)) .* ok;
        Kx = K((gv-1)*3 + 1, :);  Ky = K((gv-1)*3 + 2, :);  Kz = K((gv-1)*3 + 3, :);
        pr = @(e) e(:,1) .* Kx + e(:,2) .* Ky + e(:,3) .* Kz;
        e1 = fr.e1;  e2 = fr.e2;  e1(~ok, :) = 0;  e2(~ok, :) = 0;
        K1 = pr(e1);  K2 = pr(e2);  Kn = pr(Nv ./ vecnorm(Nv, 2, 2));
        K1C = K1 * C;  q = @(x) x .* a;
        t11 = q(sum(K1C .* K1, 2));  t12 = q(sum(K1C .* K2, 2));  t22 = q(sum((K2 * C) .* K2, 2));  tn = q(sum((Kn * C) .* Kn, 2));
        [P, ~, ids] = i_tiles(S, Tp, h, o.TensorDepth);
        s11 = P' * t11;  s22 = P' * t22;  s12 = P' * t12;  sn = P' * tn;  tr = s11 + s22;
        an = sqrt((s11 - s22).^2 + 4*s12.^2) ./ tr;  ax = 0.5 * atan2d(2*s12, s11 - s22);
        cl = arrayfun(@(k) median(fr.colat(P(:, k) > 0)) * 180/pi, (1:size(P, 2))');
        X = [X; table(repmat(h, numel(ids), 1), repmat(o.TensorDepth, numel(ids), 1), ids(:), full(sum(P, 1))', ...
             s11 ./ tr, s22 ./ tr, s12 ./ tr, an, ax, sn ./ (sn + tr), tr, cl, 'VariableNames', ...
             {'hemi','depth','node_id','n_vertices','t11','t22','t12','anisotropy','axis_deg','normal_share','trace','colat_median_deg'})]; %#ok<AGROW>
        T = [T; rheome.scale.rows("gaugetensor", ["anisotropy_median" "normal_share_median" "n_tiles"], ...
             [median(an, 'omitnan') median(sn ./ (sn + tr), 'omitnan') numel(ids)], ["ratio" "fraction" "tiles"], h)]; %#ok<AGROW>
        fprintf('[atlas %s] tensor %s: %d tiles, anisotropy %.2f, normal share %.2f\n', name, h, numel(ids), median(an), median(sn ./ (sn + tr)));
    end
end

%% ---------- connectome: Destrieux -> DK vs the dyadic roll-up; connectome-wavelet widths with a gamma sweep
function [T, X] = i_connectome(name, S, o)
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));  clear st
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', o.Band(1), ...
                    'HalfPowerFrequency2', o.Band(2), 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';  Cs = real(Fa * Fa') / size(Fa, 2);  clear F Fa
    Sf = S.Sf;  nV = Sf.nV;  K = S.Res.ImagingKernel;
    [L, M] = rheome.operators.laplace_beltrami(Sf.Vertices, Sf.Faces);  a = full(sum(M, 2));  clear L M
    Pdk = i_member(rheome.load.atlas(name, 'Desikan-Killiany'), nV);
    Pdx = i_member(rheome.load.atlas(name, 'Destrieux'), nV);
    P6 = [];  P5 = [];
    for h = o.Hemis(:)'                                     % the participant's own dyadic tree, both hemispheres
        Hh = S.B.(h);  gv = Hh.gv(:);
        Tr = rheome.geom.tree(Hh.S, L=Hh.lbo.L, M=Hh.lbo.M, MaxDepth=6);
        G6 = rheome.geom.tiles(Tr, Hh.S, 6);  G5 = rheome.geom.tiles(Tr, Hh.S, 5);
        if ~G6.heap, error('scale:atlas:heap', 'hemisphere %s: tree not full to depth 6', h); end
        P6 = [P6, sparse(gv, G6.tileOf, 1, nV, numel(G6.node_id))]; %#ok<AGROW>
        P5 = [P5, sparse(gv, G5.tileOf, 1, nV, numel(G5.node_id))]; %#ok<AGROW>
    end
    cp = @(P) i_copower(P, K, Cs);
    T = [i_zoom("destrieux_dk", Pdx, Pdk, a, cp(Pdx), cp(Pdk), "copower"); ...
         i_zoom("dyadic_6_5", P6, P5, a, cp(P6), cp(P5), "copower")];
    X = table();
    try
        [ep, prov] = rheome.connectome.resolve(name, Sf);
    catch e
        T = [T; rheome.scale.rows("connectome", "fibres_available", 0, "flag", string(e.identifier))];
        fprintf('[atlas %s] connectome: no fibres (%s)\n', name, e.message);  return
    end
    [W, keep] = rheome.operators.connectome(Sf.Vertices, Sf.Faces, ep);
    fw = @(P) full(P' * W * P);
    hemi = zeros(nV, 1);  hemi(Sf.Hemi{1}) = 1;  hemi(Sf.Hemi{2}) = 2;
    [ii, jj, ww] = find(triu(W));
    T = [T; rheome.scale.rows("connectome", ["fibres_available" "own_tractography" "n_fibres" "interhemispheric_weight"], ...
         [1 string(prov.fiberSource) == "subject-specific" prov.numFibers sum(ww(hemi(ii) ~= hemi(jj))) / sum(ww)], ...
         ["flag" "flag" "fibres" "fraction"], string(prov.fiberSource)); ...
         i_zoom("destrieux_dk", Pdx, Pdk, a, fw(Pdx), fw(Pdk), "fibres"); i_zoom("dyadic_6_5", P6, P5, a, fw(P6), fw(P5), "fibres")];

    rng(o.Seed);  seeds = keep(randperm(numel(keep), min(o.NSeeds, numel(keep))));
    [~, ~, g0] = rheome.operators.lb_connectome(Sf, W, keep);
    for gm = o.Gammas(:)'
        [A, B, g] = rheome.operators.lb_connectome(Sf, W, keep, gm * g0);
        b = rheome.eigen.modes(A, B, o.Kc);  lam = b.Lambda(:);  ab = full(sum(B, 2));
        gfb = rheome.graphfilterbank(lam, 'Wavelet', 'mexhat', 'VoicesPerOctave', 2, 'SizeLimits', [12 96]*1e-3, ...
                                     'Transform', rheome.graphtransform.eigen(b.Phi, B, lam));
        w = widths(gfb);  At = graphfilters(gfb, 'Type', 'vertex', 'Vertex', seeds);
        lam2 = lam(find(lam > 1e-6 * max(lam), 1));
        for m = 1:size(At, 3)
            E = At(:, :, m).^2 .* ab;  Et = sum(E, 1);  rw = zeros(1, numel(seeds));  ct = rw;  rm = rw;
            for p = 1:numel(seeds)
                dd = vecnorm(Sf.Vertices - Sf.Vertices(seeds(p), :), 2, 2);
                rw(p) = sqrt(sum(E(:, p) .* dd.^2) / Et(p));  ct(p) = sum(E(hemi ~= hemi(seeds(p)), p)) / Et(p);
                rm(p) = sum(E(dd > 3*w(m), p)) / Et(p);
            end
            X = [X; table(gm, g, m, 1e3*w(m), 1e3*median(rw), median(ct), median(rm), 2e3*pi/sqrt(lam2), 'VariableNames', ...
                 {'gamma_mult','gamma','member','sigma_mm','rms_width_mm','contra_frac','beyond_3sigma_frac','coupling_mm'})]; %#ok<AGROW>
        end
        T = [T; rheome.scale.rows("connectome", "coupling_mm", 2e3*pi/sqrt(lam2), "mm", "gamma_x" + gm)]; %#ok<AGROW>
        fprintf('[atlas %s] connectome gamma x%g: coupling %.0f mm\n', name, gm, 2e3*pi/sqrt(lam2));
    end
end

function C = i_copower(P, K, Cs)
% alpha co-power between units: sum over x, y, z of (P' K_k) Cs (P' K_k)' -- additive over vertices, so exact roll-up
    C = 0;  nV = size(P, 1);
    for k = 1:3, Kp = P' * K((1:nV)*3 - 3 + k, :);  C = C + Kp * Cs * Kp'; end
end

function T = i_zoom(lab, Pf, Pc, a, Cf, Cc, src)
% F11 C: each finer unit goes to the coarse unit holding most of its area
    Ov = full(Pf' * (a .* Pc));  [mx, par] = max(Ov, [], 2);  af = full(Pf' * a);
    Ag = sparse(1:numel(par), par, 1, numel(par), size(Pc, 2));
    loss = norm(Ag' * Cf * Ag - Cc, 'fro') / norm(Cc, 'fro');
    T = rheome.scale.rows("connectome", ["area_outside" "split_share" "edge_loss" "cv_area" "n_fine" "n_coarse"], ...
        [sum(af - mx) / sum(af), mean(1 - mx ./ af > 0.05), loss, std(af) / mean(af), size(Pf, 2), size(Pc, 2)], ...
        ["fraction" "fraction" "relative" "ratio" "units" "units"], lab + "_" + src);
end

%% ---------- shared helpers
function Tp = i_template(o, dmax)
    persistent cache key
    d = o.Template;
    if isempty(d), d = getenv('RHEOME_TEMPLATE'); end
    if isempty(d) && ~isempty(getenv('RHEOME_BRAINSTORM')), d = fullfile(getenv('RHEOME_BRAINSTORM'), 'defaults', 'anatomy', 'ICBM152'); end
    k = string(d) + "|" + dmax + "|" + strjoin(o.Hemis, ",");
    if isequal(key, k), Tp = cache;  return; end
    f = fullfile(d, 'tess_cortex_pial_low.mat');
    if ~isfile(f)
        error('scale:atlas:template', ['No template cortex at "%s": pass Template=, or set RHEOME_TEMPLATE, or ' ...
              'RHEOME_BRAINSTORM (a Brainstorm checkout with defaults/anatomy/ICBM152).'], f);
    end
    Tp.S = rheome.io.read.surface(f);
    if isempty(Tp.S.Sphere) || isempty(Tp.S.Hemi), error('scale:atlas:template', '%s has no Reg.Sphere or no L/R split', f); end
    A = rheome.io.read.atlas(f);  Tp.dk = i_labels(A(strcmpi({A.Name}, 'Desikan-Killiany')), Tp.S.nV);
    Tp.mri = [];  mf = fullfile(d, 'subjectimage_T1.mat');  if isfile(mf), Tp.mri = load(mf, 'SCS', 'NCS'); end
    Tp.file = f;
    for h = o.Hemis(:)'
        Sh = rheome.utils.hemisphere(Tp.S, [], char(h));
        Tr = rheome.geom.tree(Sh, MaxDepth=dmax);  G = rheome.geom.tiles(Tr, Sh, dmax);
        if ~G.heap || any(G.depth ~= dmax), error('scale:atlas:heap', 'template hemisphere %s: tree not full to depth %d', h, dmax); end
        Tp.(h) = struct('gv', Sh.GlobalVertices(:), 'u', i_unit(Tp.S.Sphere(Sh.GlobalVertices, :)), ...
                        'code', G.node_id(G.tileOf), 'dmax', dmax);
    end
    cache = Tp;  key = k;
end

function [P, A, ids] = i_tiles(S, Tp, h, l)
% the participant's vertices in the group tiles of depth l (nearest template vertex on the sphere)
    th = Tp.(h);  gv = S.B.(h).gv(:);  Fc = double(S.B.(h).S.Faces);
    c = th.code(knnsearch(th.u, i_unit(S.Sf.Sphere(gv, :))));  sh = 2^(th.dmax - l);
    ids = unique(floor(th.code / sh));  [~, ti] = ismember(floor(c / sh), ids);
    P = sparse(1:numel(gv), ti, 1, numel(gv), numel(ids));
    E = [Fc(:,[1 2]); Fc(:,[2 3]); Fc(:,[3 1])];  x = ti(E(:,1)) ~= ti(E(:,2));
    A = sparse(ti(E(x,1)), ti(E(x,2)), 1, numel(ids), numel(ids));  A = (A + A') > 0;
end

function lab = i_labels(A, nV)
% one label per vertex; unknown / medial wall / corpus callosum left empty
    keep = ~contains(lower(string(A.Label(:))), ["unknown" "medial" "corpus"]);
    L = string(A.Label(:));  L = L(keep);  [mx, k] = max(double(A.Membership(keep, :)), [], 1);
    lab = strings(nV, 1);  lab(mx > 0) = L(k(mx > 0));
end

function P = i_member(A, nV)
% [nV x U] membership, each vertex in at most one unit, medial wall excluded
    lab = i_labels(A, nV);  [u, ~, k] = unique(lab(lab ~= ""));
    P = sparse(find(lab ~= ""), k, 1, nV, numel(u));
end

function P = i_mni(V, m)
% Brainstorm SCS (m) -> MRI -> MNI (mm) by the linear NCS fit (cs_convert); [] without it
    P = [];
    if isempty(m) || ~isfield(m, 'NCS') || ~isfield(m.NCS, 'R') || isempty(m.NCS.R) || ...
       ~isfield(m, 'SCS') || ~isfield(m.SCS, 'R') || isempty(m.SCS.R), return; end
    mri = m.SCS.R \ (double(V)' - m.SCS.T(:) / 1000);
    P = 1000 * (m.NCS.R * mri + m.NCS.T(:) / 1000)';
end

function N = i_normals(Sf, gv, V, Fc)
    if isfield(Sf, 'VertNormals') && ~isempty(Sf.VertNormals), N = double(Sf.VertNormals(gv, :));  return; end
    N = zeros(size(V));  fn = cross(V(Fc(:,2),:) - V(Fc(:,1),:), V(Fc(:,3),:) - V(Fc(:,1),:), 2);
    for j = 1:3, N = N + accumarray([repmat(Fc(:,j), 3, 1) kron((1:3)', ones(size(Fc,1), 1))], fn(:), size(V)); end
    N = N ./ vecnorm(N, 2, 2);
end

function u = i_unit(X)
    u = double(X) ./ vecnorm(double(X), 2, 2);
end

% Author: Diellor Basha, 2026
