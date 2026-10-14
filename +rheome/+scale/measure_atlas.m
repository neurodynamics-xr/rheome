function [T, X, Fig] = measure_atlas(name, S, item, opts)
% SCALE.MEASURE_ATLAS  MS1 group P8 (G12) on one participant: correspondence, atlas events, the shared-gauge
% alpha tensor and the connectome roll-up -- the atlas, gauge and multimodal checks, per participant.
%
%   [T, X] = rheome.scale.measure_atlas(name, S, "correspondence")
%   [T, X] = rheome.scale.measure_atlas(name, S, "atlasevents", NSurr=20)
%   [T, X] = rheome.scale.measure_atlas(name, S, "gaugetensor")
%   [T, X] = rheome.scale.measure_atlas(name, S, "connectome", Gammas=[0 0.01 0.1 1 10])
%   [T, X, Fig] = rheome.scale.measure_atlas(name, S, "atlasone")    % one person's atlas, depths 0-6 (figure data)
%   [T, X, Fig] = rheome.scale.measure_atlas(name, S, "xlevel")      % F23: one pattern across the levels
%
% Spec: the group-analysis plan, G12. The geometry parts
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
% The group gauge is rheome.geom.sphereframe: the sphere's meridian pushed to the cortex.
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
%   connectome      the multimodal connectome per participant. (a) Destrieux -> DK against the exact dyadic roll-up (depth 6 ->
%                   5, the participant's own rheome.geom.tree), as F11 C: area outside the parent, share of
%                   finer units > 5 % outside, edge-weight loss ||Ag' C_fine Ag - C_coarse||_F / ||C_coarse||_F
%                   and the CV of the finer units' areas, on the alpha co-power connectome of the kernel and,
%                   with fibres, on the fibre-count connectome. (b) With fibres (rheome.connectome.resolve: the
%                   participant's OWN tractography first, else the HCP-1065 template): the connectome wavelet
%                   (mexhat on the LB-Connectome pencil) per member at Gammas x the auto gamma, over NSeeds
%                   seeds: area-weighted RMS width (mm, Euclidean from the seed), energy in the other
%                   hemisphere, energy beyond 3 sigma (VR1 2.1d), and the coupling scale 2 pi / sqrt(lambda_2).
%   atlasone        MS1 3.4 (F8, F12 on the paper's kernel): one person's alpha atlas on ONE tree, the group
%                   tree of rheome.scale.template, depths 0..LeafDepth. Alpha power is computed once, on the
%                   current summed over each LEAF tile (depth LeafDepth), and every coarser node is the exact
%                   sum of its leaves -- so the levels are one measure rolled up, not a separate inversion per
%                   depth (atlasevents sums the current per depth instead). dB against each node's own median;
%                   events as atlasevents (ThresholdDB, touching tiles, consecutive time tiles) at each depth,
%                   against NSurr phase-randomised surrogates of the leaf currents. Rows per hemisphere x depth
%                   as atlasevents; Fig: per hemisphere the dB tower, the event label of every cell, the
%                   event catalogue (onset, duration, tiles, peak), tile adjacency, codes, centres and areas.
%   xlevel          MS1 2.4, F23 re-run on that tree (left hemisphere, depths 0..LeafDepth). (i) PLANTED: a
%                   heat-kernel peak carried 150 deg down a meridian of the participant's registered sphere
%                   (PlantAzimuthDeg) in NTiles time tiles while it widens 7 -> 18 deg (rheome.flow.movingpeak,
%                   rheome.filters.sphere_kernel), amplitude 9 sin^2 over background 1, rolled up leaf -> root.
%                   (ii) REAL: the atlasone alpha tower, NTiles time tiles from 2 before the onset of the
%                   largest depth-TrackLevel event (or from RealT0S). Both tracked by rheome.detect.tilepath at
%                   every depth (Span "full"); the chain is the TrackLevel track, its ancestors, and the stronger
%                   child below. Leaves also in pole-to-pole spiral order (rheome.geom.spiralorder). Fig has the
%                   fields of the F23 cache (xlevel_planted / xlevel_real).
%
% ⚠ ONE TREE. atlasone and xlevel reach the tree only through rheome.scale.template (rheome.geom.tree on the
% template), so whichever partition geom.tree builds is the one every level, participant and figure uses.
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
        item (1,1) string {mustBeMember(item, ["correspondence" "atlasevents" "gaugetensor" "connectome" "atlasone" "xlevel"])} = "correspondence"
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
        opts.LeafDepth (1,1) double {mustBeInteger, mustBePositive} = 6
        opts.TrackLevel (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.NTiles (1,1) double {mustBeInteger, mustBePositive} = 24
        opts.RealT0S (1,1) double = NaN
        opts.PlantAzimuthDeg (1,1) double = 0
    end
    Fig = [];
    if isempty(S), S = rheome.scale.sensors(name); end
    dmax = max([opts.Depths opts.TensorDepth]);
    switch item
        case "correspondence", Tp = rheome.scale.template(dmax, Template=opts.Template, Hemis=opts.Hemis);  [T, X] = i_correspondence(name, S, Tp, opts);
        case "atlasevents",    Tp = rheome.scale.template(dmax, Template=opts.Template, Hemis=opts.Hemis);  [T, X] = i_events(name, S, Tp, opts);
        case "gaugetensor",    Tp = rheome.scale.template(dmax, Template=opts.Template, Hemis=opts.Hemis);  [T, X] = i_tensor(name, S, Tp, opts);
        case "connectome",     [T, X] = i_connectome(name, S, opts);
        case "atlasone",       Tp = rheome.scale.template(opts.LeafDepth, Template=opts.Template, Hemis=opts.Hemis);  [T, X, Fig] = i_atlasone(name, S, Tp, opts);
        case "xlevel",         o = opts;  o.Hemis = opts.Hemis(1);  o.NSurr = 0;
                               Tp = rheome.scale.template(opts.LeafDepth, Template=opts.Template, Hemis=o.Hemis);  [T, X, Fig] = i_xlevel(name, S, Tp, o);
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
            [P, A] = rheome.scale.grouptiles(S, Tp, h, l);  r0 = size(Y, 1);
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

function [n, dur, ntiles, lab] = i_events_of(on, A, dt)
% connected sets of on-cells: touching tiles within a time tile, same or touching tile across consecutive ones
    [nT, nb] = size(on);  id = reshape(1:nT*nb, nT, nb);
    % ⚠ e(:), t(:): with one edge (depth 1, two tiles) on(i1,:) is a row and find returns rows (MATLAB:catenate)
    [i1, j1] = find(triu(A));  [e, t] = find(on(i1, :) & on(j1, :));  e = e(:);  t = t(:);
    s = id(sub2ind([nT nb], i1(e), t));  d = id(sub2ind([nT nb], j1(e), t));
    [i2, j2] = find(A | speye(nT));  [e, t] = find(on(i2, 1:end-1) & on(j2, 2:end));  e = e(:);  t = t(:);
    s = [s; id(sub2ind([nT nb], i2(e), t))];  d = [d; id(sub2ind([nT nb], j2(e), t + 1))];
    c = conncomp(graph(sparse(s, d, 1, nT*nb, nT*nb) + sparse(d, s, 1, nT*nb, nT*nb)));
    v = find(on(:));  [~, ~, g] = unique(c(v));  n = max([g; 0]);
    lab = zeros(nT, nb, 'uint32');  lab(v) = g;
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
        [P, ~, ids] = rheome.scale.grouptiles(S, Tp, h, o.TensorDepth);
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

%% ---------- one person's atlas on one tree (MS1 3.4) and one pattern across its levels (F23, MS1 2.4)
function [T, X, Fig] = i_atlasone(name, S, Tp, o)
    [F, fs, t0] = i_rec(name, S);  mins = size(F, 2) / fs / 60;
    [bp, nb] = i_band(o, fs);  dt = nb / fs;
    T = rheome.scale.rows("atlasone", strings(0,1), [], "");  X = table();
    Fig = struct('band_hz', o.Band, 'threshold_db', o.ThresholdDB, 'dt', dt, 't_offset_s', t0, 'leaf_depth', o.LeafDepth, ...
                 'kernel', "unconstrained noise-whitened minimum norm (S.Res), current summed over each leaf tile", ...
                 'tree', "rheome.scale.template: rheome.geom.tree on the template, carried through the registered sphere");
    for h = o.Hemis(:)'
        tw = i_leaftower(S, Tp, h, o, F);
        dB = i_db(tw, i_leafpower(tw.Y, bp, nb));
        Su = cell(1, o.NSurr);
        for s = 1:o.NSurr, rng(o.Seed + s);  Su{s} = i_db(tw, i_leafpower(i_phaserand(tw.Y), bp, nb)); end
        lab = cell(1, tw.D + 1);  cg = table();
        for l = 0:tw.D
            [d, lab{l+1}] = i_stats(dB{l+1}, tw.A{l+1}, o.ThresholdDB, dt, mins);
            sr = zeros(o.NSurr, 6);  for s = 1:o.NSurr, sr(s, :) = i_stats(Su{s}{l+1}, tw.A{l+1}, o.ThresholdDB, dt, mins); end
            r = h + ".d" + l;
            T = [T; rheome.scale.rows("atlasone", ...
                 ["events_per_min" "duration_median_s" "tiles_visited_median" "frac_on" ...
                  "surr_events_per_min" "surr_duration_median_s" "surr_tiles_visited_median" "surr_frac_on" "p_rate"], ...
                 [d(2:5) mean(sr(:, 2:5), 1, 'omitnan') (1 + nnz(sr(:, 2) >= d(2))) / (1 + o.NSurr)], ...
                 ["1/min" "s" "tiles" "fraction" "1/min" "s" "tiles" "fraction" "p"], r)]; %#ok<AGROW>
            M = [d; sr];  src = ["data"; repmat("surrogate", o.NSurr, 1)];
            X = [X; table(repmat(h, o.NSurr + 1, 1), repmat(l, o.NSurr + 1, 1), src, (0:o.NSurr)', ...
                 M(:,1), M(:,2), M(:,3), M(:,4), M(:,5), M(:,6), 'VariableNames', ...
                 {'hemi','depth','source','surrogate','n_events','events_per_min','duration_median_s','tiles_visited_median','frac_on','duration_p90_s'})]; %#ok<AGROW>
            cg = [cg; i_catalogue(dB{l+1}, lab{l+1}, tw.ids{l+1}, h, l, dt, t0)]; %#ok<AGROW>
            fprintf('[atlasone %s] %s: %.1f events/min (surrogates %.1f)\n', name, r, d(2), mean(sr(:, 2)));
        end
        Fig.(h) = struct('db', {cellfun(@single, dB, 'UniformOutput', false)}, 'event', {lab}, 'catalogue', cg, ...
                         'ids', {tw.ids}, 'adjacency', {tw.A}, 'centre', {tw.C}, 'area_m2', {tw.area});
    end
end

function [T, X, Fig] = i_xlevel(name, S, Tp, o)
    D = o.LeafDepth;  L4 = o.TrackLevel;  NT = o.NTiles;  h = o.Hemis(1);
    if L4 > D, error('scale:atlas:xlevel', 'TrackLevel %d is below LeafDepth %d', L4, D); end
    SWEEP = 150;  COLAT0 = 15;  SIG = [7 18];  AMP = 9;  NWAY = 151;  LMAX = 120;  TURNS = [1 1.5 2 3 4 6 8];
    [F, fs, t0] = i_rec(name, S);  [bp, nb] = i_band(o, fs);  DT = nb / fs;
    tw = i_leaftower(S, Tp, h, o, F);  clear F
    for l = 0:D
        if ~isequal(tw.ids{l+1}(:)', 2^l:2^(l+1)-1), error('scale:atlas:heap', 'depth %d: tile codes are not the full heap', l); end
    end
    a = tw.a;  P6 = tw.P{D+1};  gv = S.B.(h).gv(:);
    u = double(S.Sf.Sphere(gv, :));  R = mean(vecnorm(u, 2, 2));  u = u ./ vecnorm(u, 2, 2);

    % (i) PLANTED: a widening heat-kernel peak down a meridian of the registered sphere
    m = [cosd(o.PlantAzimuthDeg) sind(o.PlantAzimuthDeg) 0];
    c0 = cosd(COLAT0)*[0 0 1] + sind(COLAT0)*m;  b = -sind(COLAT0)*[0 0 1] + cosd(COLAT0)*m;
    phw = linspace(0, deg2rad(SWEEP), NWAY)';  sigw = deg2rad(linspace(SIG(1), SIG(2), NWAY))';
    atoms = zeros(numel(gv), NWAY);
    for k = 1:NWAY
        c = cos(phw(k))*c0 + sin(phw(k))*b;
        at = rheome.filters.sphere_kernel(acos(max(-1, min(1, u*c'))), exp(-(0:LMAX)'.*(1:LMAX+1)' * sigw(k)^2/2));
        atoms(:, k) = at / max(at);
    end
    cums = phw * R;  v = cums(end) / ((NT-1) * DT) * (1 - 1e-9);          % (1-1e-9): exactly NT frames
    Pk = rheome.flow.movingpeak(atoms, cums, SpeedMS=v, SampleRate=1/DT);
    if size(Pk.X, 2) ~= NT, error('scale:atlas:xlevel', 'movingpeak gave %d frames, not %d', size(Pk.X, 2), NT); end
    env = sin(pi * ((0:NT-1) + 0.5) / NT).^2;  Xp = 1 + AMP * env .* Pk.X;
    phi = Pk.s / R;  centre = cos(phi)*c0 + sin(phi)*b;
    sigma = rad2deg(sigw(Pk.k0) .* (1 - Pk.w) + sigw(Pk.k0 + 1) .* Pk.w);
    leaf = P6' * (a .* Xp);  la = P6' * a;  err = 0;  plant = cell(D+1, 1);
    for l = 0:D
        up = tw.Ru{l+1} * leaf;  direct = tw.P{l+1}' * (a .* Xp);           % leaf roll-up vs the level's own tiles
        err = max(err, max(abs(up - direct), [], 'all') / max(abs(direct), [], 'all'));
        plant{l+1} = 10*log10(up ./ (tw.Ru{l+1} * la));
    end
    [~, vtrue] = max(u * centre', [], 1);  truth = zeros(D+1, NT);
    for l = 0:D, [~, ti] = max(tw.P{l+1}(vtrue, :), [], 2);  truth(l+1, :) = ti'; end

    % spiral order of the leaves
    A6 = tw.A{D+1};  cen = (P6' * (a .* u)) ./ la;  cen = cen ./ vecnorm(cen, 2, 2);  [ei, ej] = find(triu(A6));
    gap = zeros(size(TURNS));  r = zeros(1, 2^D);
    for i = 1:numel(TURNS), oi = rheome.geom.spiralorder(cen, Turns=TURNS(i));  r(oi) = 1:2^D;  gap(i) = mean(abs(r(ei) - r(ej))); end
    [~, ib] = min(gap);  turns = TURNS(ib);  [order, tspiral] = rheome.geom.spiralorder(cen, Turns=turns);
    rs = zeros(1, 2^D);  rs(order) = 1:2^D;  rt = 1:2^D;
    mono = @(r) [corr((1:NT)', r(truth(D+1,:))', 'Type', 'Spearman'), sum(diff(r(truth(D+1,:))) < 0), max(abs(diff(r(truth(D+1,:)))))];
    ms = mono(rs);  mt = mono(rt);

    % (ii) REAL: the alpha tower, NT time tiles around the largest event at TrackLevel
    dB = i_db(tw, i_leafpower(tw.Y, bp, nb));  nbin = size(dB{1}, 2);
    if isfinite(o.RealT0S)
        k0 = floor((o.RealT0S - t0) / DT) + 1;
    else
        [~, ~, ~, lab] = i_events_of(dB{L4+1} > o.ThresholdDB, tw.A{L4+1}, DT);
        k0 = 1;
        if any(lab(:))
            [~, big] = max(accumarray(double(lab(lab > 0)), 1));  [~, bins] = find(lab == big);  k0 = min(bins) - 2;
        end
    end
    k0 = max(1, min(k0, nbin - NT + 1));
    rl = cellfun(@(y) y(:, k0 + (0:NT-1)), dB, 'UniformOutput', false);

    EX = struct('name', {"planted", "real"}, 'levels', {plant, rl}, 'thr', {3, o.ThresholdDB});
    X = table();
    for e = 1:2
        Y = EX(e).levels;  own = ones(D+1, NT);  am = ones(D+1, NT);  chain = ones(D+1, NT);
        for l = 1:D
            tr = rheome.detect.tilepath(Y{l+1}, struct('A', tw.A{l+1}), Span="full", Baseline=EX(e).thr);
            own(l+1, :) = tr.tracks(1).tiles';  [~, am(l+1, :)] = max(Y{l+1}, [], 1);
        end
        chain(L4+1, :) = own(L4+1, :);
        for l = L4-1:-1:0, chain(l+1, :) = floor((chain(l+2, :) - 1) / 2) + 1; end
        for l = L4+1:D
            kids = 2*chain(l, :) - 1 + [0; 1];
            [~, w] = max([Y{l+1}(sub2ind(size(Y{l+1}), kids(1,:), 1:NT)); Y{l+1}(sub2ind(size(Y{l+1}), kids(2,:), 1:NT))]);
            chain(l+1, :) = kids(sub2ind(size(kids), w, 1:NT));
        end
        nj = @(p, l) sum(arrayfun(@(k) p(k) ~= p(k+1) && ~tw.A{l+1}(p(k), p(k+1)), 1:NT-1));
        EX(e).own = own;  EX(e).argmax = am;  EX(e).chain = chain;  EX(e).argmax_jumps = arrayfun(@(l) nj(am(l+1,:), l), 0:D);
        EX(e).on_tiles = sum(Y{L4+1}(sub2ind(size(Y{L4+1}), chain(L4+1,:), 1:NT)) > EX(e).thr);
        [lv, fr] = ndgrid(0:D, 1:NT);  val = arrayfun(@(l, f) Y{l+1}(chain(l+1, f), f), lv, fr);
        tru = nan(D+1, NT);  if e == 1, tru = truth; end
        X = [X; table(repmat(EX(e).name, numel(lv), 1), lv(:), fr(:), chain(:), own(:), am(:), tru(:), val(:), ...
             'VariableNames', {'example','level','frame','chain','track','argmax','truth','db_chain'})]; %#ok<AGROW>
    end
    hit = mean(EX(1).chain == truth, 2)';
    T = rheome.scale.rows("xlevel", ["rollup_rel_err" "planted_hit_track_level" "planted_hit_leaf" "spiral_rho" "tree_rho" ...
            "planted_leaf_argmax_jumps" "real_leaf_argmax_jumps" "real_track_level_argmax_jumps" "real_on_tiles" "real_t0_s"], ...
        [err hit(L4+1) hit(D+1) ms(1) mt(1) EX(1).argmax_jumps(D+1) EX(2).argmax_jumps(D+1) EX(2).argmax_jumps(L4+1) ...
         EX(2).on_tiles t0 + (k0-1)*DT], ...
        ["relative" "fraction" "fraction" "rho" "rho" "jumps" "jumps" "jumps" "tiles" "s"], h);
    fprintf('[xlevel %s] roll-up %.1e; planted l%d hit %.2f; real argmax jumps l%d %d vs l%d %d; on %d/%d\n', name, err, ...
            L4, hit(L4+1), D, EX(2).argmax_jumps(D+1), L4, EX(2).argmax_jumps(L4+1), EX(2).on_tiles, NT);
    Fig = struct('planted', EX(1), 'real', EX(2), 'truth', truth, 'centre', centre, 'sigma_deg', sigma, ...
        'phi_deg', rad2deg(phi), 'env', env(:), 'amp', AMP, 'dt', DT, 't_plant', (0:NT-1)'*DT, ...
        't_real', t0 + (k0 - 1 + (0:NT-1)')*DT, 'speed_sphere_ms', v, 'leaf_order', order, 'leaf_t', tspiral, ...
        'turns', turns, 'turns_tried', TURNS, 'rank_gap', gap, 'rank_gap_tree', mean(abs(rt(ei) - rt(ej))), ...
        'mono_spiral', ms, 'mono_tree', mt, 'leaf_centroid', cen, 'track_level', L4, 'adjacency', {tw.A}, ...
        'hemi', h, 'plant_azimuth_deg', o.PlantAzimuthDeg, 'rollup_rel_err', err, ...
        'kernel', "unconstrained noise-whitened minimum norm (S.Res)", 'tree', "rheome.scale.template group tree");
end

function [F, fs, t0] = i_rec(name, S)
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
    t0 = 0;  if isfield(st.rec, 'Time') && ~isempty(st.rec.Time), t0 = st.rec.Time(1); end
end

function [bp, nb] = i_band(o, fs)
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', o.Band(1), ...
                    'HalfPowerFrequency2', o.Band(2), 'SampleRate', fs);
    nb = round(o.TimeTileS * fs);
end

function tw = i_leaftower(S, Tp, h, o, F)
% the group tiles of hemisphere h at depths 0..LeafDepth, the leaf -> node roll-up, and the leaf currents
    D = o.LeafDepth;  gv = S.B.(h).gv(:);  tw.D = D;  tw.a = full(sum(S.B.(h).lbo.Mass, 2));
    for l = D:-1:0
        [tw.P{l+1}, tw.A{l+1}, tw.ids{l+1}, tw.C{l+1}] = rheome.scale.grouptiles(S, Tp, h, l);
        [~, pj] = ismember(floor(tw.ids{D+1} / 2^(D-l)), tw.ids{l+1});
        tw.Ru{l+1} = sparse(pj, 1:numel(pj), 1, numel(tw.ids{l+1}), numel(pj));
        tw.area{l+1} = tw.P{l+1}' * tw.a;
    end
    K = S.Res.ImagingKernel;  P6 = tw.P{D+1}';  tw.Y = [];
    for k = 1:3, tw.Y = [tw.Y; (P6 * K((gv - 1)*3 + k, :)) * F]; end
end

function pw = i_leafpower(Y, bp, nb)
% band power of the leaf currents (x + y + z), mean per time tile: [nLeaf x nBins]
    E = abs(hilbert(filtfilt(bp, Y.')).').^2;  n = size(Y, 1) / 3;  nbin = floor(size(E, 2) / nb);
    pw = E(1:n, :) + E(n+1:2*n, :) + E(2*n+1:end, :);
    pw = reshape(mean(reshape(pw(:, 1:nbin*nb), n, nb, nbin), 2), n, nbin);
end

function dB = i_db(tw, pw6)
% every node is the exact sum of its leaves; dB against the node's own median (a node no vertex reaches: NaN)
    dB = cell(1, tw.D + 1);
    for l = 0:tw.D
        p = tw.Ru{l+1} * pw6;  p(p <= 0) = NaN;  dB{l+1} = 10*log10(p ./ median(p, 2, 'omitnan'));
    end
end

function [d, lab] = i_stats(dB, A, thr, dt, mins)
% [n_events events_per_min duration_median_s tiles_visited_median frac_on duration_p90_s]
    on = dB > thr;  [n, dur, nt, lab] = i_events_of(on, A, dt);
    d = [n, n / mins, median(dur), median(nt), mean(on, 'all'), prctile(dur, 90)];
    if n == 0, d([3 4 6]) = NaN; end
end

function C = i_catalogue(dB, lab, ids, h, l, dt, t0)
% one row per event: onset, duration, tiles visited, peak dB and the tile code holding it
    [tile, bin] = find(lab);  tile = tile(:);  bin = bin(:);  g = double(lab(lab > 0));  g = g(:);   % (:): one tile gives rows
    if isempty(g)
        C = table(strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), 'VariableNames', ...
            {'hemi','depth','event','onset_s','duration_s','n_tiles','peak_db','peak_code'});  return
    end
    v = dB(sub2ind(size(dB), tile, bin));  v = v(:);  n = max(g);
    pk = accumarray(g, v, [], @max);  ip = arrayfun(@(e) find(g == e & v == pk(e), 1), (1:n)');
    on = accumarray(g, bin, [], @min);  off = accumarray(g, bin, [], @max);
    nt = accumarray(g, tile, [], @(x) numel(unique(x)));
    C = table(repmat(h, n, 1), repmat(l, n, 1), (1:n)', t0 + (on - 1) * dt, (off - on + 1) * dt, nt, pk, ids(tile(ip)), ...
        'VariableNames', {'hemi','depth','event','onset_s','duration_s','n_tiles','peak_db','peak_code'});
end

%% ---------- shared helpers
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
