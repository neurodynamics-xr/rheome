function [T, X, Cn] = measure_multimodal(name, S, opts)
% SCALE.MEASURE_MULTIMODAL  MEG, PET and the fibre connectome through ONE surface wavelet bank on ONE set of
% dyadic tiles, rolled up exactly across levels -- the MS1 multimodal worked example.
%
%   [T, X, Cn] = rheome.scale.measure_multimodal(name)
%   [T, X, Cn] = rheome.scale.measure_multimodal(name, S, Depth=6)
%   rheome.scale.run(name, ProtoDir=..., Sub=..., Analyses="multimodal")   % writes the three outputs below
%
% A descriptive method demonstration, not a test: no hypothesis, no p-value.
%
% THE SHARED REPRESENTATION. Every field is read on the participant's own cortex (tess_cortex_pial_low, the
% leadfield's surface, which Brainstorm's PET projection also uses) by the same mexhat rheome.graphfilterbank as
% rheome.scale.coefficients (one bank on the smaller hemisphere lambda_max, usable wavelet members only), plus the
% unfiltered field as member 0. Fields:
%   MEG    per band (Bands): the band power at each vertex of the wavelet-filtered current. The wavelet acts on
%          the source current (linear: Phi g_j Phi' M K, per orientation), power comes after:
%          p_j(v) = sum_xyz diag(X_j C_b X_j'), C_b the sensors' band cross-spectrum (Hilbert, real part),
%          K the paper's kernel S.Res (unconstrained whitened minimum norm, SnrFixed 3).    unit (A m)^2
%   PET    per tracer (pet.mat, cached by rheome.scale.importsubject from the protocol's PET study): the SUVR map
%          and its wavelet coefficients c_j = Phi g_j Phi' M x (signed).                    unit SUVR
%   fibres the endpoint-degree density of the participant's own tractography (fibres ending at a vertex per
%          m^2 of cortex) and its wavelet coefficients, as PET.                              unit fibres/m^2
% THE TILES are rheome.geom.tree on each hemisphere to Depth (default 6: 64 tiles per hemisphere), numbered
% with rheome.scale.coefficients' ladder codes (root 1, L 2, R 3, children 2c and 2c+1; finest level
% Depth+1). Only additive statistics are stored at the finest tile -- area, sum(area*x), sum(area*x^2) per
% field x member, and the fibre-count connectome between tiles -- and every coarser level is their sum.
%
% ⭐ EXACTNESS IS MEASURED, NOT ASSUMED. At every depth the roll-up is compared with the same sums taken
% directly on that depth's tiles (rheome.geom.tiles); rollup_max_relerr and connectome_rollup_max_relerr are
% the worst relative differences (round-off when the tree is a heap).
% ⚠ A SIGNED BAND-PASS COEFFICIENT AVERAGES TOWARD ZERO over a tile larger than its scale
% (rheome.geom.tilemean), so for PET and fibre members j > 0 compare tiles on ms = sum(a x^2)/sum(a), the
% coefficient energy; MEG members are powers already (compare on mean).
% ⚠ RAW MEG POWER ACROSS TILES IS MOSTLY THE KERNEL'S DEPTH SENSITIVITY: in one participant the alpha
% tile power spans ~1000x while delta/theta varies +-12 %, so rho between raw band powers (and between a raw
% power and anything deep-vs-superficial) reads the instrument. Compare a band relative to the sum of the bands
% (exact at any level from the stored sums) before reading MEG against PET or fibres.
% ⚠ The connectome is RAW endpoint counts between tiles (each end at its nearest cortex vertex, the step
% rheome.operators.connectome starts from), without its 3-hop vertex smoothing: at tile scale the smoothing
% only moves weight across a tile border, and the smoothed vertex matrix costs ~20 GB at 1e6 fibres.
% ⚠ PET is SUVR as the protocol carries it (bst-pet: no partial-volume correction unless the protocol's
% PET study is the PVC one); the representation does not depend on which.
%
% OUTPUTS
%   T   metric rows (rheome.scale.rows): tiles, members, fields, roll-up errors, fibres, tracers
%   X.tiles  one row per level x tile x field x member: level, code, hemisphere, field, modality, unit,
%            member, scale_mm (sigma = sqrt(2t), 0 = unfiltered), area_m2, sum_ax, sum_ax2, mean, ms
%   X.assoc  descriptive Spearman rho across the tiles of each level >= 3, per member, between the fields'
%            tile values (mean for MEG and for member 0, ms otherwise); n_tiles
%   Cn  struct: codes{k} and C{k} the tile connectome (fibre counts) at ladder level k-1, plus provenance
%
% See also: rheome.scale.coefficients, rheome.scale.measure_atlas, rheome.geom.tree, rheome.geom.tiles,
%           rheome.graphfilterbank, rheome.scale.importsubject, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depth (1,1) double {mustBeInteger, mustBePositive} = 6
        opts.Bands = struct('name', {"delta" "theta" "alpha" "beta"}, 'hz', {[1 4] [4 8] [8 13] [13 30]})
        opts.AssocMinLevel (1,1) double = 3
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    D = opts.Depth;  hemis = ["L" "R"];  d = fullfile(rheome.load.root(), name);
    K = S.Res.ImagingKernel;  nV = S.Sf.nV;

    gfb = rheome.graphfilterbank(min(arrayfun(@(h) max(S.B.(h).lbo.Lambda), hemis)));
    t = scales(gfb);  w = widths(gfb);  keep = find(gfb.Usable & isfinite(t) & t > 0);
    nJ = numel(keep) + 1;  scaleMM = [0, 1e3 * w(keep)];

    % ---- fields: MEG band cross-spectra, PET maps, fibre endpoints
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = S.Ref * double(st.rec.F(S.iSel, :));  clear st
    nb = numel(opts.Bands);  Cs = cell(1, nb);
    for b = 1:nb
        bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', opts.Bands(b).hz(1), ...
                        'HalfPowerFrequency2', opts.Bands(b).hz(2), 'SampleRate', fs);
        Fa = hilbert(filtfilt(bp, F.')).';  Cs{b} = real(Fa * Fa') / size(Fa, 2);  clear Fa
    end
    clear F
    pet = struct('name', {}, 'x', {}, 'units', {}, 'file', {});
    if isfile(fullfile(d, 'pet.mat')), pet = getfield(load(fullfile(d, 'pet.mat'), 'pet'), 'pet'); end
    for p = 1:numel(pet)
        if numel(pet(p).x) ~= nV
            error('scale:multimodal:pet', 'PET %s has %d values, the cortex %d vertices', pet(p).name, numel(pet(p).x), nV);
        end
    end
    ff = fullfile(d, 'fibers_subject.mat');  haveFib = isfile(ff);  deg = zeros(nV, 1);  nFib = 0;  nSelf = 0;
    if haveFib                                           % own tractography only: the example is one person's
        ep = rheome.io.read.fibers(ff);  nFib = size(ep, 1);
        v1 = knnsearch(S.Sf.Vertices, squeeze(ep(:, 1, :)));  v2 = knnsearch(S.Sf.Vertices, squeeze(ep(:, 2, :)));  clear ep
        nSelf = sum(v1 == v2);  x = v1 ~= v2;  v1 = v1(x);  v2 = v2(x);
        deg = accumarray([v1; v2], 1, [nV 1]);
    end
    fld = [compose("meg_%s", string({opts.Bands.name})), "pet_" + string({pet.name}), repmat("fibres_degree", 1, haveFib)];
    mdl = [repmat("MEG", 1, nb), repmat("PET", 1, numel(pet)), repmat("fibres", 1, haveFib)];
    uni = [repmat("(A m)^2", 1, nb), string({pet.units}), repmat("fibres/m^2", 1, haveFib)];
    nF = numel(fld);

    % ---- per hemisphere: fields x members at the vertices, sums on the finest tiles, direct sums per depth
    codes = [];  A = [];  S1 = [];  S2 = [];  tileOf = zeros(nV, 1);  errDirect = 0;  direct = cell(2, D);
    for hi = 1:2
        H = S.B.(hemis(hi));  gv = double(H.gv(:));  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  lam = H.lbo.Lambda(:);
        a = full(sum(Mm, 2));  nVh = numel(gv);
        g = zeros(numel(lam), numel(keep));
        for j = 1:numel(keep), gm = gain(gfb, keep(j)); g(:, j) = gm(lam); end
        V = zeros(nVh, nF, nJ);
        rows = reshape((gv' - 1) * 3 + (1:3)', [], 1);
        for o3 = 1:3                                     % MEG: wavelet on the current, power after
            Kd = K(rows(o3:3:end), :);  Km = Phi' * (Mm * Kd);
            for b = 1:nb, V(:, b, 1) = V(:, b, 1) + sum((Kd * Cs{b}) .* Kd, 2); end
            for j = 1:numel(keep)
                Xj = Phi * (g(:, j) .* Km);
                for b = 1:nb, V(:, b, j + 1) = V(:, b, j + 1) + sum((Xj * Cs{b}) .* Xj, 2); end
            end
        end
        sc = [arrayfun(@(p) double(pet(p).x(gv)), 1:numel(pet), 'UniformOutput', false), repmat({deg(gv) ./ a}, 1, haveFib)];
        for k = 1:numel(sc)                              % PET and fibre degree: wavelet on the map itself
            V(:, nb + k, 1) = sc{k};  c = Phi' * (Mm * sc{k});
            for j = 1:numel(keep), V(:, nb + k, j + 1) = Phi * (g(:, j) .* c); end
        end

        Tr = rheome.geom.tree(H.S, L=H.lbo.L, M=Mm, MaxDepth=D);
        G = rheome.geom.tiles(Tr, H.S, D);
        if ~G.heap || any(G.depth ~= D)
            error('scale:multimodal:ladder', '%s hemisphere: the tree stopped above depth %d. Use a smaller Depth.', hemis(hi), D);
        end
        [ch, o] = sort(G.node_id + 2^D * hi);  P = double(G.P(:, o));
        [~, tix] = max(P, [], 2);  tileOf(gv) = numel(codes) + tix;
        Vf = reshape(V, nVh, []);
        codes = [codes; ch];  A = [A; P' * a];  S1 = [S1; P' * (a .* Vf)];  S2 = [S2; P' * (a .* Vf.^2)]; %#ok<AGROW>
        for l = 1:D                                      % the same sums straight from depth l's own tiles
            Gl = rheome.geom.tiles(Tr, H.S, l);  [cl, ol] = sort(Gl.node_id + 2^l * hi);  Pl = double(Gl.P(:, ol));
            direct{hi, l} = struct('codes', cl, 'A', Pl' * a, 'S1', Pl' * (a .* Vf), 'P', Pl, 'gv', gv);
        end
    end
    if any(tileOf == 0), error('scale:multimodal:cover', '%d vertices in no tile', sum(tileOf == 0)); end

    % ---- roll-up to every ladder level, checked against the direct sums
    nT = numel(codes);  L1 = D + 1;  tiles = table();  Cn = struct('codes', {cell(1, L1 + 1)}, 'C', {cell(1, L1 + 1)});
    Cfine = sparse(nT, nT);
    if haveFib, Cfine = sparse(tileOf([v1; v2]), tileOf([v2; v1]), 1, nT, nT); end
    errConn = 0;
    for k = 0:L1
        anc = floor(codes / 2^(L1 - k));  [uc, ~, gi] = unique(anc);  Ag = sparse(1:nT, gi, 1, nT, numel(uc));
        Ak = Ag' * A;  S1k = Ag' * S1;  S2k = Ag' * S2;  Ck = Ag' * Cfine * Ag;
        Cn.codes{k + 1} = uc;  Cn.C{k + 1} = full(Ck);
        if k >= 2                                        % levels 2..L1 are depths 1..D of the hemisphere trees
            for hi = 1:2
                dr = direct{hi, k - 1};  [~, ix] = ismember(dr.codes, uc);
                errDirect = max([errDirect, max(abs(Ak(ix) - dr.A)) / max(dr.A), ...
                                 max(abs(S1k(ix, :) - dr.S1), [], 'all') / max(abs(dr.S1), [], 'all')]);
            end
            if haveFib                                   % direct tile connectome at this depth
                tl = zeros(nV, 1);
                for hi = 1:2
                    dr = direct{hi, k - 1};  [~, tx] = max(dr.P, [], 2);  [~, ix] = ismember(dr.codes, uc);  tl(dr.gv) = ix(tx);
                end
                Cd = sparse(tl([v1; v2]), tl([v2; v1]), 1, numel(uc), numel(uc));
                errConn = max(errConn, full(max(abs(Cd - Ck), [], 'all')) / max(1, full(max(Cd, [], 'all'))));
            end
        end
        hemi = zeros(numel(uc), 1);  if k >= 1, hemi = floor(uc / 2^(k - 1)) - 1; end
        [ti, fj] = ndgrid(1:numel(uc), 1:nF * nJ);  fi = mod(fj - 1, nF) + 1;  mj = floor((fj - 1) / nF);
        s1 = S1k(:);  s2 = S2k(:);  ar = Ak(ti(:));
        tiles = [tiles; table(repmat(k, numel(ti), 1), uc(ti(:)), hemi(ti(:)), fld(fi(:))', mdl(fi(:))', uni(fi(:))', ...
                 mj(:), scaleMM(mj(:) + 1)', ar, s1, s2, s1 ./ ar, s2 ./ ar, 'VariableNames', ...
                 {'level','code','hemisphere','field','modality','unit','member','scale_mm','area_m2','sum_ax','sum_ax2','mean','ms'})]; %#ok<AGROW>
    end
    Cn.level = 0:L1;  Cn.n_fibres = nFib;  Cn.n_self = nSelf;  Cn.depth = L1;

    % ---- descriptive association across tiles (no test)
    assoc = table();
    for k = max(opts.AssocMinLevel, 0):L1
        for m = 0:nJ - 1
            Y = zeros(2^k, nF);
            for f = 1:nF
                r = tiles(tiles.level == k & tiles.member == m & tiles.field == fld(f), :);
                if mdl(f) == "MEG" || m == 0, Y(:, f) = r.mean; else, Y(:, f) = r.ms; end
            end
            R = corr(Y, 'Type', 'Spearman');
            [i1, i2] = find(triu(true(nF), 1));
            assoc = [assoc; table(repmat(k, numel(i1), 1), repmat(m, numel(i1), 1), repmat(scaleMM(m + 1), numel(i1), 1), ...
                     fld(i1)', fld(i2)', R(sub2ind([nF nF], i1, i2)), repmat(2^k, numel(i1), 1), ...
                     'VariableNames', {'level','member','scale_mm','field_a','field_b','rho','n_tiles'})]; %#ok<AGROW>
        end
    end
    X = struct('tiles', tiles, 'assoc', assoc);

    hl = zeros(nV, 1);  hl(double(S.B.L.gv)) = 1;  hl(double(S.B.R.gv)) = 2;  inter = NaN;
    if haveFib, inter = mean(hl(v1) ~= hl(v2)); end
    T = rheome.scale.rows("multimodal", ["n_tiles_finest" "n_levels" "n_members" "n_fields" "n_meg_bands" "n_pet" ...
            "fibres_available" "n_fibres" "fibres_self_tile_vertex" "interhemispheric_fraction" ...
            "rollup_max_relerr" "connectome_rollup_max_relerr" "connectome_total"], ...
            [nT L1 + 1 nJ nF nb numel(pet) haveFib nFib nSelf inter errDirect errConn full(sum(Cfine, 'all'))], ...
            ["tiles" "levels" "members" "fields" "bands" "tracers" "flag" "fibres" "fibres" "fraction" "relative" "relative" "fibre ends"], ...
            strjoin(fld, " "));
end

% Author: Diellor Basha, 2026
