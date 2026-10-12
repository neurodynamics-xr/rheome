function [Z, Y, Op] = jointfields(name, S, opts)
% MODES.JOINTFIELDS  One participant's joint graph-time wavelet field W(v,s,t,f) on group tiles, its nulls and T(s).
%
%   [Z, Y] = rheome.modes.jointfields(name)
%   [Z, Y] = rheome.modes.jointfields(name, S, Depths=0:5, Octaves=[1 2; 2 4; 4 8; 8 16; 16 32], Seed=1)
%   [Z, Y] = rheome.modes.jointfields(name, S, Canary=true)     % J5: planted patterns on F, no real-data series
%
% Protocol section 14 (H4 ATLAS-1), J1, J3, J4 and the J5 plants. Same recording, clean span and work rate as
% rheome.modes.fields; the band-limited field of section 4.1 is untouched.
%
% J1 FIELD. The separable members of the joint filterbank (MD24, rheome.jointfilterbank: graph factor x time
% factor) applied factor by factor -- wt would materialise [K x nOmega x Nf]:
%   graph scale s  one dyadic depth d of the group tree (Depths, 0-5). The tight logitersine bank's members
%                  (2 voices, as rheome.modes.operator) whose wavelength lies in the depth's window are kept:
%                  windows are cut at the geometric midpoints between the depths' tile diameters (the
%                  participant's median tile area, 2 sqrt(A/pi)); the coarsest depth is open above. The scale
%                  is pooled to the group tiles AT ITS OWN DEPTH (area-weighted mean): one row per tile x scale
%                  (Z.tiles{1}: hemi, tile, depth, scale_mm, ...; tile codes are unique across depths)
%   octave f       the jointfilterbank's raised-cosine band member for each row of Octaves, applied to the
%                  one-sided spectrum: an analytic signal, decimated to the largest rate fs/m >= 3 x the upper
%                  edge (rheome.modes.octave). valid as rheome.modes.analytic
% The field is computed on the PERIODIC part (MS1 rule 8): each tile x scale row is split by its own Welch
% background fit (rheome.modes.split) before the octave.
%
% J3 NULLS (each through the participant's own gain, kernel, bank and split):
%   f_fieldspread  F, field spread: independent normal sources at every vertex, one spectrum per generator
%                  tile (depth GenDepth = 6), variance ~ area, plus the participant's own empty room on the
%                  channels it shares. The generator spectra are the nonnegative least-squares solution that
%                  makes F's tile power spectra at GenDepth equal the data's minus the empty room's
%                  (M(i,j) = W_i C_j W_i', C_j the sensor covariance of tile j's sources: Op.Cgen), per
%                  0.5 Hz bin; then one scalar recalibration of total power. So each location keeps its own
%                  spectrum and power, and every coupling between locations is removed (what remains is the
%                  instrument's spread). Z.fspread: fit error, recalibration factor, empty room used
%   n4_aperiodic   N4: the aperiodic complement of each row's split
%   n5_phaserand   N5: the periodic rows with one random phase per frequency bin shared by every row (keeps
%                  every auto- and cross-spectrum, across scales too)
% N2 (spin) and N3 act on fitted atoms, not here.
%
% J4 RESOLUTION. Z.transfer = rheome.modes.scaletransfer(S, Z.scales): T(s) per depth and the scale edge (mm).
%
% J5 CANARY (Canary=true). No real-data series is written. Four patterns are planted (normal current on the
% participant's cortex, through its gain, then the same pipeline) on F as the background, in turn (states of
% 0.2-0.4 s, 20 ms ramps) at CanaryHz, each, while on, at the sensor power of F in the plant's octave (0 dB):
%   wave200_ap    a 200 mm plane wave travelling A->P (SCS x) over the whole cortex
%   pair_fo       two lagged 30 mm blobs, frontal leading occipital by 90 deg
%   standing200   cos(2 pi z / 200 mm), a standing (zero-lag) dorsal-ventral pattern
%   fine50        a 50 mm wave (A->P) in a 40 mm-sigma posterior-dorsal patch
% Sources: canary, f_fieldspread, f_fieldspread2 (a second F realisation, seed + 1; made whenever F is: the F null
% is the union of both realisations' atoms over every K, design §14 F union amendment, 2026-10-12).
% Y.canary_truth [rows x 4] the noise-free joint coefficients of each plant; Z.canary the plant table.
%
% Z: .name .workRate .bands (band, lo, hi, peak_hz, rate) .scales (depth, n_tiles, diameter_mm, lo_mm, hi_mm,
%    scale_mm) .variants .tiles {1} .A .R {per depth} .sources .series .transfer .fspread .seed .spanS .channels
%    .kind = "joint"
% Y: Y.(z_<source>_<band>_v1) single complex [rows x T], Y.(valid_<band>), Y.mask, Y.yper (the data's periodic
%    rows at WorkRate: rheome.modes.fit draws its N5 surrogates from it)
%
% See also: rheome.modes.fields, rheome.modes.operator, rheome.modes.scaletransfer, rheome.modes.discover,
%   rheome.modes.jointcheck, rheome.jointfilterbank
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Depths double = 0:5
        opts.Octaves double = [1 2; 2 4; 4 8; 8 16; 16 32]
        opts.WorkRate (1,1) double = 200
        opts.FitRange (1,2) double = [1 45]
        opts.GenDepth (1,1) double = 6
        opts.Seed (1,1) double = 1
        opts.Template char = ''
        opts.Nulls (1,1) logical = true
        opts.Canary (1,1) logical = false
        opts.CanaryHz (1,1) double = 10
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
    chn = string(st.chan.Name(S.iSel));  clear st
    CS = rheome.scale.cleanspan(F, fs);  F = F(:, CS.first:CS.last);  mk = CS.mask(CS.first:CS.last);
    r = max(1, round(fs / opts.WorkRate));  fw = fs / r;
    F = resample(F.', 1, r).';  N = size(F, 2);
    mask = any(reshape([mk(:)', false(1, r*N - numel(mk))], r, N), 1);  mask = mask(1:N);

    %% scales: one window of bank members per depth, pooled at that depth
    dep = sort(opts.Depths(:));  Tp = rheome.scale.template(max([dep; opts.GenDepth]), Template=opts.Template);
    D = arrayfun(@(d) i_diameter(S, Tp, d), dep);
    hi = [Inf; sqrt(D(1:end-1) .* D(2:end))];  lo = [hi(2:end); D(end)^2 / hi(end)];
    if numel(dep) == 1, hi = Inf; lo = D * 2^-0.25; end
    Op = arrayfun(@(k) rheome.modes.operator(S, Tp, CutoffMM=lo(k), MaxWavelengthMM=hi(k), Depth=dep(k)), (1:numel(dep))', 'uni', 0);
    Op = vertcat(Op{:});
    sc = arrayfun(@(o) exp(mean(log(o.keptWavelengthMM))), Op);
    Z = struct('name', name, 'kind', "joint", 'workRate', fw, 'seed', opts.Seed, 'spanS', N / fw, 'channels', chn);
    Z.options = struct('Depths', dep', 'Octaves', opts.Octaves, 'GenDepth', opts.GenDepth, 'WorkRate', opts.WorkRate, ...
                       'FitRange', opts.FitRange);                       % rheome.modes.fit recomputes with these
    Z.scales = table(dep, arrayfun(@(o) size(o.W, 1), Op), D, lo, hi, sc, arrayfun(@(o) string(mat2str(round(o.keptWavelengthMM))), Op), ...
                     'VariableNames', {'depth','n_tiles','diameter_mm','lo_mm','hi_mm','scale_mm','members_mm'});
    tl = table();  A = [];
    for k = 1:numel(Op)
        t = Op(k).tiles;  t.depth = repmat(dep(k), height(t), 1);  t.scale_mm = repmat(sc(k), height(t), 1);
        tl = [tl; t]; %#ok<AGROW>
        A = blkdiag(A, Op(k).A);
    end
    Wj = vertcat(Op.W);
    Z.tiles = {tl};  Z.A = logical(A);  Z.R = {Op.R};
    Z.variants = table(1, NaN, NaN, height(tl), "joint", 'VariableNames', {'variant','cutoff_mm','depth','n_tiles','kind'});
    Z.transfer = rheome.modes.scaletransfer(S, Z.scales);

    %% octaves: the jointfilterbank's time members (rheome.modes.octave)
    nb = size(opts.Octaves, 1);  bn = compose("o%g_%g", opts.Octaves(:, 1), opts.Octaves(:, 2));
    bands = table(bn, opts.Octaves(:, 1), opts.Octaves(:, 2), nan(nb, 1), nan(nb, 1), ...
                  'VariableNames', {'band','lo','hi','peak_hz','rate'});
    Y = struct('mask', mask);  rows = cell(0, 6);  spD = [];

    function add(src, Xs, which)
        if isstruct(Xs), sp = Xs; else, sp = rheome.modes.split(Wj * Xs, fw, mask, FitRange=opts.FitRange); end
        if which == "ap", Yt = sp.ap; elseif which == "rand", Yt = rheome.modes.phaserand(sp.per, ...
                RandStream('mt19937ar', 'Seed', opts.Seed + 5)); else, Yt = sp.per; end
        if src == "data", Y.yper = single(sp.per);  spD = sp; end           % N4, N5 reuse the split; fit's N5 surrogates
        for bb = 1:nb
            [z, rt, valid] = rheome.modes.octave(Yt, fw, opts.Octaves(bb, :), mask);
            vn = sprintf('z_%s_%s_v1', src, bn(bb));  Y.(vn) = z;  Y.(sprintf('valid_%s', bn(bb))) = valid;
            bands.rate(bb) = rt;
            rows(end+1, :) = {string(vn), string(src), bn(bb), 1, nnz(valid), rheome.modes.hash(z)}; %#ok<AGROW>
        end
    end

    srcs = strings(0);
    if ~opts.Canary, add("data", F, "per");  srcs = "data"; end
    if opts.Nulls || opts.Canary
        E = i_emptyroom(name, chn, fw, N);
        [Xf, Z.fspread] = i_fieldspread(S, Tp, F, fw, mask, E, opts, RandStream('mt19937ar', 'Seed', opts.Seed + 17));
        add("f_fieldspread", Xf, "per");  srcs(end+1) = "f_fieldspread";
        Xf2 = i_fieldspread(S, Tp, F, fw, mask, E, opts, RandStream('mt19937ar', 'Seed', opts.Seed + 18));
        add("f_fieldspread2", Xf2, "per");  srcs(end+1) = "f_fieldspread2";   % F union: both realisations are the null
    end
    if opts.Canary
        [Xc, Z.canary, B] = i_canary(S, Xf, fw, opts, RandStream('mt19937ar', 'Seed', opts.Seed + 23));
        Y.canary_truth = single(Wj * B);
        add("canary", Xc, "per");  srcs = ["canary" srcs];
    elseif opts.Nulls
        add("n4_aperiodic", spD, "ap");  add("n5_phaserand", spD, "rand");  srcs = [srcs "n4_aperiodic" "n5_phaserand"];
    end
    Z.bands = bands;  Z.sources = srcs;
    Z.series = cell2table(rows, 'VariableNames', {'var','source','band','variant','n_valid','sha256'});
end

function d = i_diameter(S, Tp, depth)
% the participant's median tile diameter at a depth, 2 sqrt(A / pi), mm
    a = [];
    for h = ["L" "R"]
        P = rheome.scale.grouptiles(S, Tp, h, depth);  m = full(sum(S.B.(h).lbo.Mass, 2));
        a = [a, full(sum(P .* m, 1))]; %#ok<AGROW>
    end
    d = 2e3 * sqrt(median(a) / pi);
end

function [X, Q] = i_fieldspread(S, Tp, F, fw, mask, E, o, rs)
% null F: independent sources, each generator tile with its own spectrum and power (NNLS on the data's tile
% spectra minus the empty room's), through the gain; plus the empty room on its channels
    Og = rheome.modes.operator(S, Tp, CutoffMM=1, Depth=o.GenDepth, GenCov=true);
    N = size(F, 2);  nf = 2^nextpow2(8 * fw);  ng = size(Og.Cgen, 3);
    [P, f] = rheome.scale.cleanwelch(Og.W * F, fw, 1:N, mask, nf);
    Pe = 0;  if ~isempty(E), Pe = rheome.scale.cleanwelch(Og.W(:, E.cols) * E.F, fw, 1:N, false(1, N), nf); end
    Tg = max(P - Pe, 0);
    Mx = zeros(size(Og.W, 1), ng);  for j = 1:ng, Mx(:, j) = sum((Og.W * Og.Cgen(:, :, j)) .* Og.W, 2); end
    eb = 0:0.5:fw/2;  bi = discretize(f, eb);  fb = eb(1:end-1)' + 0.25;  nbin = numel(fb);
    Sg = zeros(ng, nbin);  res = nan(nbin, 1);  ok = false(nbin, 1);
    for k = 1:nbin
        t = mean(Tg(bi == k, :), 1)';  if isempty(t) || ~all(isfinite(t)) || ~any(t), continue; end
        Mn = Mx / max(Mx(:));  tn = t / max(t);  s = lsqnonneg(Mn, tn);     % both scaled: lsqnonneg's tolerance is absolute
        Sg(:, k) = s * max(t) / max(Mx(:));  res(k) = norm(Mn * s - tn) / norm(tn);  ok(k) = true;
    end
    fk = (0:N-1) * fw / N;  fk = min(fk, fw - fk);
    X = zeros(size(F));
    for j = 1:ng
        [U, d] = eig(Og.Cgen(:, :, j), 'vector');  [d, o2] = sort(max(d, 0), 'descend');
        kk = find(cumsum(d) >= 0.999 * sum(d), 1);  if isempty(kk) || ~any(d), continue; end
        g = sqrt(interp1(fb(ok), Sg(j, ok), fk, 'linear', 0) * fw / 2);     % unit white noise: one-sided PSD 2/fw
        w = real(ifft(fft(randn(rs, kk, N), [], 2) .* g, [], 2));
        X = X + (U(:, o2(1:kk)) .* sqrt(d(1:kk))') * w;
    end
    Pf = rheome.scale.cleanwelch(Og.W * X, fw, 1:N, false(1, N), nf);
    in = f >= o.FitRange(1) & f <= o.FitRange(2);  c = sum(Tg(in, :), 'all') / sum(Pf(in, :), 'all');
    X = X * sqrt(c);
    if ~isempty(E), X(E.cols, :) = X(E.cols, :) + E.F; end
    ib = fb >= o.FitRange(1) & fb <= o.FitRange(2);
    Q = struct('gen_depth', o.GenDepth, 'n_generators', ng, 'nnls_rel_err_median', median(res(ib), 'omitnan'), ...
               'recalibration', c, 'emptyroom', ~isempty(E), 'tile_logpower_err_median', ...
               median(abs(log(c * Pf(in, :) ./ Tg(in, :))), 'all', 'omitnan'));
    if ~isempty(E), Q.emptyroomS = E.durationS;  Q.emptyroomChannels = numel(E.cols); end
end

function [X, T, B] = i_canary(S, Xf, fw, o, rs)
% J5: four patterns planted in turn on F, each at F's sensor power in the plant's octave
    V = double(S.Sf.Vertices);  nr = double(S.Sf.VertNormals);  nr = nr ./ vecnorm(nr, 2, 2);  G = S.G;
    Gn = nr(:,1)' .* G(:, 1:3:end) + nr(:,2)' .* G(:, 2:3:end) + nr(:,3)' .* G(:, 3:3:end);
    c0 = mean(V(V(:,1) <= prctile(V(:,1), 10), :), 1);  c1 = mean(V(V(:,1) >= prctile(V(:,1), 90), :), 1);
    c2 = mean(V(V(:,3) >= prctile(V(:,3), 80) & V(:,1) <= median(V(:,1)), :), 1);
    g = @(c, s) exp(-vecnorm(V - c, 2, 2).^2 / (2 * s^2));
    Q = [exp(1i*2*pi/0.200 * V(:,1)), g(c1, 0.03) + exp(-1i*pi/2) * g(c0, 0.03), cos(2*pi/0.200 * V(:,3)), ...
         exp(1i*2*pi/0.050 * V(:,1)) .* g(c2, 0.04)];
    names = ["wave200_ap"; "pair_fo"; "standing200"; "fine50"];
    B = Gn * Q;  B = B ./ vecnorm(B);  N = size(Xf, 2);  t = (0:N-1) / fw;
    oc = find(o.Octaves(:, 1) <= o.CanaryHz & o.Octaves(:, 2) > o.CanaryHz, 1);
    fk = (0:N-1) * fw / N;  fk = min(fk, fw - fk);  in = fk >= o.Octaves(oc, 1) & fk < o.Octaves(oc, 2);
    Pb = sum(abs(fft(Xf, [], 2)).^2 .* in, 'all') / N^2;                   % F's sensor power in the octave
    st = zeros(1, N);  i0 = 1;  last = 0;                                    % states: 0.2-0.4 s, no repeat
    while i0 <= N
        k = randi(rs, 4);  while k == last, k = randi(rs, 4); end
        n = round((0.2 + 0.2 * rand(rs)) * fw);  st(i0:min(N, i0+n-1)) = k;  i0 = i0 + n;  last = k;
    end
    ramp = ones(1, max(1, round(0.02 * fw))) / max(1, round(0.02 * fw));
    X = Xf;  ph = 2*pi*rand(rs, 1, 4);
    for k = 1:4
        a = conv(double(st == k), ramp, 'same');
        X = X + sqrt(2 * Pb) * real(B(:, k) * (a .* exp(1i * (2*pi*o.CanaryHz*t + ph(k)))));
    end
    T = table(names, repmat(o.CanaryHz, 4, 1), repmat(compose("o%g_%g", o.Octaves(oc, 1), o.Octaves(oc, 2)), 4, 1), ...
              arrayfun(@(k) mean(st == k), (1:4)'), 'VariableNames', {'plant','hz','band','occupancy'});
end

function E = i_emptyroom(name, chn, fs, N)
% the participant's own empty room (noise.mat of rheome.scale.importsubject) on the shared good channels
% (as rheome.modes.fields: mirror-padded to the span when shorter)
    E = [];  f = fullfile(rheome.load.root(), name, 'noise.mat');
    if ~isfile(f), return; end
    L = builtin('load', f, 'nrec');  nrec = L.nrec;
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;
    [~, cols, jn] = intersect(chn, nn(ok), 'stable');  jo = find(ok);  jn = jo(jn);
    if isempty(cols), return; end
    X = resample(double(nrec.F(jn, :)).', round(fs), round(nrec.sfreq)).';  n0 = size(X, 2);
    while size(X, 2) < N, X = [X, fliplr(X)]; end %#ok<AGROW>
    E = struct('F', X(:, 1:N), 'cols', cols, 'durationS', n0 / fs);
end

% Author: Diellor Basha, 2026
