function [Z, Y, Op] = fields(name, S, opts)
% MODES.FIELDS  One participant's complex band tile series z(t) and the nulls N1, N4, N5 (stage 1 of the atlas).
%
%   [Z, Y] = rheome.modes.fields(name)
%   [Z, Y] = rheome.modes.fields(name, S, CutoffMM=100, Depth=3, Seed=1)
%   [Z, Y, Op] = ...                    % Op: the rheome.modes.operator of each variant (cell, primary first)
%
% Steps (protocol section 4.1): the recording (good MEG channels, its clean span: rheome.scale.cleanspan) is
% resampled to WorkRate and mapped to the group tiles by rheome.modes.operator (normal current, low-passed at
% a wavelength cutoff, tile mean). Each tile series is split into periodic and aperiodic parts
% (rheome.modes.split); the periodic part gives the analytic band signal (rheome.modes.analytic) per band.
% Band edges: the fixed band, re-centred on the participant's periodic peak where the median periodic PSD
% over primary tiles has an interior maximum above its aperiodic level (same width; one set of edges per
% participant, used for every variant and every null). Only bands no wider than MaxRecentreHz (5 Hz: delta,
% theta, alpha) are re-centred: the maximum of a 17 Hz beta band is not an edge, and re-centring it on a
% 20 Hz maximum moved its lower edge into alpha (11.6 Hz in one participant).
%
% Variants: the primary (CutoffMM, Depth) and the sensitivity cutoffs SensCutoffMM at Depth and depths
% SensDepths at CutoffMM. Sources (primary variant only, each through the SAME pipeline):
%   data          the recording
%   n1_sources    N1 leakage: independent normal sources at every vertex (variance ~ area), coloured by the
%                 recording's mean sensor spectrum, scaled to its sensor power, through the participant's own
%                 gain and the pipeline. Generated as chol(Gn diag(a) Gn')' * white, which has exactly that
%                 sensor covariance
%   n1_emptyroom  N1 sensor noise: the participant's own empty room (noise.mat) on the channels it shares,
%                 through the same W (absent channels zero); missing when the cache holds no noise run.
%                 A run shorter than the rest span is mirror-padded to its length (Z.emptyroomS = its own).
%                 A channel whose empty-room std exceeds 10x the rest's on that channel is dropped (a faulty
%                 sensor flagged good; MLC11-17 in 60 OMEGA-D empty rooms) and named in Z.emptyroomDropped
%   n4_aperiodic  N4: the aperiodic part of the data tile series (the complement of the periodic split)
%   n5_phaserand  N5: the data's periodic tile series with one random phase per frequency bin shared by all
%                 tiles (keeps every auto- and cross-spectrum; exact at tile level because everything before
%                 is linear and spatial and everything after is linear and time-invariant). One realisation
%                 here; rheome.modes.fit draws its own NSurr from Y.yper
%
% Z (metadata): .name .workRate .bands table (band, lo, hi, peak_hz, rate) .variants table (variant, cutoff_mm,
%   depth, n_tiles) .tiles {variant} tables .A .R (primary adjacency, tile resolution) .sources .series table
%   (var, source, band, variant, n_valid, sha256) .periodic table (tile x band periodic fraction, primary)
%   .seed .spanS .channels
% Y (arrays, each a top-level variable of the saved file so MatFile reads one without the rest):
%   Y.(z_<source>_<band>_v<k>) single complex [nTiles x T], Y.(valid_<band>) logical [1 x T],
%   Y.yper [nTiles x N] single (primary periodic series at WorkRate), Y.mask [1 x N]
%
% See also: rheome.modes.operator, rheome.modes.discover, rheome.modes.fit, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.CutoffMM (1,1) double = 100
        opts.Depth (1,1) double = 3
        opts.SensCutoffMM double = [92 130]
        opts.SensDepths double = [2 4]
        opts.BandNames string = ["delta" "theta" "alpha" "beta"]
        opts.Bands double = [1 4; 4 8; 8 13; 13 30]
        opts.WorkRate (1,1) double = 200
        opts.FitRange (1,2) double = [1 45]
        opts.Seed (1,1) double = 1
        opts.Template char = ''
        opts.Nulls (1,1) logical = true
        opts.MaxRecentreHz (1,1) double = 5
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
    chn = string(st.chan.Name(S.iSel));  clear st
    CS = rheome.scale.cleanspan(F, fs);  F = F(:, CS.first:CS.last);  mk = CS.mask(CS.first:CS.last);
    r = max(1, round(fs / opts.WorkRate));  fw = fs / r;
    F = resample(F.', 1, r).';  N = size(F, 2);
    mask = any(reshape([mk(:)', false(1, r*N - numel(mk))], r, N), 1);  mask = mask(1:N);

    var = [opts.CutoffMM opts.Depth; opts.SensCutoffMM(:) repmat(opts.Depth, numel(opts.SensCutoffMM), 1); ...
           repmat(opts.CutoffMM, numel(opts.SensDepths), 1) opts.SensDepths(:)];
    Tp = rheome.scale.template(max(var(:, 2)), Template=opts.Template);
    Op = cell(size(var, 1), 1);
    for v = 1:size(var, 1), Op{v} = rheome.modes.operator(S, Tp, CutoffMM=var(v, 1), Depth=var(v, 2)); end

    % primary data: split, band edges from the periodic peak
    Yd = Op{1}.W * F;  sp = rheome.modes.split(Yd, fw, mask, FitRange=opts.FitRange);
    bands = i_edges(sp, opts);
    Z = struct('name', name, 'workRate', fw, 'bands', bands, 'seed', opts.Seed, 'spanS', N / fw, 'channels', chn);
    Z.variants = table((1:size(var, 1))', var(:, 1), var(:, 2), cellfun(@(o) size(o.W, 1), Op), ...
                       'VariableNames', {'variant','cutoff_mm','depth','n_tiles'});
    Z.tiles = cellfun(@(o) o.tiles, Op, 'uni', 0);  Z.A = Op{1}.A;  Z.R = Op{1}.R;
    Z.periodic = i_pfrac(sp, bands, Op{1}.tiles);
    Y = struct('yper', single(sp.per), 'mask', mask);
    rows = cell(0, 6);

    function add(src, Yt, v, isPer)
        if ~isPer, s2 = rheome.modes.split(Yt, fw, mask, FitRange=opts.FitRange);  Yt = s2.per; end
        for b = 1:height(bands)
            [z, rt, valid] = rheome.modes.analytic(Yt, fw, [bands.lo(b) bands.hi(b)], mask);
            vn = sprintf('z_%s_%s_v%d', src, bands.band(b), v);  Y.(vn) = z;
            Y.(sprintf('valid_%s', bands.band(b))) = valid;  bands.rate(b) = rt;
            rows(end+1, :) = {string(vn), string(src), bands.band(b), v, nnz(valid), rheome.modes.hash(z)}; %#ok<AGROW>
        end
    end

    add("data", sp.per, 1, true);
    for v = 2:size(var, 1), add("data", Op{v}.W * F, v, false); end
    srcs = "data";
    if opts.Nulls
        rs = RandStream('mt19937ar', 'Seed', opts.Seed);
        add("n1_sources", Op{1}.W * i_sources(Op{1}.C1, F, fw, rs), 1, false);  srcs(end+1) = "n1_sources";
        E = i_emptyroom(name, chn, fw, N, F);
        if ~isempty(E)
            add("n1_emptyroom", Op{1}.W(:, E.cols) * E.F, 1, false);  srcs(end+1) = "n1_emptyroom";
            Z.emptyroomS = E.durationS;  Z.emptyroomChannels = numel(E.cols);  Z.emptyroomDropped = E.dropped;
        end
        add("n4_aperiodic", sp.ap, 1, true);  srcs(end+1) = "n4_aperiodic";
        add("n5_phaserand", rheome.modes.phaserand(sp.per, RandStream('mt19937ar', 'Seed', opts.Seed + 5)), 1, true);
        srcs(end+1) = "n5_phaserand";
    end
    Z.bands = bands;  Z.sources = srcs;
    Z.series = cell2table(rows, 'VariableNames', {'var','source','band','variant','n_valid','sha256'});
end

function B = i_edges(sp, o)
% fixed band, re-centred (same width) on an interior periodic peak above the aperiodic level
    Pper = median(max(sp.P - sp.Pap, 0), 2);  ratio = median(sp.P ./ sp.Pap, 2);  f = sp.f;
    n = numel(o.BandNames);  lo = o.Bands(:, 1);  hi = o.Bands(:, 2);  pk = nan(n, 1);
    for b = 1:n
        in = find(f > lo(b) & f < hi(b));
        if numel(in) < 3 || hi(b) - lo(b) > o.MaxRecentreHz, continue; end
        [~, j] = max(Pper(in));  j = in(j);
        if j > in(1) && j < in(end) && ratio(j) > 1 && Pper(j) > Pper(j-1) && Pper(j) > Pper(j+1), pk(b) = f(j); end
    end
    w = hi - lo;  L = lo;  H = hi;  ok = isfinite(pk);
    L(ok) = max(o.FitRange(1), pk(ok) - w(ok)/2);  H(ok) = min(o.FitRange(2), pk(ok) + w(ok)/2);
    B = table(o.BandNames(:), L, H, pk, nan(n, 1), 'VariableNames', {'band','lo','hi','peak_hz','rate'});
end

function T = i_pfrac(sp, B, tiles)
    T = table();
    for b = 1:height(B)
        m = sp.f >= B.lo(b) & sp.f < B.hi(b);
        pf = sum(max(sp.P(m, :) - sp.Pap(m, :), 0), 1) ./ sum(sp.P(m, :), 1);
        T = [T; table(repmat(B.band(b), height(tiles), 1), tiles.hemi, tiles.tile, pf(:), ...
             'VariableNames', {'band','hemi','tile','periodic_fraction'})]; %#ok<AGROW>
    end
end

function X = i_sources(C1, F, fs, rs)
% independent normal sources (variance ~ vertex area) through the gain, coloured by the mean sensor spectrum
    [U, D] = eig(C1, 'vector');  D = max(D, 0);
    N = size(F, 2);  X = (U .* sqrt(D)') * randn(rs, size(F, 1), N);
    nf = 2^nextpow2(8 * fs);  [P, f] = pwelch(F.', hann(nf), nf/2, nf, fs);  p = mean(P, 2);
    fk = (0:N-1)' * fs / N;  fk = min(fk, fs - fk);
    X = real(ifft(fft(X, [], 2) .* sqrt(interp1(f, p, fk, 'linear', 0)).', [], 2));
    X = X * sqrt(sum(F.^2, 'all') / sum(X.^2, 'all'));
end

function E = i_emptyroom(name, chn, fs, N, F)
% the participant's own empty room (noise.mat of rheome.scale.importsubject) on the shared good channels, without
% the channels whose level exceeds 10x the rest run's F on that channel (as rheome.modes.jointfields; an empty
% room is quieter than a head: ~0.5x)
    E = [];  f = fullfile(rheome.load.root(), name, 'noise.mat');
    if ~isfile(f), return; end
    L = builtin('load', f, 'nrec');  nrec = L.nrec;
    nn = string(nrec.ChannelName(:));  ok = nrec.ChannelFlag(:) == 1;
    [~, cols, jn] = intersect(chn, nn(ok), 'stable');  jo = find(ok);  jn = jo(jn);
    if isempty(cols), return; end
    X = resample(double(nrec.F(jn, :)).', round(fs), round(nrec.sfreq)).';  n0 = size(X, 2);
    bad = std(X, 0, 2) > 10 * std(F(cols, :), 0, 2);
    dropped = strjoin(chn(cols(bad)).', " ");  X = X(~bad, :);  cols = cols(~bad);
    if isempty(cols), return; end
    % shorter than the rest span: mirror-padded (no jump at the joins); N1 tests spatial structure, not dwell
    while size(X, 2) < N, X = [X, fliplr(X)]; end %#ok<AGROW>
    E = struct('F', X(:, 1:N), 'cols', cols, 'durationS', n0 / fs, 'dropped', dropped);
end

% Author: Diellor Basha, 2026
