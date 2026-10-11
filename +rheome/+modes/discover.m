function M = discover(files, outDir, opts)
% MODES.DISCOVER  Group discovery: pooled phase-invariant k-means, K by split-half reproducibility, frozen + hashed.
%
%   M = rheome.modes.discover(files, outDir)           % files: fields .mat paths, or a text file of them
%   M = rheome.modes.discover(files, outDir, NPer=1000, Ks=2:15, NSplits=100, Workers=16)
%
% Protocol sections 4.2-4.4 on the discovery participants (whoever is listed: the fixed OMEGA-D half is the
% caller's subject list). Per band, on one variant (default the primary: >= 100 mm, depth 3):
%   snapshots   each participant gives NPer snapshots of its data series, equally spaced over the valid
%               samples whose norm is at or above the participant's MinAmpQuantile (the GFP-peak convention of
%               microstates: low-amplitude snapshots are mostly noise once normalised) -- equal samples per
%               participant, so no one dominates the pool
%   modes       rheome.modes.cluster on the pool: centroids c_k, |c_k| the amplitude topography, arg c_k the
%               phase layout
%   K           for every K in Ks: the reference clustering of the whole pool, then NSplits random participant
%               halves, each half clustered (NPerSplit snapshots per participant) and both Hungarian-matched
%               to the reference, so mode k is the same reference mode in every split; r_k = |<A_k, B_k>|.
%               K* = the largest K whose median r_k over splits is >= RThresh for EVERY k. No K passing is
%               recorded (k_reproducible = false) and K* falls back to min(Ks), flagged
%   sensitivity (a) exact DMD per participant (top DMDModes modes by amplitude), pooled and clustered at K*,
%               Hungarian |r| to each centroid; (b) complex PCA of the mean trace-normalised cross-spectral
%               matrix, best |r| of each centroid to its top K* components (a k-means mode equal to a cPCA
%               component carries no dynamics beyond second order). robust = r_dmd >= RobustR or
%               r_cpca >= RobustR
%   mixtures    mixture_share: the largest share of c_k's energy in the span of two other centroids. ⚠ The K rule
%               counts reproducible MIXTURES: with three planted states (mean dwell 300 ms) the alpha band-pass
%               smears every switch, the transition snapshots form their own reproducible clusters, and K* reached
%               6 with three modes >= 0.998 inside span{planted} (request 30825727 evidence). Reported, not acted on
%   null modes  the same clustering at K* on each null source's snapshots (N1 sources, N1 empty room, N4):
%               the centroid's best |r| to any of them (Level-1 (i) needs < 0.7 to every N1 mode)
% The centroids are frozen: centroids.mat (M) and centroids.sha256 = rheome.modes.hash over every band's C in
% band order, which rheome.modes.fit checks before it reads them.
%
% Writes to outDir: centroids.mat, centroids.sha256, modes_kcurve.csv, modes_centroids.csv,
% modes_discovery.csv, discover_sources.csv, discover_provenance.json.
%
% See also: rheome.modes.fields, rheome.modes.cluster, rheome.modes.match, rheome.modes.fit
%
% Author: Diellor Basha, 2026

    arguments
        files string
        outDir (1,:) char
        opts.Variant (1,1) double = 1
        opts.Bands string = string.empty
        opts.NPer (1,1) double = 1000
        opts.NPerSplit (1,1) double = 300
        opts.MinAmpQuantile (1,1) double = 0.5
        opts.Ks double = 2:15
        opts.NSplits (1,1) double = 100
        opts.RThresh (1,1) double = 0.8
        opts.RobustR (1,1) double = 0.7
        opts.Replicates (1,1) double = 5
        opts.SplitReplicates (1,1) double = 3
        opts.DMDModes (1,1) double = 10
        opts.NullSources string = ["n1_sources" "n1_emptyroom" "n4_aperiodic"]
        opts.Seed (1,1) double = 1
        opts.Workers (1,1) double = 0
    end
    if isscalar(files) && ~endsWith(files, ".mat"), files = strip(splitlines(strtrim(fileread(files)))); files(files == "") = []; end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    if opts.Workers > 0 && isempty(gcp('nocreate')), parpool('Processes', opts.Workers); end
    t00 = tic;

    %% 1. per participant: snapshots, cross-spectra, DMD modes (one file open at a time)
    P = numel(files);  m0 = matfile(files(1));  Z0 = m0.Z;  bands = Z0.bands.band;
    if ~isempty(opts.Bands), bands = intersect(bands, opts.Bands, 'stable'); end
    tiles = Z0.tiles{opts.Variant};  d = height(tiles);  nb = numel(bands);  src = ["data" opts.NullSources];
    X = cell(P, nb, numel(src));  CSD = zeros(d, d, nb);  DM = cell(P, nb);  info = table();
    for p = 1:P
        m = matfile(files(p));  Z = m.Z;  tp = Z.tiles{opts.Variant};
        if height(tp) ~= d || any(tp.tile ~= tiles.tile) || any(tp.hemi ~= tiles.hemi)
            error('modes:discover:tiles', '%s: its tiles differ from %s (variant %d)', files(p), files(1), opts.Variant);
        end
        for b = 1:nb
            valid = m.(sprintf('valid_%s', bands(b)));
            for s = 1:numel(src)
                if s > 1 && opts.Variant ~= 1, continue; end          % nulls exist on the primary variant only
                vn = sprintf('z_%s_%s_v%d', src(s), bands(b), opts.Variant);
                if ~ismember(vn, string(who(m))), continue; end
                z = double(m.(vn));  [X{p, b, s}, nAv] = i_snap(z, valid, opts);
                if s == 1
                    zv = z(:, valid);  C = zv * zv' / size(zv, 2);  CSD(:, :, b) = CSD(:, :, b) + C / real(trace(C)) / P;
                    DM{p, b} = i_dmd(z, valid, opts.DMDModes);
                    info = [info; table(string(Z.name), bands(b), nAv, size(X{p, b, s}, 2), 'VariableNames', ...
                            {'participant','band','n_available','n_snapshots'})]; %#ok<AGROW>
                end
            end
        end
        fprintf('[modes.discover] %d/%d %s\n', p, P, Z.name);
    end

    %% 2. per band: K curve, centroids, sensitivity, null modes
    kc = table();  cent = table();  disc = table();  Bd = struct([]);
    for b = 1:nb
        Xd = [X{:, b, 1}];
        Ks = opts.Ks(opts.Ks <= size(Xd, 2));  med = nan(numel(Ks), max(Ks));  ref = cell(numel(Ks), 1);
        for ik = 1:numel(Ks)
            K = Ks(ik);  ref{ik} = rheome.modes.cluster(Xd, K, Replicates=opts.Replicates, Seed=opts.Seed);
            r = nan(opts.NSplits, K);  sub = i_sub(X(:, b, 1), opts.NPerSplit);  C0 = ref{ik};
            sd = opts.Seed;  rp = arrayfun(@(s) randperm(RandStream('mt19937ar', 'Seed', sd + 1000*K + s), P), 1:opts.NSplits, 'uni', 0);
            parfor (s = 1:opts.NSplits, opts.Workers)
                h = rp{s};  A = [sub{h(1:floor(P/2))}];  B = [sub{h(floor(P/2)+1:end)}];
                CA = rheome.modes.cluster(A, K, Replicates=opts.SplitReplicates, Seed=sd + s);
                CB = rheome.modes.cluster(B, K, Replicates=opts.SplitReplicates, Seed=sd + s + 7);
                [~, pa] = rheome.modes.match(C0, CA);  [~, pb] = rheome.modes.match(C0, CB);
                r(s, :) = abs(sum(conj(CA(:, pa)) .* CB(:, pb), 1));
            end
            med(ik, 1:K) = median(r, 1);
            kc = [kc; table(repmat(bands(b), K, 1), repmat(K, K, 1), (1:K)', med(ik, 1:K)', repmat(min(med(ik, 1:K)), K, 1), ...
                  prctile(r, 25, 1)', 'VariableNames', {'band','K','mode','median_r','min_median_r','q25_r'})]; %#ok<AGROW>
            fprintf('[modes.discover] %s K=%2d  min median |r| %.3f  (%.0f s)\n', bands(b), K, min(med(ik, 1:K)), toc(t00));
        end
        ok = find(min(med, [], 2) >= opts.RThresh, 1, 'last');  repro = ~isempty(ok);
        if ~repro, ok = 1; end
        K = Ks(ok);  C = ref{ok};
        [~, lab] = max(abs(C' * (Xd ./ vecnorm(Xd))).^2, [], 1);  occ = accumarray(lab(:), 1, [K 1]) / numel(lab);
        % sensitivity: DMD pooled and clustered at K; cPCA of the mean cross-spectrum
        Dm = [DM{:, b}];  rd = nan(K, 1);
        if size(Dm, 2) >= K, Cd = rheome.modes.cluster(Dm, K, Replicates=opts.Replicates, Seed=opts.Seed);  rd = rheome.modes.match(C, Cd); end
        [V, L] = eig((CSD(:, :, b) + CSD(:, :, b)') / 2, 'vector');  [~, o] = sort(real(L), 'descend');  V = V(:, o(1:min(K, d)));
        rc = max(abs(C' * V), [], 2);
        mix = zeros(K, 1);                         % share of c_k in the span of two other modes (1 = a mixture)
        for k = 1:K
            o = setdiff(1:K, k);
            for i = o, for j = o(o > i), mix(k) = max(mix(k), norm(orth(C(:, [i j]))' * C(:, k))^2); end, end
        end
        nullr = struct();  nullC = struct();
        for s = 2:numel(src)
            Xn = [X{:, b, s}];  nullr.(src(s)) = nan(K, 1);
            if size(Xn, 2) < K, continue; end
            nullC.(src(s)) = rheome.modes.cluster(Xn, K, Replicates=opts.Replicates, Seed=opts.Seed);
            nullr.(src(s)) = max(abs(C' * nullC.(src(s))), [], 2);
        end
        Bd(b).band = bands(b);  Bd(b).K = K;  Bd(b).k_reproducible = repro;  Bd(b).C = C;  Bd(b).nullC = nullC;
        Bd(b).occupancy = occ;  Bd(b).r_dmd = rd;  Bd(b).r_cpca = rc;  Bd(b).cpca = V;  Bd(b).mixture_share = mix;
        T = table(repmat(bands(b), K, 1), (1:K)', repmat(K, K, 1), repmat(repro, K, 1), occ, rd, rc, ...
                  (rd >= opts.RobustR) | (rc >= opts.RobustR), mix, 'VariableNames', ...
                  {'band','mode','K','k_reproducible','occupancy_pooled','r_dmd','r_cpca','robust','mixture_share'});
        for s = 2:numel(src), T.("r_" + src(s) + "_max") = nullr.(src(s)); end
        disc = [disc; T]; %#ok<AGROW>
        for k = 1:K
            cent = [cent; table(repmat(bands(b), d, 1), repmat(k, d, 1), tiles.hemi, tiles.tile, abs(C(:, k)), angle(C(:, k)), ...
                    real(C(:, k)), imag(C(:, k)), 'VariableNames', {'band','mode','hemi','tile','amp','phase_rad','re','im'})]; %#ok<AGROW>
        end
        fprintf('[modes.discover] %s: K* = %d (%s)\n', bands(b), K, string(repro));
    end

    %% 3. freeze
    M = struct('bands', Bd, 'tiles', tiles, 'variant', opts.Variant, 'variantSpec', Z0.variants(opts.Variant, :), ...
               'participants', string(files(:)), 'options', opts, 'created', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
    M.sha256 = rheome.modes.hash(Bd.C);
    save(fullfile(outDir, 'centroids.mat'), 'M', '-v7');
    fid = fopen(fullfile(outDir, 'centroids.sha256'), 'w');  fprintf(fid, '%s  centroids.mat:M.bands.C\n', M.sha256);  fclose(fid);
    writetable(kc, fullfile(outDir, 'modes_kcurve.csv'));  writetable(cent, fullfile(outDir, 'modes_centroids.csv'));
    writetable(disc, fullfile(outDir, 'modes_discovery.csv'));  writetable(info, fullfile(outDir, 'discover_sources.csv'));
    pv = struct('participants', P, 'bands', bands, 'sha256', M.sha256, 'seconds', toc(t00), 'options', opts, ...
                'commit', i_git(), 'matlab', version);
    fid = fopen(fullfile(outDir, 'discover_provenance.json'), 'w');  fprintf(fid, '%s', jsonencode(pv, PrettyPrint=true));  fclose(fid);
end

function [S, nAv] = i_snap(z, valid, o)
% NPer snapshots, equally spaced over the valid samples at or above the participant's amplitude quantile
    a = vecnorm(z, 2, 1);  ok = valid & a >= quantile(a(valid), o.MinAmpQuantile);  idx = find(ok);  nAv = numel(idx);
    idx = idx(unique(round(linspace(1, nAv, min(o.NPer, nAv)))));
    S = z(:, idx) ./ a(idx);
end

function sub = i_sub(Xc, n)
    sub = cellfun(@(x) x(:, unique(round(linspace(1, size(x, 2), min(n, size(x, 2)))))), Xc, 'uni', 0);
end

function Phi = i_dmd(z, valid, nm)
% exact DMD on consecutive valid pairs; the nm modes of largest mean amplitude, unit norm
    t = find(valid(1:end-1) & valid(2:end));  X1 = z(:, t);  X2 = z(:, t + 1);
    [U, Sg, V] = svd(X1, 'econ');  r = nnz(diag(Sg) > max(size(X1)) * eps(Sg(1)));
    U = U(:, 1:r);  Sg = Sg(1:r, 1:r);  V = V(:, 1:r);
    [Wv, ~] = eig(U' * X2 * V / Sg);  Phi = X2 * V / Sg * Wv;  Phi = Phi ./ vecnorm(Phi);
    amp = mean(abs(Phi \ X1), 2);  [~, o] = sort(amp, 'descend');  Phi = Phi(:, o(1:min(nm, r)));
end

function s = i_git()
    [rc, s] = system(sprintf('git -C "%s" rev-parse HEAD', fileparts(fileparts(fileparts(mfilename('fullpath'))))));
    s = strtrim(s);  if rc ~= 0, s = ''; end
end

% Author: Diellor Basha, 2026
