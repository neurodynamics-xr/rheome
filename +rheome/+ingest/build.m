function file = build(src, cfg, opts)
% INGEST.BUILD  One recording -> one store of constant-Q tile summaries, in +data/<name>/.
%
%   file = rheome.ingest.build(name, cfg)
%   file = rheome.ingest.build(name, cfg, Channels="MEG", Store="default", Mask=[], ...
%                       ChannelChunk=32, Overwrite=false, Verbose=true)
%   file = rheome.ingest.build('/path/to/recording.mat', cfg, ...)      % store beside the file
%
% Opens the recording through pagedrecording, selects the channels, runs rheome.ingest.bank,
% rheome.ingest.grid, rheome.ingest.reduce (in channel chunks, to bound memory) and rheome.ingest.rollup, and
% writes a MAT v7.3 store WITHOUT compression whose per-level arrays are TOP-LEVEL
% variables -- matfile slices only those (paged-io notes). The file is keyed by its
% parameters the way the repo keys its caches, with the full config hash in meta:
%
%     +data/<name>/ingest__F0.25__V10__morse.mat            whole record, default store
%     +data/<name>/ingest__F0.25__V10__morse__native__paged.mat   Store="native", Paged="always"
%
% The path is chosen by cfg.Paged: "auto" pages when the whole-record CWT buffer would
% exceed MaxBytes (the 2400 Hz native store), "always"/"never" force it. A paged store
% carries "__paged" in its name only when forced, so the auto choice is recorded in
% meta.paged and meta.plan rather than the name.
%
% Variables:  meta, bands, frame, grid, and per level L (zero-padded NN):
%     LNN_n [K x 1] uint32       LNN_sumX, LNN_sumX2 [K x C] double
%     LNN_absMax, LNN_min, LNN_max [K x C] single
%     LNN_energy [K x C x B_L] double   LNN_envMax [K x C x B_L] single   LNN_nCoi [K x B_L] uint32
%     with sensor positions (Groups): LNN_g_* over the sensor tree's internal nodes, and
%     the wavelength axis LNN_s_energy / LNN_s_envMax [K x C x S] and LNN_gs_* for groups
%     (meta.sbands, meta.sframe; rheome.ingest.spatial)
%
% THE STORE IS THE DIAGONAL ON EVERY AXIS. Time (design §3.13): level L carries only the
% temporal bands whose natural level is <= L. Position (select design §2b-c): channel rows
% only from cfg.ChannelMinTile up (grid.channelLevel); finer tiles are described by the
% sensor tree; a wavelength band is stored on a node only if the node is at least one of
% its wavelengths across (grid.spaceSlots, packed [K x nSlots]; grid.chanSpace for
% sensors). In memory every unit has every band at every level; the diagonal is taken
% when writing. Below the diagonal a query goes to raw samples (phase 2).
%
% ORIGINAL NOTE -- level L carries only the bands whose natural
% level (first tile >= the band's support, bands.naturalLevel) is <= L, listed in
% grid.bandsAt{L+1}. A band's third index at level L is its position in that list, not
% its j. Below its natural level a band's tile is a smoothed reading of a wider
% neighbourhood, not a measurement of its own extent (the residual is negative there);
% queries finer than the natural level go to raw samples. In memory (rheome.ingest.reduce,
% rheome.ingest.rollup) every band exists at every level; the diagonal is taken when writing.
%
% ⚠ A store with this name and a DIFFERENT config hash is never overwritten silently:
% error ingest:build:hash unless Overwrite=true. The same hash returns the existing file
% without rebuilding (Overwrite=true rebuilds).
%
% ⚠ Whole-record: the CWT runs over the record per channel (0.45 GB single at 600 Hz,
% 156 filters). Paged: rheome.ingest.pages gives each band a level and rheome.ingest.reducepaged runs
% one sub-bank per job and page (design §3.12); the 2400 Hz native store takes this path.
%
% ⚠ Deterministic: same store, same cfg, same release give identical level variables.
% Only meta.created differs between two builds.
%
% INPUTS:
%   src   dataset name under +data, or a path to a recording store
%   cfg   rheome.ingest.config
%   Channels      "MEG" (MEG type AND ChannelFlag == 1, the repo's selection idiom),
%                 another type string, or explicit row indices into the store
%   Store         "default" (recording.mat) | "native" (recording_native.mat)
%   Mask          [nT x 1] logical, true = valid sample; [] = all valid
%   ChannelChunk  channels per reduce call (32)
%   Overwrite     rebuild even if a store exists (false)
%   Verbose       print progress (true)
% OUTPUT:
%   file  path of the store
%
% See also: rheome.ingest.config, rheome.ingest.reduce, rheome.ingest.rollup, rheome.ingest.totable, rheome.pagedrecording
%
% Author: Diellor Basha, 2026

    arguments
        src (1,:) char
        cfg (1,1) struct
        opts.Channels     = "MEG"
        opts.Store        (1,1) string {mustBeMember(opts.Store, ["default","native"])} = "default"
        opts.Mask         logical = logical([])
        opts.ChannelChunk (1,1) double {mustBeInteger, mustBePositive} = 32
        opts.Overwrite    (1,1) logical = false
        opts.Verbose      (1,1) logical = true
        opts.Groups       = "auto"     % "auto" (the dataset's helmet), false, or an [nCh x 3] positions matrix
    end

    t0 = tic;
    pr0 = rheome.pagedrecording(src, 'Store', char(opts.Store));       % header only
    nT = pr0.NumSamples;  fs = pr0.SamplingFrequency;

    % channel selection
    if isnumeric(opts.Channels)
        iSel = double(opts.Channels(:));
    else
        typ = char(string(opts.Channels));
        iSel = find(strcmp(pr0.ChannelType(:), typ) & pr0.ChannelFlag(:) == 1);
        if isempty(iSel)
            error('ingest:build:channels', 'No good channels of type ''%s'' in %s.', typ, pr0.File);
        end
    end
    C = numel(iSel);

    % where the store goes, and whether one is already there
    outDir = fileparts(pr0.File);
    tag = '';
    if strcmp(opts.Store, "native"), tag = [tag '__native']; end
    if strcmp(cfg.Paged, 'always') && ~strcmp(cfg.Bank, 'frame'), tag = [tag '__paged']; end
    if strcmp(cfg.Bank, 'frame'), kernel = 'frame'; else, kernel = cfg.Wavelet; end
    file = fullfile(outDir, sprintf('ingest__F%g__V%d__%s%s.mat', cfg.FrameFloor, cfg.VoicesPerOctave, kernel, tag));
    if exist(file, 'file') == 2
        old = load(file, 'meta');
        if isfield(old.meta, 'hash') && strcmp(old.meta.hash, cfg.hash) && ~opts.Overwrite
            if opts.Verbose, fprintf('rheome.ingest.build: %s exists with the same hash; not rebuilt.\n', file); end
            return
        elseif ~opts.Overwrite
            error('ingest:build:hash', ...
                  'A store exists at %s with hash %s; this config hashes %s. Pass Overwrite=true to replace it.', ...
                  file, old.meta.hash, cfg.hash);
        end
    end

    [fb, bands, frame] = rheome.ingest.bank(nT, fs, cfg);
    g = rheome.ingest.grid(nT, fs, cfg);
    bytesPer = 8 + 8 * strcmp(cfg.Precision, 'double');
    WT_FACTOR = 5;   % measured: wt's working set is ~5x its output buffer (peak RSS 4.2 GB for a 0.45 GB buffer)
    isFrame = isa(fb, 'rheome.timefilterbank');
    paged = ~isFrame && (strcmp(cfg.Paged, 'always') || ...
            (strcmp(cfg.Paged, 'auto') && WT_FACTOR * frame.nF * nT * bytesPer > cfg.MaxBytes));
    plan = table();
    if paged, plan = rheome.ingest.pages(g, bands, cfg); end
    if opts.Verbose
        fprintf('rheome.ingest.build: %s -- %d ch, %d samples @ %g Hz, %d filters in %d bands, %d levels, %s\n', ...
                src, C, nT, fs, frame.nF, height(bands), g.Lmax + 1, ...
                i_mode(paged));
        if paged
            fprintf('  paged: %d jobs; levels %s; halo up to %d samples; buffer up to %.2f GB per channel\n', ...
                    height(plan), mat2str(unique(plan.level)'), max(plan.haloSamples), max(plan.bufferBytes)/1e9);
        end
    end

    % level 0 in channel chunks, concatenated along the channel axis
    T0 = [];  PK = [];  CP = [];
    for a = 1:opts.ChannelChunk:C
        idx = iSel(a:min(a + opts.ChannelChunk - 1, C));
        pr = rheome.pagedrecording(src, 'Store', char(opts.Store), 'Channels', idx, 'Precision', 'double');
        X  = read(pr, 1, nT).';                                   % [nT x chunk], read ONCE
        if paged
            % pages slice the chunk in memory: reading through pagedrecording per page
            % pulls the hyperslab of every channel each time (measured 6 GB peak RSS on
            % the 600 Hz store against 4.2 GB for the whole-record path)
            readfcn = @(s1, s2) X(s1:s2, :);
            Tc = rheome.ingest.reducepaged(readfcn, numel(idx), bands, g, frame, plan, Mask=opts.Mask, ...
                                    Precision=cfg.Precision, MaxBytes=cfg.MaxBytes);
        else
            Tc = rheome.ingest.reduce(X, fb, bands, g, frame, Mask=opts.Mask, ...
                               Precision=cfg.Precision, MaxBytes=cfg.MaxBytes);
        end
        % ⭐ THE SPECTRAL PEAK OF EACH BAND, AT THE BAND'S OWN LEVEL. Not mergeable and not
        % needed to be: the diagonal gives every band one home level, so this is computed
        % once, here, on the samples already in memory. Measured: 0.2 s per channel for all
        % 15 bands of a 600 s record.
        if cfg.Peaks
            Pc = rheome.ingest.peaks(X(:, 1), bands, g);
            for cc = 2:size(X, 2)
                Q = rheome.ingest.peaks(X(:, cc), bands, g);
                for L = 1:numel(Pc)
                    Pc(L).freq = cat(3, Pc(L).freq, Q(L).freq);
                    Pc(L).amp  = cat(3, Pc(L).amp,  Q(L).amp);
                    Pc(L).power = cat(3, Pc(L).power, Q(L).power);
                end
            end
            if isempty(PK)
                PK = Pc;
            else
                for L = 1:numel(PK)
                    PK(L).freq = cat(3, PK(L).freq, Pc(L).freq);
                    PK(L).amp  = cat(3, PK(L).amp,  Pc(L).amp);
                    PK(L).power = cat(3, PK(L).power, Pc(L).power);
                end
            end
        end
        % ⭐ COUPLING, AND IT MERGES. Accumulated at the SLOW band's level, where every band
        % gets the same 22.6 of its own cycles; the sums roll up, so a coarser estimate is
        % exact and free (rheome.ingest.couple).
        % ⚠ COUPLING IS A FRAME-BANK FEATURE. It reads the bank's per-member bins and gains
        % to build each band's analytic signal on a shared grid; the Morse reference path has
        % no such structure, so it is skipped there rather than faked.
        if ~isequal(cfg.Couple, false) && isa(fb, 'rheome.timefilterbank')
            Cc1 = rheome.ingest.couple(X(:, 1), fb, bands, g, Pairs=cfg.Couple, PhaseBins=cfg.PhaseBins);
            for cc = 2:size(X, 2)
                Q = rheome.ingest.couple(X(:, cc), fb, bands, g, Pairs=cfg.Couple, PhaseBins=cfg.PhaseBins);
                for L = 1:numel(Cc1)
                    Cc1(L).vec = cat(3, Cc1(L).vec, Q(L).vec);
                    Cc1(L).amp = cat(3, Cc1(L).amp, Q(L).amp);
                end
            end
            if isempty(CP)
                CP = Cc1;
            else
                for L = 1:numel(CP)
                    CP(L).vec = cat(3, CP(L).vec, Cc1(L).vec);
                    CP(L).amp = cat(3, CP(L).amp, Cc1(L).amp);
                end
            end
        end
        clear X
        if isempty(T0)
            T0 = Tc;
        else
            for f = ["sumX","sumX2","absMax","min","max","energy","envMax", ...
                     "sumAbs","sumSqrt","sumX3","sumX4"]
                if isfield(Tc, f), T0.(f) = cat(2, T0.(f), Tc.(f)); end
            end
        end
        if opts.Verbose
            fprintf('  channels %d-%d of %d  (%.0f s)\n', a, a + numel(idx) - 1, C, toc(t0));
        end
    end
    T = rheome.ingest.rollup(T0, g);

    % the position pyramid: sensor groups by recursive spectral bisection of the helmet,
    % and the wavelength axis: a tight log-itersine frame on the same graph
    tree = table();  Tg = {};  chanRows = [];  sbands = table();  sframe = struct();  G = struct();
    if ~isequal(opts.Groups, false)
        [tree, chanRows, G] = i_tree(src, opts.Groups, pr0, iSel);
        if ~isempty(tree) && cfg.Space
            tS = tic;
            prAll = rheome.pagedrecording(src, 'Store', char(opts.Store), 'Channels', iSel, 'Precision', 'double');
            [S0, sbands, sframe] = rheome.ingest.spatial(@(a, b) read(prAll, a, b).', G, g, cfg, Mask=opts.Mask);
            T{1}.senergy = S0.senergy;  T{1}.senvMax = S0.senvMax;
            T = rheome.ingest.rollup(T{1}, g);                                % re-roll with the spatial fields
            if opts.Verbose, fprintf('  wavelength axis: %d members (A = %.3f, B = %.3f), %.0f s\n', height(sbands), sframe.A, sframe.B, toc(tS)); end
        end
        if ~isempty(tree)
            Tg = rheome.ingest.groups(T, tree, chanRows);
            if opts.Verbose, fprintf('  sensor tree: %d nodes, %d internal, depth %d\n', height(tree), nnz(~tree.is_leaf), max(tree.depth)); end
        end
    end

    % meta
    meta = struct();
    meta.source      = pr0.File;
    if rheome.load.has(src)
        meta.name = src;                                          % a dataset name is the recording id
    else
        [d0, stem] = fileparts(pr0.File);  [~, dn] = fileparts(d0);
        meta.name = [dn '_' stem];                                % folder_stem, never a full path
    end
    meta.fs          = fs;
    meta.nT          = nT;
    meta.duration    = nT / fs;
    meta.C           = C;
    meta.iChannel    = iSel(:)';
    meta.ChannelName = pr0.ChannelName(iSel);
    meta.ChannelType = pr0.ChannelType(iSel);
    meta.units       = 'T';
    meta.cfg         = cfg;
    meta.hash        = cfg.hash;
    meta.nValid      = double(sum(T0.n));
    meta.masked      = ~isempty(opts.Mask);
    meta.paged       = paged;
    meta.plan        = plan;
    meta.release     = version;
    meta.toolboxes   = i_toolboxes();
    meta.commit      = i_commit();
    meta.created     = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    meta.buildSeconds = [];                                       % filled below

    % the diagonal: level L carries the bands whose natural level is <= L (design §3.13)
    bandsAt = arrayfun(@(L) find(bands.naturalLevel <= L)', 0:g.Lmax, 'UniformOutput', false);
    g.bandsAt = bandsAt;
    % the HOME level of a band: the one level where a non-mergeable per-band measurement
    % (the spectral peak) is computed and stored.
    g.homeAt = arrayfun(@(L) find(bands.naturalLevel == L)', 0:g.Lmax, 'UniformOutput', false);
    meta.tree = tree;                                             % [] when no positions
    meta.groupChannelRows = chanRows;
    if ~isempty(Tg), meta.groupNodes = Tg{1}.nodes; else, meta.groupNodes = []; end
    meta.sbands = sbands;  meta.sframe = sframe;                  % empty when no wavelength axis
    % the position diagonal, part 1: channel rows only from ChannelMinTile up (when a tree exists)
    if ~isempty(Tg), Lchan = max(0, ceil(log2(cfg.ChannelMinTile / cfg.FrameFloor))); else, Lchan = 0; end
    Lchan = min(Lchan, g.Lmax);
    g.channelLevel = Lchan;
    % ⭐ THE ENVELOPE EXEMPTION. min and max merge by extremum, so the pair on a tile CONTAINS
    % every sample in it: stored at every level they are a min/max mipmap, and an exact
    % level-of-detail envelope of one channel's trace costs one array read at any zoom. They
    % are also the two cheapest columns in the store, so the diagonal buys nothing by
    % dropping them -- 10.4 MB for the whole pyramid of a 600 s CTF recording against 41 MB for the group
    % energies the diagonal was really for.
    if isfield(cfg, 'ChannelEnvelope'), Lenv = 0 * cfg.ChannelEnvelope + Lchan * ~cfg.ChannelEnvelope; else, Lenv = Lchan; end
    g.envelopeLevel = Lenv;
    % part 2: a wavelength band on a node only if the node is at least one of its wavelengths across
    spaceSlots = table();  chanSpace = [];
    if ~isempty(Tg) && isfield(Tg{1}, 'senergy')
        [spaceSlots, chanSpace] = i_spaceslots(tree, Tg{1}.nodes, sbands, cfg.SpaceDiagonal);
    elseif isfield(T{1}, 'senergy')
        chanSpace = 1:height(sbands);
    end
    g.spaceSlots = spaceSlots;  g.chanSpace = chanSpace;
    S = struct('meta', meta, 'bands', bands, 'frame', frame, 'grid', g);
    for L = 0:g.Lmax
        P = T{L+1};  bb = bandsAt{L+1};
        S.(sprintf('L%02d_n', L))    = P.n;
        S.(sprintf('L%02d_nCoi', L)) = P.nCoi(:, bb);
        if L >= Lenv && L < Lchan
            for f = ["min","max"]                                 % the envelope exemption
                S.(sprintf('L%02d_%s', L, f)) = P.(f);
            end
        end
        if L >= Lchan
            for f = ["sumX","sumX2","absMax","min","max","sumAbs","sumSqrt","sumX3","sumX4"]
                if isfield(P, f), S.(sprintf('L%02d_%s', L, f)) = P.(f); end
            end
            S.(sprintf('L%02d_energy', L)) = P.energy(:, :, bb);
            S.(sprintf('L%02d_envMax', L)) = P.envMax(:, :, bb);
            if isfield(P, 'senergy')
                S.(sprintf('L%02d_s_energy', L)) = P.senergy(:, :, chanSpace);
                S.(sprintf('L%02d_s_envMax', L)) = P.senvMax(:, :, chanSpace);
            end
        end
        if ~isempty(Tg)
            Q = Tg{L+1};
            for f = ["sumX","sumX2","absMax","min","max","sumAbs","sumSqrt","sumX3","sumX4"]
                if isfield(Q, f), S.(sprintf('L%02d_g_%s', L, f)) = Q.(f); end
            end
            S.(sprintf('L%02d_g_energy', L)) = Q.energy(:, :, bb);
            S.(sprintf('L%02d_g_envMax', L)) = Q.envMax(:, :, bb);
            if isfield(Q, 'senergy')
                % packed [K x nSlots]: slot -> (sband_id, group_id) in grid.spaceSlots
                K = size(Q.senergy, 1);  nSl = height(spaceSlots);
                E = zeros(K, nSl);  V = zeros(K, nSl, 'single');
                for sl = 1:nSl
                    E(:, sl) = Q.senergy(:, spaceSlots.col(sl), spaceSlots.sband_id(sl));
                    V(:, sl) = Q.senvMax(:, spaceSlots.col(sl), spaceSlots.sband_id(sl));
                end
                S.(sprintf('L%02d_gs_energy', L)) = E;
                S.(sprintf('L%02d_gs_envMax', L)) = V;
            end
        end
    end
    % the coupling sums, at each pair's slow level
    if ~isempty(CP)
        prs = cell(1, g.Lmax + 1);
        for i = 1:numel(CP)
            L = CP(i).level;
            S.(sprintf('L%02d_cpVec', L)) = single(permute(CP(i).vec, [1 3 2]));
            S.(sprintf('L%02d_cpAmp', L)) = single(permute(CP(i).amp, [1 3 2]));
            prs{L+1} = CP(i).pairs;
        end
        g.pairsAt = prs;
        S.grid = g;
    end

    % the peaks, on the diagonal: one array per level that is some band's home
    if ~isempty(PK)
        for i = 1:numel(PK)
            L = PK(i).level;
            S.(sprintf('L%02d_pkFreq', L))  = single(permute(PK(i).freq,  [1 3 2]));
            S.(sprintf('L%02d_pkAmp', L))   = single(permute(PK(i).amp,   [1 3 2]));
            S.(sprintf('L%02d_pkPower', L)) = single(permute(PK(i).power, [1 3 2]));
        end
    end
    S.meta.buildSeconds = toc(t0);
    save(file, '-struct', 'S', '-v7.3', '-nocompression');

    if opts.Verbose
        d = dir(file);
        fprintf('rheome.ingest.build: wrote %s (%.0f MB) in %.0f s\n', file, d.bytes / 1e6, S.meta.buildSeconds);
    end
end

function [slots, chanSpace] = i_spaceslots(tree, nodes, sbands, diagonal)
% which (spatial band, node) pairs are stored: node diameter >= the band's shortest
% wavelength, plus the root for every band; without calibrated wavelengths, all of them.
% chanSpace: the bands a single sensor (diameter 0) carries.
    nS = height(sbands);  nG = numel(nodes);
    wl = sbands.wavelength_lo;
    if ~diagonal || any(isnan(wl)), wl = zeros(nS, 1); end
    rows = {};
    for b = 1:nS
        for i = 1:nG
            n = nodes(i);
            if tree.parent_id(n) == 0 || tree.diameter(n) >= wl(b)
                rows(end+1, :) = {size(rows, 1) + 1, b, n, i}; %#ok<AGROW>
            end
        end
    end
    slots = cell2table(rows, 'VariableNames', {'slot','sband_id','group_id','col'});
    chanSpace = find(wl <= 0)';
end

function [tree, chanRows, G] = i_tree(src, groupsOpt, pr0, iSel)
% the helmet graph ON THE STORE'S CHANNELS, in the store's order, from the dataset's
% channel file or from given positions; the tree and the wavelength axis share it
    tree = table();  chanRows = [];  G = struct();
    if isnumeric(groupsOpt)
        P = groupsOpt;
        if size(P, 1) ~= numel(iSel) || size(P, 2) ~= 3
            error('ingest:build:groups', 'Groups positions must be [%d x 3].', numel(iSel));
        end
    elseif rheome.load.has(src)
        try
            arr0 = rheome.sensors.meg(src);
        catch
            return                                                  % no channel file: no groups
        end
        [tf, loc] = ismember(pr0.ChannelName(iSel), arr0.Labels);
        if ~all(tf), return; end
        P = arr0.Pos(loc, :);
    else
        return
    end
    % distances without pdist/squareform (Statistics Toolbox is not on every machine)
    G2 = P * P';  Dm = sqrt(max(diag(G2) + diag(G2)' - 2*G2, 0));  aperture = max(Dm(:));
    Dm(1:size(P,1)+1:end) = Inf;
    arr = struct('Name', 'store channels', 'Kind', 'meg', 'Pos', P, 'Dim', 2, 'Pitch', median(min(Dm, [], 2)), ...
                 'Aperture', aperture, 'Labels', {pr0.ChannelName(iSel)}, 'nCh', size(P, 1));
    chanRows = 1:numel(iSel);
    G = rheome.sensors.graph(arr, 'Faces', false);
    tree = rheome.sensors.tree(G);
end

function m = i_mode(paged)
    if paged, m = 'paged'; else, m = 'whole record'; end
end

function tb = i_toolboxes()
    v = ver;
    keep = ismember({v.Name}, {'Signal Processing Toolbox', 'Wavelet Toolbox'});
    tb = struct('Name', {v(keep).Name}, 'Version', {v(keep).Version});
end

function c = i_commit()
    c = '';
    here = fileparts(fileparts(mfilename('fullpath')));
    [st, out] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', here));
    if st == 0, c = strtrim(out); end
end
% Author: Diellor Basha, 2026
