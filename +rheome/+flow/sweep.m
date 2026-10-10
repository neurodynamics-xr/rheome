function [M, R, T] = sweep(name, varargin)
% FLOW.SWEEP  Sweep one LEVEL of the store's time tiling, measuring on every axis, and write the
% results into the database as measurements.
%
%   [M, R, T] = rheome.flow.sweep('sub01', Level=3)
%   M = rheome.flow.sweep(name, Level=3, BandIds=6, Kinds=["meanAct","curlRMS"], Write=false)
%
% One pass over the tiles at one level. Three families of measurement come back per tile:
%
%   TIME DOMAIN      rms, power, peak, ptp, envMax                 <- ALREADY IN THE STORE
%   FREQUENCY DOMAIN bandPower, share, spectralCentroid/Spread/Entropy/Flatness/Crest,
%                    dominantBand, dominantFreq, residual           <- ALREADY IN THE STORE
%   CORTICAL TILES   meanAct, energy, divRMS, curlRMS, normalShare, nVortex   <- computed here
%   EIGENMODE OCTAVES  spaceOctaveEnergy_sbJ, spaceOctaveShare_sbJ, GLOBAL (Scope="none")
%                      <- computed here. ⭐ A mode has no node: its support is the whole hemisphere,
%                      so per-tile eigenmode power needs a convention and has none. Summed over the
%                      hemisphere it is exact, so this is written with unit_id = 0.
%
% ⚠⚠ crest, kurtosis, skewness and the three factor statistics are DELIBERATELY ABSENT. They are
% only defined per channel, and channel rows start at db.channelLevel (level 6 = 16 s tiles here),
% so at the 2 s level this sweep targets they do not exist. Computing them from the group moments
% gives numbers -- crest 0.443-3.027 where absMax/rms is >= 1 by definition, kurtosis to 2.8e262 --
% and they are meaningless. See the note added to rheome.select.derive.
%
% ⭐⭐ THE FIRST TWO ARE NOT RE-MEASURED AND ARE NOT WRITTEN. They are exact functions of the
% mergeable moments the store already holds, so `rheome.select.derive` produces them at ANY level for the
% cost of an array read -- 34 statistics, listed by `rheome.select.derive()` with their formulas. Computing
% them again from samples would be slower and would give a second, drifting copy of a number the
% store already defines. This function only *joins* them.
%
% ⭐ ONLY THE CORTICAL MEASUREMENTS ARE NEW, and they go into the `measure` relation rather than a
% feature relation because they DO NOT MERGE: the curl over a patch is not the sum of the curls of
% its halves, and a time-averaged magnitude under a node is not a sum either. That is exactly the
% distinction `rheome.select.schema` draws between its two matrices -- 'mergeable' (dense, a parent is the
% merge of its children, prunable by descent) and 'note' (sparse, measured where it was measured,
% never descended). Putting a curl in a feature relation would make it look prunable.
%
% ⭐ THE WINDOWS ARE THE STORE'S OWN TILES. `rheome.select.derive` returns t_lo/t_hi per (level, k) and
% those spans drive the cortical pass, so every row joins on (level, k) with no interpolation and no
% second definition of "window". At level 3 on this store that is 300 tiles of 2 s.
% ⚠ Level 3 is also the measured sweet spot: the activation decorrelates in ~2 s and adjacent 2 s
% tiles are already 90% decorrelated, so rows are near-independent. Finer tiles are correlated;
% coarser ones average the variation away.
%
% ⚠ WHAT IS SAFE TO WRITE. Magnitude, divergence, curl and the gauge split under a cortical node
% are measurements. A SIZE is not: planting a source of known scale through this chain returns
% ~100 mm regardless of the truth below ~20 dB (measured by planting: slope +0.02, R^2 0.002 over a
% 70-267 mm range), while LOCATION recovers to 43-52 mm. So this writes where and how much, never
% how big, and `cortex_node.wavelength_mm` carries the admissibility of each node.
%
% ⚠ WRITES GO TO THE SIDECAR, NOT THE INGEST STORE. Labels and measurements share
% <store>__labels.mat; the ingest file is never modified. The write is one transaction, atomic and
% idempotent by content hash -- submitting the same rows twice returns the existing txn without
% writing -- and removable with rheome.select.rollback(db, txn_id). The txn_id is returned in M.
%
% NAME-VALUE
%   Level (3)        the time level to sweep
%   BandIds ([])     the store's band ids to measure; default = every band the level carries in the
%                    4-64 Hz range. ⚠ A level only carries bands whose natural level is at or below
%                    it, so asking for a low-frequency band at a fine level gets nothing.
%   Depth (3)        cortical tree depth; nodes at 1..Depth are measured
%   Kinds            which cortical measurements to write
%   GroupDepth (2)   sensor-tree depth for the derived sensor-side rows
%   Write (true)     write the cortical rows into the database
%   Store ("")       an explicit store file; default = the frame-bank store for this recording
%
% OUTPUTS
%   M  cortical measurements, long form: recording_id scope unit_id level k band_id kind value
%   R  the sensor-side derived statistics at the same level (time + frequency), from rheome.select.derive
%   T  the wide per-tile cortical table (rheome.flow.windowtable), if you want it before melting
%
% See also: rheome.select.measure, rheome.select.measures, rheome.select.derive, rheome.select.rollback, rheome.flow.windowtable,
%           rheome.select.schema
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Level', 3, @isscalar);
    p.addParameter('BandIds', [], @isnumeric);
    p.addParameter('Depth', 3, @isscalar);
    p.addParameter('GroupDepth', 2, @isscalar);
    p.addParameter('Kinds', ["meanAct","energy","crest","normalShare","divRMS","curlRMS","nVortex"]);
    p.addParameter('Write', true, @islogical);
    p.addParameter('Store', "", @(x) isstring(x) || ischar(x));
    p.addParameter('MaxWindows', Inf, @isscalar);
    p.addParameter('Verbose', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;  vb = o.Verbose;  L = o.Level;

    %% the store
    file = string(o.Store);
    if file == ""
        C = rheome.select.catalog();
        hit = find(C.recording_id == string(name) & C.bank == "frame" & C.store == "default", 1);
        if isempty(hit), error('flow:sweep:store', 'No default frame store for %s.', name); end
        file = C.file(hit);
    end
    db = rheome.select.open(char(file));
    if vb, fprintf('store %s\n  level %d: %d tiles\n', file, L, db.grid.K(L+1)); end

    %% bands this level carries, in the flow range
    bd = db.bands;
    cand = find(bd.fLo >= 4 & bd.fHi <= 64);
    if ~isempty(o.BandIds), cand = intersect(cand, o.BandIds(:)); end
    if isempty(cand), error('flow:sweep:bands', 'No bands in 4-64 Hz for this store.'); end
    bands = [bd.fLo(cand), bd.fHi(cand)];
    if vb, fprintf('  bands %s -> %s\n', mat2str(cand(:)'), ...
            strjoin(arrayfun(@(k) sprintf('%g-%g Hz', bands(k,1), bands(k,2)), 1:numel(cand), ...
            'UniformOutput', false), ', ')); end

    %% TIME + FREQUENCY: already in the store, so derive rather than measure
    % ⚠⚠ NOT Stats="all". In group scope every moment is a sum (or a maximum) over the node's
    % member channels, so a statistic that divides a MAXIMUM by a SUM is not interpretable --
    % measured crest 0.443-3.027 when absMax/rms is >= 1 by definition, and kurtosis to 2.8e262 by
    % catastrophic cancellation. Only the sums, the maxima and the band-energy ratios survive.
    % ⚠ The time-domain shape statistics (crest, kurtosis, skewness, shapeFactor, impulseFactor,
    % clearanceFactor) need Scope="channel", which needs Level >= db.channelLevel -- 6 here, a 16 s
    % tile. They are therefore NOT available at the 2 s level this sweep is built for, and
    % rheome.select.derive refuses rather than inventing them.
    GOOD = ["mean","power","rms","energy","envMax","min","max","ptp","bandPower","share", ...
            "residual","coneFrac","spectralEntropy","spectralFlatness","spectralCrest", ...
            "spectralCentroid","spectralSpread","dominantBand","dominantFreq"];
    R = rheome.select.derive(db, Level=L, Scope="group", Depth=o.GroupDepth, Stats=GOOD);
    tiles = unique(R(:, {'k','t_lo','t_hi'}), 'rows');
    tiles = sortrows(tiles, 'k');
    if isfinite(o.MaxWindows), tiles = tiles(1:min(height(tiles), o.MaxWindows), :); end
    if vb, fprintf('  derived %d sensor-side rows x %d statistics (nothing written)\n', ...
            height(R), width(R)-12); end

    %% CORTICAL: the only new measurement
    [T, X] = rheome.flow.windowtable(name, Tiles=[tiles.t_lo tiles.t_hi], Bands=bands, ...
        BandIds=cand(:)', MaxDepth=o.Depth, Flow=true, Verbose=vb);
    T.k = reshape(tiles.k(T.window), [], 1);
    T.level = repmat(L, height(T), 1);

    %% melt to the measurement matrix
    kinds = string(o.Kinds);
    kinds = kinds(ismember(kinds, string(T.Properties.VariableNames)));
    n = height(T) * numel(kinds);
    M = table('Size', [n 8], ...
        'VariableTypes', {'string','string','double','double','double','double','string','double'}, ...
        'VariableNames', {'recording_id','scope','unit_id','level','k','band_id','kind','value'});
    r = 0;
    for kk = kinds
        idx = r + (1:height(T));
        M.recording_id(idx) = string(name);  M.scope(idx) = "cortex";
        M.unit_id(idx) = T.node;  M.level(idx) = T.level;  M.k(idx) = T.k;
        M.band_id(idx) = T.band_id;  M.kind(idx) = kk;  M.value(idx) = T.(kk);
        r = r + height(T);
    end
    %% GLOBAL EIGENMODE BANDPOWER -- the only form in which it is well defined
    % ⭐⭐ AN EIGENMODE HAS NO NODE. Its support is the whole hemisphere, so "how much of mode m lives
    % in tile n" has no answer without a convention and different conventions disagree. Summed over
    % the hemisphere it is exact, so spatial-octave energy is written with Scope="none" and
    % unit_id = 0 -- a property of (time tile, temporal band, spatial octave) and of nothing smaller.
    % ⭐ The localisable twin already exists and is a DIFFERENT quantity: band-limit the current to
    % the octave (linear and exact) and integrate |J|^2 over the tile. That is a per-node number; it
    % is just not "the mode's power there".
    % ⚠ rheome.geom.ladder says the (tile x octave) plane collapses to at most ONE usable cell on this
    % cortex -- the 94-189 mm octave, which misses its own tile requirement by 0.1% -- so the
    % localised version is not worth an axis here. The global version is, and it is cheap.
    G = zeros(0, 7);
    if ~isempty(X.GroupPower)
        for bj = 1:numel(cand)
            P = X.GroupPower(:,:,bj);  P(~X.okGroup,:) = 0;
            tot = sum(P, 1);
            for j = 1:max(X.octave)
                sel = X.octave == j & X.okGroup;
                if ~any(sel), continue; end
                e = sum(P(sel,:), 1);
                for w = 1:numel(e)
                    G = [G; L, tiles.k(w), 0, cand(bj), j, e(w), e(w)/max(tot(w),realmin)]; %#ok<AGROW>
                end
            end
        end
    end
    if vb, fprintf('  cortical: %d rows = %d tiles x %d bands x %d nodes x %d kinds\n', ...
            height(M), height(tiles), numel(cand), numel(unique(T.node)), numel(kinds));
        fprintf('  global eigenmode octaves: %d rows (scope="none": a mode has no node)\n', size(G,1)*2);
    end
    for j = unique(G(:,5))'
        for nmv = ["spaceOctaveEnergy" "spaceOctaveShare"]
            col = 6 + (nmv == "spaceOctaveShare");
            sel = G(:,5) == j;
            add = table(repmat(string(name), sum(sel), 1), repmat("none", sum(sel), 1), ...
                G(sel,3), G(sel,1), G(sel,2), G(sel,4), ...
                repmat(sprintf('%s_sb%d', nmv, j), sum(sel), 1), G(sel,col), ...
                'VariableNames', M.Properties.VariableNames);
            M = [M; add]; %#ok<AGROW>
        end
    end

    %% write, one transaction per kind so a single measurement can be rolled back alone
    if o.Write
        % ⚠ ITERATE THE KINDS PRESENT IN M, NOT THE REQUESTED CORTICAL LIST. An earlier version
        % looped over o.Kinds, so the global eigenmode-octave rows were appended to the returned
        % table and silently never written -- 1800 rows that looked saved and were not. The scope
        % comes from the rows too, since the two families differ ("cortex" vs "none").
        allk = unique(M.kind);
        txn = zeros(numel(allk),1);  existed = false(numel(allk),1);
        for i = 1:numel(allk)
            m = M.kind == allk(i);
            sc = unique(M.scope(m));
            assert(isscalar(sc), 'kind %s spans two scopes', allk(i));
            rows = M(m, {'level','k','unit_id','band_id','value'});
            t = rheome.select.measure(db, rows, Kind=allk(i), Scope=sc, Source="rheome.flow.sweep");
            txn(i) = t.txn_id;  existed(i) = isfield(t,'existed') && t.existed;
        end
        kinds = allk;
        if vb
            fprintf('  wrote txn %s%s\n', mat2str(txn'), ...
                repmat(' (already present)', 1, any(existed)));
            fprintf('  roll back with: rheome.select.rollback(db, %d)\n', txn(1));
        end
        M.Properties.UserData = struct('txn_id', txn, 'store', file);
    end
end

% Author: Diellor Basha, 2026
