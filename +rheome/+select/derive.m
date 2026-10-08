function R = derive(db, opts)
% SELECT.DERIVE  Per-tile derived statistics at one level, computed from the stored moments.
%
%   V = rheome.select.derive()                                    % the vocabulary, as a table
%   R = rheome.select.derive(db, Level=4)                         % every scalar statistic, all channels
%   R = rheome.select.derive(db, Level=4, Channels=7, Stats=["rms","spectralCentroid"])
%   R = rheome.select.derive(db, Level=6, Scope="group", Depth=2, Stats="all")
%   R = rheome.select.derive(db, Level=4, Window=[120 180], Stats=["share","bandPower"])
%
% One row per (unit, tile), in key order, with the requested statistics as columns. NOTHING
% here is stored: every column is a function of the tile's mergeable moments (n, sumX,
% sumX2, absMax, min, max) or of its band-energy vector (energy_j, envMax_j, nCoi_j). This
% is the "scalarization" layer -- what MATLAB's feature extractors compute per frame, except
% that here the frame can be any level of the pyramid and the moments were merged, not
% recomputed. Band-resolved columns (bandPower, share, bandEnergy, bandEnvMax) are
% [nRows x nBands] matrix variables; R.Properties.CustomProperties.bands names their order.
%
% ⚠ A NODE IS A SUM, NOT AN AVERAGE. In Scope="group" every moment is the exact sum (or
% maximum) over the node's member channels, so `energy` is the patch's total energy and
% `rms` is sqrt(n_sensors) times a typical channel's. The `n_sensors` column comes back with
% the rows for that reason; `share`, `spectralEntropy` and the other SPECTRAL shape statistics
% are ratios of band energies -- all sums over the same channels -- and are unaffected.
%
% ⚠⚠ BUT THE TIME-DOMAIN SHAPE STATISTICS ARE NOT, and being ratios does not save them: they
% divide a MAXIMUM over channels by a SUM over channels, which is not a ratio of like things.
% Measured on the reference frame store at level 3, Depth=2:
%     crest     0.443 - 3.027   <- absMax/rms is >= 1 BY DEFINITION; below 1 is proof it is broken
%     kurtosis  -1.1e264 - 2.8e262     skewness  -3.4e272
% The kurtosis and skewness blow-ups are catastrophic cancellation: sumX3 and sumX4 are summed
% over channels while the centre they are referred to is the sum of the channel means, which is
% not the patch's mean of anything.
% ⭐ VALID IN GROUP SCOPE: mean power rms energy envMax min max ptp bandEnergy bandEnvMax
%   bandPower share residual coneFrac spectralEntropy spectralFlatness spectralCrest
%   spectralCentroid spectralSpread dominantBand dominantFreq peakFreq peakAmp fftPower.
% ⚠ NOT VALID IN GROUP SCOPE: crest kurtosis skewness shapeFactor impulseFactor
%   clearanceFactor std var meanAbs -- use Scope="channel", which needs level >= channelLevel.
%   Asking for "all" in group scope returns them anyway; this is a caveat, not a guard, because
%   the right subset depends on what a patch is being taken to mean.
%
% ⚠ A PREVIEW STORE HOLDS ONLY min AND max (rheome.ingest.preview: no bank, no bands, no sums), so
% only min, max and ptp can be derived from it. "all" quietly narrows to those three there;
% naming any other statistic says what the store is instead of failing on a missing array.
%
% ⚠ BELOW THE CHANNEL FLOOR ONLY min, max AND ptp EXIST PER CHANNEL. Those three come from
% the min/max mipmap (rheome.ingest.config ChannelEnvelope), which reaches every level because the
% pair merges by extremum; asking for rms or a band statistic there says so and names the
% floor. Group scope has no floor at all.
%
% ⚠ A LEVEL CARRIES ONLY THE BANDS ON ITS SIDE OF THE DIAGONAL. The spectral statistics are
% computed over the bands the level actually carries, so spectralEntropy at level 2 and at
% level 8 are entropies of different-length vectors -- comparable down a column, not across
% levels. `nBands` is returned so a caller can see it, and `residual` says how much of the
% tile's power the carried bands miss.
%
% ⚠ RESIDUAL CAN BE NEGATIVE ON ONE TILE. The frame is tight over the WHOLE record, so the
% bands account for sumX2 exactly when summed over every tile and every band; within one
% tile a transient's wavelet response spreads across the tile edge, so a tile can hold more
% band energy than time-domain energy. Read it as "power not carried here", signed.
%
% ⚠ WINDOW KEEPS TILES THAT INTERSECT IT, not tiles inside it. A 64 s tile overlapping a 10 s
% window is the tile that covers that window; dropping it would leave the window empty.
%
% OPTIONS
%   Level     (required) pyramid level
%   Stats     names, or "all" (every scalar), "time", "spectral", "band", "fft"
%   Scope     "channel" (default) | "group" (sensor-tree nodes)
%   Channels  channel_ids ([] = all)          Groups/Depth  as in rheome.select.query
%   Bands     band_ids for the band-resolved and spectral columns ([] = all carried)
%   Window    [t0 t1] seconds ([] = whole record)
%
% See also: rheome.select.query, rheome.select.rows, rheome.select.frames, rheome.recordingbrowser
%
% Author: Diellor Basha, 2026

    arguments
        db = []
        opts.Level    double = []
        opts.Stats    string = "all"
        opts.Scope    (1,1) string {mustBeMember(opts.Scope, ["channel","group"])} = "channel"
        opts.Channels double = []
        opts.Groups   double = []
        opts.Depth    double = []
        opts.Bands    double = []
        opts.Window   double = []
    end

    V = i_vocab();
    if isempty(db), R = V; return; end
    if isempty(opts.Level)
        error('select:derive:level', 'Level is required.');
    end
    L = opts.Level;  g = db.grid;  scope = char(opts.Scope);
    if L < 0 || L > g.Lmax
        error('select:derive:level', 'Level %d is outside 0..%d.', L, g.Lmax);
    end
    want = i_expand(opts.Stats, V);
    if isfield(db, 'preview') && db.preview
        want = i_previewonly(want, V, opts.Stats);             % a preview holds min and max
    end
    need = unique([V.needs{ismember(V.name, want)}]);          % the moment arrays to read

    % ---- read only the moments the requested statistics need ----
    K = g.K(L+1);
    n = double(rheome.select.level(db, L, 'n'));  n = n(:);
    A = struct();  bandsAt = [];  units = [];
    for s = ["absMax","min","max","sumX2","mean","rms"]
        if ismember(s, need)
            [A.(char(s)), units] = sel_values(db, L, char(s), opts.Channels, [], scope, opts.Groups, opts.Depth);
        end
    end
    for s = ["energy","envMax","nCoi","peakFreq","peakAmp","fftPower"]
        if ismember(s, need)
            [A.(char(s)), units, bandsAt] = sel_values(db, L, char(s), opts.Channels, opts.Bands, scope, opts.Groups, opts.Depth);
        end
    end
    for s = ["sumAbs","sumSqrt","sumX3","sumX4"]
        if ismember(s, need)
            [A.(char(s)), units] = sel_values(db, L, char(s), opts.Channels, [], scope, opts.Groups, opts.Depth);
        end
    end
    if isempty(units)                                          % no moment needed a unit axis
        probe = 'absMax';
        if L < db.channelLevel && ~strcmp(scope, 'group'), probe = 'min'; end   % the envelope exemption
        [~, units] = sel_values(db, L, probe, opts.Channels, [], scope, opts.Groups, opts.Depth);
    end
    units = units(:);  U = numel(units);
    if isempty(bandsAt), bandsAt = g.bandsAt{L+1}; end
    nB = numel(bandsAt);
    fc = db.bands.fCenter(bandsAt);  fc = reshape(double(fc), 1, 1, nB);

    % ---- the window keeps every tile that INTERSECTS it ----
    ext = g.tExtent(L+1);  tc = g.tCenter{L+1}(:);
    lo = (0:K-1)' * ext;  hi = min(lo + ext, db.meta.duration);
    keep = true(K, 1);
    if ~isempty(opts.Window)
        w = sort(double(opts.Window(:)'));
        keep = hi > w(1) & lo < w(2);
        if ~any(keep)                                          % a window past the record: nearest tile
            [~, i] = min(abs(tc - mean(w)));  keep(i) = true;
        end
    end
    kk = find(keep);  nK = numel(kk);

    % ⚠ A NODE'S MOMENTS ARE SUMS OVER ITS SENSORS. g_sumX2 adds every member channel's
    % squares while n counts SAMPLES, so a node's rms is sqrt(nSensors) times a typical
    % channel's and its crest is that much smaller -- a node's "crest" is not a crest. The
    % sensor count travels with the rows so a caller can normalise rather than be misled.
    ns = ones(numel(units), 1);
    if strcmp(scope, 'group') && ~isempty(db.tree)
        ns = arrayfun(@(u) db.tree.n_sensors(db.tree.node_id == u), units(:));
    end

    % ---- the key columns ----
    [kg, ug] = ndgrid(kk, 1:U);
    R = table(repmat(string(db.recording_id), nK*U, 1), repmat(string(scope), nK*U, 1), ...
              units(ug(:)), repmat(L, nK*U, 1), kg(:), tc(kg(:)), repmat(ext, nK*U, 1), ...
              lo(kg(:)), hi(kg(:)), n(kg(:)), repmat(nB, nK*U, 1), ...
              'VariableNames', {'recording_id','scope','unit_id','level','k','t_center','t_extent','t_lo','t_hi','n','nBands'});
    R.n_sensors = ns(ug(:));

    nn = max(n(kk), 1);                                        % a fully masked tile has n = 0
    sub = @(M) reshape(M(kk, :), nK*U, 1);                     % [K x U] -> rows
    subB = @(M) reshape(M(kk, :, :), nK*U, nB);                % [K x U x nB] -> rows x bands

    % ---- band-energy quantities, shared by every spectral statistic ----
    if isfield(A, 'energy')
        E = double(A.energy);                                  % [K x U x nB]
        tot = sum(E, 3);
        sh  = E ./ max(tot, realmin);
    end

    for s = want
        switch s
            case 'mean',        R.mean = sub(A.mean);
            case 'rms',         R.rms  = sub(A.rms);
            case 'std'
                v = max(A.rms.^2 - A.mean.^2, 0) .* (n(:) ./ max(n(:)-1, 1));
                R.std = sqrt(sub(v));
            case 'var'
                R.var = sub(max(A.rms.^2 - A.mean.^2, 0) .* (n(:) ./ max(n(:)-1, 1)));
            case 'power',       R.power = sub(A.sumX2) ./ repmat(nn, U, 1);
            case 'peak',        R.peak = sub(A.absMax);
            case 'min',         R.min  = sub(A.min);
            case 'max',         R.max  = sub(A.max);
            case 'ptp',         R.ptp  = sub(A.max) - sub(A.min);
            case 'crest',       R.crest = sub(A.absMax) ./ max(sub(A.rms), realmin);
            case 'energy',      R.energy = sub(tot);
            case 'envMax',      R.envMax = sub(max(A.envMax, [], 3));
            case 'residual',    R.residual = 1 - sub(tot) ./ max(sub(A.sumX2), realmin);
            case 'coneFrac',    R.coneFrac = max(subB(double(A.nCoi)), [], 2) ./ repmat(nn, U, 1);
            case 'meanAbs',     R.meanAbs = sub(A.sumAbs) ./ repmat(nn, U, 1);
            case 'shapeFactor'
                R.shapeFactor = sub(A.rms) ./ max(sub(A.sumAbs) ./ repmat(nn, U, 1), realmin);
            case 'impulseFactor'
                R.impulseFactor = sub(A.absMax) ./ max(sub(A.sumAbs) ./ repmat(nn, U, 1), realmin);
            case 'clearanceFactor'
                R.clearanceFactor = sub(A.absMax) ./ max((sub(A.sumSqrt) ./ repmat(nn, U, 1)).^2, realmin);
            case 'skewness'
                [m3, sg] = i_central(A, n, kk, U, nn, sub);
                R.skewness = m3 ./ max(sg.^3, realmin);
            case 'kurtosis'
                [~, sg, m4] = i_central(A, n, kk, U, nn, sub);
                R.kurtosis = m4 ./ max(sg.^4, realmin);
            case 'peakFreq',   R.peakFreq  = subB(double(A.peakFreq));
            case 'peakAmp',    R.peakAmp   = subB(double(A.peakAmp));
            case 'fftPower',   R.fftPower  = subB(double(A.fftPower));
            case 'bandEnergy',  R.bandEnergy  = subB(E);
            case 'bandEnvMax',  R.bandEnvMax  = subB(double(A.envMax));
            case 'bandPower',   R.bandPower   = subB(E) ./ repmat(nn, U, nB);
            case 'share',       R.share       = subB(sh);
            case 'spectralEntropy'
                H = -sum(sh .* log2(max(sh, realmin)), 3);
                if nB < 2, H(:) = NaN; else, H = H / log2(nB); end
                R.spectralEntropy = sub(H);
            case 'spectralFlatness'
                R.spectralFlatness = sub(exp(mean(log(max(E, realmin)), 3)) ./ max(mean(E, 3), realmin));
            case 'spectralCrest'
                R.spectralCrest = sub(max(E, [], 3) ./ max(mean(E, 3), realmin));
            case 'spectralCentroid'
                R.spectralCentroid = 2 .^ sub(sum(sh .* log2(fc), 3));
            case 'spectralSpread'
                c = sum(sh .* log2(fc), 3);
                R.spectralSpread = sqrt(sub(sum(sh .* (log2(fc) - c).^2, 3)));
            case 'dominantBand'
                [~, i] = max(E, [], 3);
                R.dominantBand = reshape(bandsAt(i(kk, :)), nK*U, 1);
            case 'dominantFreq'
                [~, i] = max(E, [], 3);
                R.dominantFreq = reshape(double(db.bands.fCenter(bandsAt(i(kk, :)))), nK*U, 1);
        end
    end
    R = sortrows(R, {'unit_id','k'});
    R = addprop(R, {'bands','statNames'}, {'table','table'});
    R.Properties.CustomProperties.bands = bandsAt(:)';
    R.Properties.CustomProperties.statNames = want;
end

% A preview store can answer min, max and ptp. Asking for "all" there means those three;
% asking for rms by name is an error that names the store, not a missing variable.
function want = i_previewonly(want, V, asked)
    ok = ["min","max","ptp"];
    asked = string(asked(:)');
    named = setdiff(intersect(asked, V.name), ok);
    if ~isempty(named)
        error('select:derive:preview', ...
              ['This is a preview store: it holds min and max only, so %s cannot be derived ' ...
               'from it. Use min, max or ptp, or open the recording''s tile store.'], named(1));
    end
    want = intersect(want, ok, 'stable');
    if isempty(want), want = "ptp"; end
end

% The third and fourth CENTRAL moments from the raw sums, and the population sigma with them.
% ⚠ Central moments are differences of large like-sized terms; at these magnitudes (T^4 is
% ~1e-52) the cancellation is benign, but the order matters -- expand about the mean rather
% than subtracting at the end.
function [m3, sg, m4] = i_central(A, n, kk, U, nn, sub)
    mu = sub(A.mean);  r2 = sub(A.rms).^2;
    s3 = sub(A.sumX3) ./ repmat(nn, U, 1);
    v  = max(r2 - mu.^2, 0);
    sg = sqrt(v);
    m3 = s3 - 3*mu.*r2 + 2*mu.^3;
    m4 = [];
    if isfield(A, 'sumX4')
        s4 = sub(A.sumX4) ./ repmat(nn, U, 1);
        m4 = s4 - 4*mu.*s3 + 6*mu.^2.*r2 - 3*mu.^4;
    end
end

% ---- the vocabulary: name, kind, what it reads, and the formula, all in one place ----
% 'kind' drives the Stats shorthands and tells a caller whether a column is a vector per
% band; 'needs' is what sel_values must read, and nothing more.
function V = i_vocab()
    r = { 'mean',             'time',     "mean",                 'sumX / n'
          'rms',              'time',     "rms",                  'sqrt(sumX2 / n)'
          'std',              'time',     ["mean","rms"],         'sqrt(sumX2/n - mean^2) * sqrt(n/(n-1))'
          'var',              'time',     ["mean","rms"],         'std^2'
          'power',            'time',     "sumX2",                'sumX2 / n'
          'peak',             'time',     "absMax",               'absMax'
          'min',              'time',     "min",                  'min'
          'max',              'time',     "max",                  'max'
          'ptp',              'time',     ["min","max"],          'max - min'
          'crest',            'time',     ["absMax","rms"],       'absMax / rms'
          'energy',           'spectral', "energy",               'sum_j energy_j'
          'envMax',           'spectral', "envMax",               'max_j envMax_j'
          'residual',         'spectral', ["energy","sumX2"],     '1 - sum_j energy_j / sumX2'
          'coneFrac',         'spectral', "nCoi",                 'max_j nCoi_j / n'
          'spectralEntropy',  'spectral', "energy",               '-sum p log2 p / log2 nBands,  p = share'
          'spectralFlatness', 'spectral', "energy",               'geometric mean / arithmetic mean of energy_j'
          'spectralCrest',    'spectral', "energy",               'max energy_j / mean energy_j'
          'spectralCentroid', 'spectral', "energy",               '2^(sum p log2 fCenter_j)  (Hz, geometric)'
          'spectralSpread',   'spectral', "energy",               'sqrt(sum p (log2 fCenter_j - log2 centroid)^2)  (octaves)'
          'dominantBand',     'spectral', "energy",               'argmax_j energy_j'
          'dominantFreq',     'spectral', "energy",               'fCenter of the dominant band'
          'bandEnergy',       'band',     "energy",               'energy_j  (stored)'
          'bandEnvMax',       'band',     "envMax",               'envMax_j  (stored)'
          'bandPower',        'band',     "energy",               'energy_j / n'
          'share',            'band',     "energy",               'energy_j / sum_j energy_j'
          'meanAbs',          'time',     "sumAbs",               'sum|x| / n'
          'shapeFactor',      'time',     ["rms","sumAbs"],       'rms / meanAbs'
          'impulseFactor',    'time',     ["absMax","sumAbs"],    'absMax / meanAbs'
          'clearanceFactor',  'time',     ["absMax","sumSqrt"],   'absMax / (sum sqrt|x| / n)^2'
          'skewness',         'time',     ["mean","rms","sumX3"], 'third central moment / sigma^3'
          'kurtosis',         'time',     ["mean","rms","sumX3","sumX4"], 'fourth central moment / sigma^4'
          'peakFreq',         'fft',      "peakFreq",             'argmax of the tile''s spectrum in the band (Hz)'
          'peakAmp',          'fft',      "peakAmp",              'that peak''s amplitude, window-corrected'
          'fftPower',         'fft',      "fftPower",             'the band''s power from the same spectrum' };
    V = table(string(r(:,1)), string(r(:,2)), r(:,3), string(r(:,4)), ...
              'VariableNames', {'name','kind','needs','formula'});
end

% "all" is every SCALAR statistic: the band-resolved ones are matrix columns and a caller
% asks for them by name or by kind.
function want = i_expand(stats, V)
    stats = string(stats(:)');
    want = string.empty(1, 0);
    for s = stats
        switch s
            % ⚠ THE KIND NAMES MUST NOT COLLIDE WITH STAT NAMES. 'peak' is already a
            % statistic (absMax), so the per-tile FFT family is 'fft': asking for "peak"
            % must keep meaning the peak value, not expand to the spectral family.
            case "all",  want = [want, V.name(~ismember(V.kind, ["band","fft"]))'];  %#ok<AGROW>
            case {"time","spectral","band","fft"}, want = [want, V.name(V.kind == s)']; %#ok<AGROW>
            otherwise
                if ~ismember(s, V.name)
                    error('select:derive:stat', 'Unknown statistic "%s". rheome.select.derive() lists the vocabulary.', s);
                end
                want = [want, s];                                                   %#ok<AGROW>
        end
    end
    [~, i] = ismember(want, V.name);                 % vocabulary order, deduplicated
    want = V.name(unique(i))';
end
% Author: Diellor Basha, 2026
