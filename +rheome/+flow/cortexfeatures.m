function F = cortexfeatures(name, opts)
% FLOW.CORTEXFEATURES  The mergeable cortical features of a recording, written beside its store.
%
%   F = rheome.flow.cortexfeatures('sub01')                 % alpha (store band 6), depth-6 leaves
%   F = rheome.flow.cortexfeatures(name, BandId=6, Depth=6, Voices=2, Write=true)
%
% Per finest cell -- a depth-D cortical node (~47 mm at D=6) x a level-0 store tile (0.25 s) -- four
% SUMS, then everything coarser by rheome.select.rollupcortex (exact, both axes):
%   area_s       node area x tile length (m^2 s): the denominator
%   on_area_s    the same, counted only while the node is ON (rheome.detect.onstate on its tile-mean envelope)
%                -> occupancy = on_area_s / area_s at any node and level
%                ⚠ a leaf whose envelope is not bimodal (F.leafBimodal false) still gets labels, but
%                they are not a state; they are summed in all the same, as rheome.detect.onstate warns
%   energy       integral of envelope^2 over the node's area and the tile ((A m)^2 m^2 s)
%   energy_cb<m>  integral of (g_m * envelope)^2 per member m of a tight, log-spaced spatial bank
%                (logitersine, Voices per octave; sum_m g_m^2 = 1)
%
% ⭐ SUMS, SO A PARENT IS ITS CHILDREN. Occupancy is stored as its numerator and denominator because a
% mean does not add; everything here does, over nodes (area) and tiles (time).
% ⚠ the cband energies partition the envelope's energy only GLOBALLY: a filter spreads across node
% boundaries, so per node the octave energies need not sum to energy. And it is a FILTERED VIEW per
% octave, never a size (rheome.select.labelkinds).
% ⭐ ONE BANK FOR BOTH HEMISPHERES. The bank is designed on the union of the two hemispheres' spectral
% ranges and evaluated on each hemisphere's own eigenvalues, so cband m means the same octave on both
% sides (separate banks gave 11 and 12 members).
% ⚠ The envelope is averaged over the 0.25 s tile before squaring: these are energies of the level-0
% envelope, which is what removes the band's beating (the 8-16 Hz envelope varies up to ~4 Hz).
%
% Writes <store>__cortex.mat beside the ingest store (never the store itself): the cortex_node rows, and
% per level the four stats over every node. rheome.select.rows(db, "feature_cortex" | "feature_cortex_space" |
% "cortex_node", Level=L) reads it.
%
% See also: rheome.select.rollupcortex, rheome.geom.cortexnodes, rheome.flow.envelopemodes, rheome.detect.onstate, rheome.select.rows
%
% Author: Diellor Basha, 2026

    arguments
        name (1,1) string
        opts.BandId (1,1) double = 6
        opts.Depth (1,1) double {mustBeInteger, mustBePositive} = 6
        opts.Voices (1,1) double {mustBePositive} = 2
        opts.Write (1,1) logical = true
        opts.Store (1,1) string = ""
    end
    db = i_store(name, opts.Store);
    g = db.grid;  K = g.K(:)';  t0 = g.tExtent(1);
    E = rheome.flow.envelopemodes(name);  fr = E.L.fs;  bs = round(t0 * fr);  n0 = K(1);
    if size(E.L.C, 2) < n0 * bs
        error('flow:cortexfeatures:length', 'The envelope has %d frames; the grid needs %d.', size(E.L.C,2), n0*bs);
    end
    C = rheome.geom.cortexnodes(name, MaxDepth=opts.Depth);
    B = rheome.load.bases(name);
    nz = @(l) min(l(l > 1e-9 * max(l)));                          % first nonzero eigenvalue
    rng = [min(nz(B.L.lbo.Lambda), nz(B.R.lbo.Lambda)) max(max(B.L.lbo.Lambda), max(B.R.lbo.Lambda))];
    X = [];  ids = [];  sb = [];  bim = [];
    ws = warning;  warning('off', 'stats:gmdistribution:FailedToConverge');  warning('off', 'stats:gmdistribution:FailedToConvergeReps');  cw = onCleanup(@() warning(ws));   % 128 fits per run; .bimodal is the verdict
    for hh = ["L" "R"]
        H = B.(char(hh));  S = H.S;  lam = H.lbo.Lambda(:);  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;
        av = full(sum(Mm, 2));  T = C.T.(char(hh));  G = rheome.geom.tiles(T, S, opts.Depth);
        P = double(G.P);  aT = full(P' * av);  h = 2 + (hh == "R");
        gid = G.node_id + (h - 1) .* 2.^G.depth;                      % whole-cortex ids
        gfb = rheome.graphfilterbank(rng, 'Wavelet', 'logitersine', 'VoicesPerOctave', opts.Voices, ...
            'Transform', rheome.graphtransform.eigen(Phi, Mm, lam));
        Gm = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0));
        Cb = squeeze(mean(reshape(double(E.(char(hh)).C(:, 1:n0*bs)), size(Phi,2), bs, n0), 2));
        Yv = Phi * Cb;                                                % level-0 envelope per vertex
        yT = (P' * (av .* Yv)) ./ aT;                                 % tile-mean envelope per node
        on = false(numel(aT), n0);
        bh = false(numel(aT), 1);
        for k = 1:numel(aT), s = rheome.detect.onstate(max(yT(k,:), realmin)); on(k,:) = s.on; bh(k) = s.bimodal; end
        Xh = zeros(numel(aT), n0, 3 + gfb.NumMembers);
        Xh(:,:,1) = repmat(aT * t0, 1, n0);
        Xh(:,:,2) = double(on) .* (aT * t0);
        Xh(:,:,3) = (P' * (av .* Yv.^2)) * t0;
        for m = 1:gfb.NumMembers
            xm = Phi * (Gm(:, m) .* Cb);
            Xh(:,:,3+m) = (P' * (av .* xm.^2)) * t0;
        end
        X = [X; Xh]; ids = [ids; gid]; bim = [bim; bh]; %#ok<AGROW>
        sb = 1e3 * wavelengths(gfb);
    end
    R = rheome.select.rollupcortex(X, ids, K);
    F.recording_id = string(db.recording_id);  F.band_id = opts.BandId;  F.K = K;  F.tExtent = g.tExtent;
    F.ids = R.ids;  F.cbandWavelengthMM = sb(:);  F.depthLeaf = opts.Depth;
    F.leafIds = ids;  F.leafBimodal = bim;                            % which leaves' on/off is a state
    F.lev = R.lev;                                                    % [nNode x K x (3 + nSb)]
    F.stats = ["area_s" "on_area_s" "energy" compose("energy_cb%d", 1:numel(sb))];
    nodes = C.nodes(ismember(C.nodes.node_id, R.ids), :);
    nodes.recording_id = repmat(F.recording_id, height(nodes), 1);
    F.nodes = movevars(nodes, 'recording_id', 'Before', 1);
    F.file = rheome.select.cortexfile(db);
    if opts.Write, save(F.file, '-struct', 'F', '-v7.3'); end
end

function db = i_store(name, file)
    if file == ""
        Ct = rheome.select.catalog();
        hit = find(Ct.recording_id == name & Ct.bank == "frame" & Ct.store == "default", 1);
        if isempty(hit), error('flow:cortexfeatures:store', 'No default frame store for %s.', name); end
        file = Ct.file(hit);
    end
    db = rheome.select.open(char(file));
end

% Author: Diellor Basha, 2026
