function L = labelalpha(name, opts)
% FLOW.LABELALPHA  Label a recording's alpha with its joint dyadic structure, into the select store.
%
%   L = rheome.flow.labelalpha('sub01')                              % all three families, written
%   L = rheome.flow.labelalpha(name, Families="class", Write=false)          % rows only
%
% Three note families, each written through rheome.select.measure (Scope="cortex", Source="rheome.flow.labelalpha")
% into the store's label sidecar, one transaction per kind -- idempotent by content, so a rerun returns
% the existing transactions (.existed), and removable with rheome.select.rollback:
%
%   CLASS       node = depth-6 leaf (~47 mm), tile = level 3 (2 s)
%     alpha_onFrac   the leaf's on-occupancy over the tile
%     alpha_class    modal class of the leaf's on-cells (level-0 tiles) inside the tile:
%                    1 focal-brief  2 focal-sustained  3 widespread-brief  4 widespread-sustained
%                    sustained  = the leaf's own occupancy >= Theta over an enclosing tile >= SustainS
%                    widespread = its depth-WideDepth ancestor (~130 mm at 3) has occupancy >= Theta
%                                 in the same level-0 tile
%   EPISODE     node = hemisphere (2 left, 3 right), tile = the level-0 tile holding tOn
%     alpha_episode_tOn _tOff _tauR _tauF _plateauS _r2
%                    rheome.detect.spindle (envelope mode) candidates on the hemisphere's PEAK amplitude
%                    (max over depth-7 tiles after a 0.1 s mean), located by rheome.detect.risefall
%   TRAJECTORY  node = start tile (depth 7, ~33 mm), tile = the level-3 window (2 s)
%     alpha_traj_netMM _speedMS _straight _durS _endNode _score
%                    rheome.detect.tilepath per window at Rate frames/s; baseline and threshold from the
%                    per-mode circular-shift surrogate at Alpha false paths per second
%                    (rheome.detect.tilepathnull); only windows whose best path beats the threshold
%
% ⭐ CLASS IS READ FROM THE MERGEABLE SIDECAR, not recomputed: occupancy at any node and level is
% on_area_s / area_s of rheome.flow.cortexfeatures' exact roll-up, so the labels and the features cannot
% disagree. The sidecar is built first if absent.
% ⚠ ONLY BIMODAL LEAVES ARE CLASSED. A leaf whose envelope has no two states (the sidecar's
% leafBimodal) has on/off labels that are not a state; it still contributes to its ancestors'
% occupancy (as in alpha_jointstate_omega.m) but gets no class or onFrac rows of its own.
% ⚠ "widespread" is registered PENDING: extent through the inverse is not leakage-calibrated, so a
% class-3/4 label says the ~130 mm node was on, not that the brain's alpha covered 130 mm.
% ⚠ Trajectory speed is COMPRESSED near the floor (0.02 m/s reads +65%) and a 40 mm amplitude swap
% reads 0.012-0.024 m/s -- see rheome.select.labelkinds before reading a speed as a velocity.
%
% OUTPUT (struct L)
%   .rows   struct of tables, one per kind (unit_id level k value)
%   .txn    struct of rheome.select.measure results, one per kind (empty when Write=false)
%   .null   per hemisphere: the trajectory null (rheome.detect.tilepathnull)
%
% See also: rheome.flow.cortexfeatures, rheome.select.measure, rheome.select.measures, rheome.select.rollback, rheome.select.labelkinds
%
% Author: Diellor Basha, 2026

    arguments
        name (1,1) string
        opts.Families (1,:) string {mustBeMember(opts.Families, ["class","episode","trajectory"])} = ["class","episode","trajectory"]
        opts.Write (1,1) logical = true
        opts.Theta (1,1) double = 0.75
        opts.SustainS (1,1) double = 8
        opts.WideDepth (1,1) double {mustBeInteger} = 3
        opts.TrajDepth (1,1) double {mustBeInteger} = 7
        opts.RateLevel (1,1) double {mustBeInteger} = 3               % 300 / 2^3 = 37.5 frames/s
        opts.Alpha (1,1) double = 0.1
        opts.MinR2 (1,1) double = 0.7
        opts.Seed (1,1) double = 17
        opts.Store (1,1) string = ""
    end
    db = i_store(name, opts.Store);  g = db.grid;  t0 = g.tExtent(1);
    L.rows = struct();  L.txn = struct();  L.null = struct();
    if any(opts.Families == "class"), L.rows = i_class(db, name, opts, L.rows); end
    if any(opts.Families ~= "class")
        E = rheome.flow.envelopemodes(name);  B = rheome.load.bases(name);  C = rheome.geom.cortexnodes(name, MaxDepth=opts.TrajDepth);
        fr = E.L.fs;  rng(opts.Seed);
        for hh = ["L" "R"]
            H = B.(char(hh));  S = H.S;  av = full(sum(H.lbo.Mass, 2));  h = 2 + (hh == "R");
            G = rheome.geom.tiles(C.T.(char(hh)), S, opts.TrajDepth, Ruler=rheome.geom.edgegraph(S));
            P = double(G.P);  Wd = (P' * (av .* H.lbo.Phi)) ./ (P' * av);
            gid = G.node_id + (h - 1) .* 2.^G.depth;
            Y = max(Wd * double(E.(char(hh)).C), 0);
            if any(opts.Families == "episode"), L.rows = i_episodes(L.rows, Y, fr, t0, g, h, opts); end
            if any(opts.Families == "trajectory")
                [L.rows, L.null.(char(hh))] = i_trajectories(L.rows, Y, double(E.(char(hh)).C), Wd, G, gid, fr, g, opts);
            end
        end
    end
    if opts.Write
        for kind = string(fieldnames(L.rows))'
            if height(L.rows.(kind)) == 0, continue; end
            L.txn.(kind) = rheome.select.measure(db, L.rows.(kind), Kind=kind, Scope="cortex", Source="rheome.flow.labelalpha");
        end
    end
end

% Class and onFrac per bimodal depth-6 leaf x level-3 tile, from the sidecar's occupancy pyramid.
function R = i_class(db, name, opts, R)
    f = rheome.select.cortexfile(db);
    if ~isfile(f), rheome.flow.cortexfeatures(name, Store=string(db.file)); end
    F = load(f, 'ids', 'lev', 'K', 'tExtent', 'leafIds', 'leafBimodal');
    occ = @(Lv) F.lev{Lv+1}(:,:,2) ./ F.lev{Lv+1}(:,:,1);
    [~, iLeaf] = ismember(F.leafIds, F.ids);
    on0 = occ(0);  on0 = on0(iLeaf, :) > 0.5;                         % binary at the leaf's own cell
    n0 = size(on0, 2);  k0 = 0:n0-1;
    sus = false(size(on0));
    for Lv = find(F.tExtent >= opts.SustainS) - 1
        O = occ(Lv);  kk = min(floor(k0 / 2^Lv) + 1, F.K(Lv+1));
        sus = sus | O(iLeaf, kk) >= opts.Theta;
    end
    dLeaf = floor(log2(F.leafIds(1))) - 1;                             % per-hemisphere depth
    [~, iAnc] = ismember(floor(F.leafIds / 2^(dLeaf - opts.WideDepth)), F.ids);
    O0 = occ(0);  wide = O0(iAnc, :) >= opts.Theta;
    cls = (1 + sus + 2 * wide) .* on0;                                 % 0 where off
    keep = find(F.leafBimodal);  L3 = 3;  bs = 2^L3;  K3 = F.K(L3+1);
    O3 = occ(L3);  u = [];  k = [];  vf = [];  vc = [];  uc = [];  kc = [];
    for i = keep(:)'
        c = reshape([cls(i, :) zeros(1, K3*bs - n0)], bs, K3);
        u = [u; repmat(F.leafIds(i), K3, 1)]; k = [k; (1:K3)']; vf = [vf; O3(iLeaf(i), :)']; %#ok<AGROW>
        c(c == 0) = NaN;  m = mode(c, 1);  has = isfinite(m);
        uc = [uc; repmat(F.leafIds(i), nnz(has), 1)]; kc = [kc; find(has)']; vc = [vc; m(has)']; %#ok<AGROW>
    end
    R.alpha_onFrac = table(u, repmat(L3, numel(u), 1), k, vf, 'VariableNames', {'unit_id','level','k','value'});
    R.alpha_class = table(uc, repmat(L3, numel(uc), 1), kc, vc, 'VariableNames', {'unit_id','level','k','value'});
end

% Episodes on the hemisphere's peak amplitude: spindle candidates, double-sigmoid fits.
function R = i_episodes(R, Y, fr, t0, g, h, opts)
    A = max(movmean(Y, max(1, round(0.1*fr)), 2), [], 1);
    Eb = rheome.detect.spindle(A, fr, struct('envelope', true, 'minDur', 0.4, 'maxDur', Inf, 'hiSD', 1.5, ...
        'loSD', 1.0, 'mergeGap', 0.1, 'smoothSec', 0.2, 'workRate', 100));
    v = zeros(0, 6);  kk = zeros(0, 1);
    for i = 1:height(Eb)
        Q = rheome.detect.risefall(A, fr, Eb.startSec(i), Eb.startSec(i) + Eb.durationSec(i));
        if ~Q.ok || Q.r2 < opts.MinR2 || Q.tOn < 0 || Q.tOn >= g.K(1) * t0, continue; end
        pl = 0;  if ~isempty(Q.plateau), pl = diff(Q.plateau); end
        v(end+1, :) = [Q.tOn Q.tOff Q.tauR Q.tauF pl Q.r2]; %#ok<AGROW>
        kk(end+1, 1) = floor(Q.tOn / t0) + 1; %#ok<AGROW>
    end
    nm = ["tOn" "tOff" "tauR" "tauF" "plateauS" "r2"];
    for j = 1:6
        R = i_append(R, "alpha_episode_" + nm(j), table(repmat(h, numel(kk), 1), zeros(numel(kk), 1), kk, v(:, j), ...
            'VariableNames', {'unit_id','level','k','value'}));
    end
end

% Trajectories per level-3 window: tilepath above the surrogate null.
function [R, N] = i_trajectories(R, Y, C, Wd, G, gid, fr, g, opts)
    L3 = 3;  nW = g.K(L3+1);  W = round(g.tExtent(L3+1) * fr);  bs = 2^opts.RateLevel;
    rate = fr / bs;  nF = floor(W / bs);  nW = min(nW, floor(size(Y, 2) / W));
    sh = randi([W, size(C, 2) - W], size(C, 1), 1);
    for m = 1:size(C, 1), C(m, :) = circshift(C(m, :), sh(m)); end
    Yr = i_winblock(Y, bs, W, nW);  Ys = i_winblock(max(Wd * C, 0), bs, W, nW);
    N = rheome.detect.tilepathnull(Ys, G, nF, rate, Alpha=opts.Alpha, WindowS=g.tExtent(L3+1));
    rows = zeros(0, 8);
    for w = 1:nW
        p = rheome.detect.tilepath(Yr(:, (w-1)*nF + (1:nF)), G, Baseline=N.baseline, SampleRate=rate);
        if p.nTracks == 0 || p.tracks(1).score <= N.threshold, continue; end
        p.tracks = p.tracks(1);  s = rheome.detect.trackstats(p, G, rate);  tl = p.tracks(1).tiles;
        rows(end+1, :) = [gid(tl(1)) w s.netMM s.speedMS s.straight s.durationS gid(tl(end)) p.tracks(1).score]; %#ok<AGROW>
    end
    nm = ["netMM" "speedMS" "straight" "durS" "endNode" "score"];
    for j = 1:6
        R = i_append(R, "alpha_traj_" + nm(j), table(rows(:, 1), repmat(L3, size(rows, 1), 1), rows(:, 2), rows(:, 2+j), ...
            'VariableNames', {'unit_id','level','k','value'}));
    end
end

function Y = i_winblock(Y0, bs, nW, nWin)
    nF = floor(nW / bs);  nT = size(Y0, 1);
    Z = reshape(Y0(:, 1:nW*nWin), nT, nW, nWin);  Z = Z(:, 1:nF*bs, :);
    Y = reshape(mean(reshape(Z, nT, bs, nF, nWin), 2), nT, nF*nWin);
end

function R = i_append(R, kind, T)
    T = T(~isnan(T.value), :);
    if isfield(R, kind), R.(kind) = [R.(kind); T]; else, R.(kind) = T; end
end

function db = i_store(name, file)
    if file == ""
        Ct = rheome.select.catalog();
        hit = find(Ct.recording_id == name & Ct.bank == "frame" & Ct.store == "default", 1);
        if isempty(hit), error('flow:labelalpha:store', 'No default frame store for %s.', name); end
        file = Ct.file(hit);
    end
    db = rheome.select.open(char(file));
end

% Author: Diellor Basha, 2026
