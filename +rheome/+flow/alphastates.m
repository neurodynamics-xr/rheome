function Q = alphastates(name, opts)
% FLOW.ALPHASTATES  Query the alpha labels: states of one class in one hemisphere, with their episodes.
%
%   Q = rheome.flow.alphastates('sub01')                          % class 4 (widespread-sustained), R
%   Q = rheome.flow.alphastates(name, Class=[2 4], Hemi="L", Export=false)
%
% Reads what rheome.flow.labelalpha wrote (rheome.select.measures, Source "rheome.flow.labelalpha") and joins it on the
% store's own (level, k) -- no interpolation:
%   STATE     a run of consecutive level-3 tiles (2 s) in which one depth-6 leaf's alpha_class is in
%             Class; its span is the tiles' span, its onFrac the mean alpha_onFrac over them
%   EPISODE   an alpha_episode_* row of the same hemisphere (unit 2 / 3) whose level-0 tile (the one
%             holding tOn) falls inside a state's level-3 tiles: level-0 tile k0 is inside level-3
%             tile floor((k0-1)/8)+1
%
% ⚠ A STATE IS PER LEAF. Neighbouring leaves on at once give overlapping states; .union merges them
% into the hemisphere's time line of "some leaf in Class" runs, which is the one to count durations on.
% ⚠ The classes inherit rheome.select.labelkinds' calibration: "widespread" (3, 4) is PENDING a leakage
% calibration, so a class-4 state is "the ~130 mm ancestor was on", not "alpha covered 130 mm".
%
% OUTPUT (struct Q)
%   .states    table: state_id node_id hemi k_start k_end t_start t_end dur_s on_frac n_episodes
%   .episodes  table: state_id t_on t_off tau_r tau_f plateau_s r2 (one row per state x episode)
%   .union     table: k_start k_end t_start t_end dur_s n_leaves (hemisphere time line)
%   .files     exported CSVs under results/<subject>/labels/ (Export=true)
%
% See also: rheome.flow.labelalpha, rheome.select.measures, rheome.select.labelkinds, rheome.load.outpath
%
% Author: Diellor Basha, 2026

    arguments
        name (1,1) string
        opts.Class (1,:) double = 4
        opts.Hemi (1,1) string {mustBeMember(opts.Hemi, ["L","R"])} = "R"
        opts.Export (1,1) logical = true
        opts.Store (1,1) string = ""
    end
    db = i_store(name, opts.Store);  g = db.grid;  L3 = 3;  te = g.tExtent(L3+1);  bs = 2^L3;
    hNode = 2 + (opts.Hemi == "R");
    M = rheome.select.measures(db);  M = M(M.source == "rheome.flow.labelalpha" & M.scope == "cortex", :);
    if isempty(M), error('flow:alphastates:none', 'No rheome.flow.labelalpha rows for %s; run it first.', name); end
    isHemi = @(u) floor(u ./ 2.^(floor(log2(u)) - 1)) == hNode;       % whole-cortex heap: 2 L, 3 R
    C = M(M.kind == "alpha_class" & M.level == L3 & isHemi(M.unit_id) & ismember(M.value, opts.Class), :);
    O = M(M.kind == "alpha_onFrac" & M.level == L3, :);
    % states: runs of consecutive tiles per leaf
    st = zeros(0, 4);
    for u = unique(C.unit_id)'
        k = sort(C.k(C.unit_id == u));  br = [0; find(diff(k) > 1); numel(k)];
        for r = 1:numel(br) - 1
            ks = k(br(r)+1);  ke = k(br(r+1));
            of = O.value(O.unit_id == u & O.k >= ks & O.k <= ke);
            st(end+1, :) = [u ks ke mean(of)]; %#ok<AGROW>
        end
    end
    n = size(st, 1);
    S = table((1:n)', st(:,1), repmat(opts.Hemi, n, 1), st(:,2), st(:,3), (st(:,2)-1)*te, st(:,3)*te, ...
        (st(:,3)-st(:,2)+1)*te, st(:,4), zeros(n, 1), 'VariableNames', ...
        {'state_id','node_id','hemi','k_start','k_end','t_start','t_end','dur_s','on_frac','n_episodes'});
    % episodes of the hemisphere, each kind aligned by (k, write order)
    nm = ["tOn" "tOff" "tauR" "tauF" "plateauS" "r2"];
    Ek = M(M.kind == "alpha_episode_tOn" & M.unit_id == hNode, :);  Ek = sortrows(Ek, {'k','measure_id'});
    Ev = zeros(height(Ek), numel(nm));
    for j = 1:numel(nm)
        x = sortrows(M(M.kind == "alpha_episode_" + nm(j) & M.unit_id == hNode, :), {'k','measure_id'});
        Ev(:, j) = x.value;
    end
    k3 = floor((Ek.k - 1) / bs) + 1;
    ep = zeros(0, 7);
    for i = 1:n
        in = find(k3 >= S.k_start(i) & k3 <= S.k_end(i));
        S.n_episodes(i) = numel(in);
        ep = [ep; repmat(i, numel(in), 1) Ev(in, :)]; %#ok<AGROW>
    end
    Q.states = S;
    Q.episodes = array2table(ep, 'VariableNames', {'state_id','t_on','t_off','tau_r','tau_f','plateau_s','r2'});
    % hemisphere time line: tiles in which any leaf is in Class
    on = false(1, g.K(L3+1));  cnt = zeros(1, g.K(L3+1));
    for i = 1:n, on(S.k_start(i):S.k_end(i)) = true;  cnt(S.k_start(i):S.k_end(i)) = cnt(S.k_start(i):S.k_end(i)) + 1; end
    d = diff([0 on 0]);  ks = find(d == 1)';  ke = find(d == -1)' - 1;
    nl = arrayfun(@(a, b) max(cnt(a:b)), ks, ke);
    Q.union = table(ks, ke, (ks-1)*te, ke*te, (ke-ks+1)*te, nl, 'VariableNames', ...
        {'k_start','k_end','t_start','t_end','dur_s','n_leaves'});
    Q.files = strings(0, 1);
    if opts.Export
        tag = sprintf('alpha_states_%s_class%s', opts.Hemi, strjoin(string(opts.Class), ''));
        parts = ["states" "episodes" "union"];
        for p = parts
            f = string(rheome.load.outpath("labels", tag + "_" + p + ".csv", name));
            writetable(Q.(p), f);  Q.files(end+1, 1) = f;
        end
    end
end

function db = i_store(name, file)
    if file == ""
        Ct = rheome.select.catalog();
        hit = find(Ct.recording_id == name & Ct.bank == "frame" & Ct.store == "default", 1);
        if isempty(hit), error('flow:alphastates:store', 'No default frame store for %s.', name); end
        file = Ct.file(hit);
    end
    db = rheome.select.open(char(file));
end

% Author: Diellor Basha, 2026
