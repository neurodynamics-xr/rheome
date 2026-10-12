function R = jointcheck(files, outDir, opts)
% MODES.JOINTCHECK  The J5 canary validation of joint discovery (protocol section 14), over the canary heads.
%
%   R = rheome.modes.jointcheck(files, outDir)     % files: rheome.modes.jointfields(..., Canary=true) .mat paths
%   R = rheome.modes.jointcheck(files, outDir, NSplits=100, Workers=16)
%
% Joint discovery (rheome.modes.discover with the §4.2 K rule, NMF and cPCA sensitivity) on the CANARY series:
% four planted patterns on null F in each head (rheome.modes.jointfields, Canary=true; no real-data series exists
% in those files, so no real-data atom can be inspected). Then the J5 rule:
%   recovered  for wave200_ap, pair_fo and standing200: the best |r| of any atom of the plant's octave to the
%              head's noise-free planted joint coefficients (Y.canary_truth) is >= RMin (0.7) in EVERY head
%   fine50     reported, not gated: its best |r| and its T-bar -- the energy-weighted scale transfer of the
%              planted coefficients over the depths (the head's own T(s), rheome.modes.scaletransfer) and of
%              its best-matching atom
%   F          no atom of null F passes Level-1 (i): F's own atoms (discovered at K* on f_fieldspread) are each
%              compared with the atoms of a second, independent F realisation (f_fieldspread2); an F atom would
%              pass (i) only with |r| < RMin to every one of them. Pass = none does
% pass = all three recovered AND no F atom passes. Not gated, to read a failure: joint_canary_byK.csv, the best |r|
% of each plant at EVERY K's reference clustering -- a K rule that merged two plants (K* too small: with three
% heads a split half is one or two heads) reads differently from an instrument that cannot separate them.
% Rehearsal (two OMEGA heads that are not the canary, 20 splits): K* = 2 merged the wave and the pair (wave |r|
% 0.52, 0.63); at K = 6 the wave, pair and standing plants read 0.97-1.00 and fine50 0.25-0.33 (T-bar 0.42-0.61).
%
% Writes, besides discover's files: joint_canary.csv (plant x head), joint_canary_f.csv (F atoms),
% joint_canary_transfer.csv (T(s) per head), joint_canary_byK.csv, joint_canary.json (the verdict and its parts).
%
% See also: rheome.modes.jointfields, rheome.modes.discover, rheome.modes.scaletransfer
%
% Author: Diellor Basha, 2026

    arguments
        files string
        outDir (1,:) char
        opts.RMin (1,1) double = 0.7
        opts.NSplits (1,1) double = 100
        opts.Ks double = 2:15
        opts.Workers (1,1) double = 0
        opts.Seed (1,1) double = 1
    end
    if isscalar(files) && ~endsWith(files, ".mat"), files = string(strip(splitlines(strtrim(fileread(files))))); files(files == "") = []; end
    m0 = matfile(files(1));  Z0 = m0.Z;  pb = unique(Z0.canary.band, 'stable');
    M = rheome.modes.discover(files, outDir, Source="canary", NullSources=["f_fieldspread" "f_fieldspread2"], ...
                              Bands=pb, NMF=true, NSplits=opts.NSplits, Ks=opts.Ks, Workers=opts.Workers, Seed=opts.Seed);
    rows = table();  trf = table();  heads = strings(0);  byK = table();
    for p = 1:numel(files)
        m = matfile(files(p));  Z = m.Z;  tr = m.canary_truth;  tl = Z.tiles{1};  heads(end+1) = string(Z.name); %#ok<AGROW>
        T = Z.transfer.table;  t = T.T(arrayfun(@(d) find(T.depth == d), tl.depth));
        trf = [trf; addvars(T, repmat(string(Z.name), height(T), 1), repmat(Z.transfer.edge_mm, height(T), 1), ...
               'Before', 1, 'NewVariableNames', {'participant','edge_mm'})]; %#ok<AGROW>
        for k = 1:height(Z.canary)
            C = M.bands([M.bands.band] == Z.canary.band(k)).C;  x = double(tr(:, k));
            rr = abs(C' * (x / norm(x)));  [r, a] = max(rr);
            Bb = M.bands([M.bands.band] == Z.canary.band(k));
            byK = [byK; table(repmat(string(Z.name), numel(Bb.refK), 1), repmat(Z.canary.plant(k), numel(Bb.refK), 1), Bb.refK(:), ...
                   cellfun(@(c) max(abs(c' * (x / norm(x)))), Bb.refC(:)), 'VariableNames', {'participant','plant','K','r_best_atom'})]; %#ok<AGROW>
            rows = [rows; table(string(Z.name), Z.canary.plant(k), Z.canary.band(k), r, a, i_tbar(x, t), i_tbar(C(:, a), t), ...
                    'VariableNames', {'participant','plant','band','r_best_atom','atom','t_bar_planted','t_bar_atom'})]; %#ok<AGROW>
        end
    end
    gated = ["wave200_ap" "pair_fo" "standing200"];
    rec = arrayfun(@(g) all(rows.r_best_atom(rows.plant == g) >= opts.RMin), gated);
    fr = table();
    for b = 1:numel(M.bands)
        nc = M.bands(b).nullC;
        if ~isfield(nc, 'f_fieldspread') || ~isfield(nc, 'f_fieldspread2'), continue; end
        r = max(abs(nc.f_fieldspread' * nc.f_fieldspread2), [], 2);
        fr = [fr; table(repmat(M.bands(b).band, numel(r), 1), (1:numel(r))', r, r < opts.RMin, ...
              'VariableNames', {'band','f_atom','r_max_to_f2','passes_level1_i'})]; %#ok<AGROW>
    end
    fok = ~isempty(fr) && ~any(fr.passes_level1_i);
    R = struct('pass', all(rec) && fok, 'recovered', table(gated', rec', 'VariableNames', {'plant','recovered'}), ...
               'f_atoms_passing_level1_i', sum(fr.passes_level1_i), 'f_checked', ~isempty(fr), ...
               'r_min', opts.RMin, 'heads', heads, 'sha256', M.sha256);
    R.fine50 = rows(rows.plant == "fine50", :);
    writetable(rows, fullfile(outDir, 'joint_canary.csv'));  writetable(fr, fullfile(outDir, 'joint_canary_f.csv'));
    writetable(trf, fullfile(outDir, 'joint_canary_transfer.csv'));  writetable(byK, fullfile(outDir, 'joint_canary_byK.csv'));
    fid = fopen(fullfile(outDir, 'joint_canary.json'), 'w');  fprintf(fid, '%s', jsonencode(R, PrettyPrint=true));  fclose(fid);
    fprintf('[modes.jointcheck] J5 %s: recovered %s, F atoms passing Level-1 (i) %d\n', ...
            string(R.pass), mat2str(rec), R.f_atoms_passing_level1_i);
end

function tb = i_tbar(x, t)
% energy-weighted scale transfer of a joint coefficient vector (rows = tile x scale, t = T of each row's depth)
    e = abs(x).^2;  tb = sum(e .* t) / sum(e);
end

% Author: Diellor Basha, 2026
