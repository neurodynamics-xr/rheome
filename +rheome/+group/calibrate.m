function [C, D] = calibrate(inDir, outDir, opts)
% GROUP.CALIBRATE  Does tile-wise group inference hold its error rate and find planted effects on the real cohort tile maps? Null labels and plants.
%
%   [C, D] = rheome.group.calibrate(inDir, outDir)
%   [C, D] = rheome.group.calibrate(inDir, outDir, NDraw=500, NPerm=2000, NSpin=1000, Effects=[0 0.2 0.35 0.5 0.8])
%
% inDir holds per-participant gaugetensor.csv (rheome.scale.run's <sub>/ or nsp cf-atlas's flat <sub>_ files,
% any depth) and metrics.csv, whose `dataset` column is the cohort; no draw mixes cohorts. The tile maps are
% the participants' own values per group tile (Measures: anisotropy, normal_share, log10 trace). No real outcome
% is read: the labels and covariates are random, so age and AD risk stay locked (AGE-1, ATLAS-1).
%
% Each draw, per cohort x measure, takes a random half (kind "group") or a standard-normal covariate (kind
% "covariate") and a random seed tile; the plant is w = exp(-theta^2 / (2 PlantDeg^2)) on the seed's hemisphere
% (theta: angle between tile centres on the template sphere), added as effect * sd(tile) * w to the labelled half
% (group: Cohen's d at the seed) or as effect * sd(tile) * w * z (covariate: a standardised slope). Effect 0 is
% the null-label condition. Each map goes through rheome.group.glm (Freedman-Lane, max-|t| FWER) and the t map
% through rheome.group.spin against w. Writes to outDir:
%   draws.csv         one row per cohort x measure x kind x effect x draw: any_fwer (any tile p_fwer < Alpha),
%                     fp_outside (any tile with w < 0.05 significant), detected (any tile with w >= 0.5),
%                     sensitivity (share of w >= 0.5 tiles), unc_rate (share of tiles with p < Alpha), spin_r, spin_p
%   calibration.csv   per cohort x measure x kind x effect: n_draws and the rates with Clopper-Pearson 95 % CIs;
%                     ⭐ calibrated = no false rate is above Alpha (each CI's lower bound <= Alpha): any_fwer and
%                     spin_p < Alpha at effect 0, fp_outside at every effect. Power (detected) is reported, not gated
%   tiles.csv         the tiles tested per cohort (hemi, depth, node_id, centre u1..u3)
%   participants.csv  the participants per cohort and measure (complete cases)
%
% Tiles missing in more than MaxMissing of a cohort's participants are dropped, then participants missing any
% kept tile (shortcut: complete cases; per-tile missingness would break the permutation's shared row order).
% Tile centres: Centres (a table hemi, node_id, u1, u2, u3) or rheome.scale.template(depth, Template=Template).
%
% See also: rheome.group.glm, rheome.group.spin, rheome.scale.reducegauge
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char
        opts.Measures string = ["anisotropy" "normal_share" "log10_trace"]
        opts.Kinds string = ["group" "covariate"]
        opts.Effects double = [0 0.2 0.35 0.5 0.8]
        opts.NDraw (1,1) double {mustBeInteger, mustBePositive} = 500
        opts.NPerm (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.NSpin (1,1) double {mustBeInteger, mustBePositive} = 1000
        opts.Alpha (1,1) double = 0.05
        opts.PlantDeg (1,1) double = 30
        opts.MaxMissing (1,1) double = 0.1
        opts.Seed (1,1) double = 7
        opts.Centres = []
        opts.Template char = ''
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    X = i_read(inDir);  dep = max(X.depth);  X = X(X.depth == dep, :);
    X.log10_trace = log10(X.trace);
    T = i_centres(opts, dep);
    rs = RandStream('mt19937ar', 'Seed', opts.Seed);
    D = {};  Tl = table();  Pp = table();
    for ds = unique(X.dataset)'
        Xd = X(X.dataset == ds, :);
        tk = unique(Xd(:, {'hemi','node_id'}));  [ok, ic] = ismember(tk, T(:, {'hemi','node_id'}));
        if ~all(ok), error('group:calibrate:centres', '%s: %d tiles have no centre', ds, nnz(~ok)); end
        U = T{ic, {'u1','u2','u3'}};
        Tl = [Tl; [table(repmat(ds, height(tk), 1), repmat(dep, height(tk), 1), 'VariableNames', {'dataset','depth'}) tk ...
                   array2table(U, 'VariableNames', {'u1','u2','u3'})]]; %#ok<AGROW>
        for ms = opts.Measures
            [Y, sub, keep] = i_matrix(Xd, tk, ms, opts.MaxMissing);
            n = size(Y, 1);  Uk = U(keep, :);  hk = tk.hemi(keep);
            Pp = [Pp; table(repmat(ds, n, 1), repmat(ms, n, 1), sub, 'VariableNames', {'dataset','measure','subject'})]; %#ok<AGROW>
            if n < 8 || size(Y, 2) < 2, warning('group:calibrate:small', '%s %s: %d participants, %d tiles: skipped', ds, ms, n, size(Y, 2)); continue, end
            sd = std(Y, 0, 1);  cosu = min(1, max(-1, Uk * Uk'));
            Sx = rheome.group.spin(zeros(size(Y, 2), 1), ones(size(Y, 2), 1), Uk, hk, NSpin=opts.NSpin, Seed=randi(rs, 2^31 - 1));
            fprintf('[group.calibrate] %s %s: %d participants x %d tiles\n', ds, ms, n, size(Y, 2));
            for d = 1:opts.NDraw
                g = zeros(n, 1);  g(randperm(rs, n, floor(n/2))) = 1;  z = randn(rs, n, 1);  z = (z - mean(z)) / std(z);
                s0 = randi(rs, size(Y, 2));  w = exp(-acosd(cosu(:, s0)).^2 / (2 * opts.PlantDeg^2)) .* (hk == hk(s0));
                sp = randi(rs, 2^31 - 1);
                for kind = opts.Kinds
                    if kind == "group", v = g; else, v = z; end
                    for e = opts.Effects
                        Ye = Y + e * (v .* (sd .* w'));
                        R = rheome.group.glm(Ye, [ones(n, 1) v], 2, NPerm=opts.NPerm, Seed=sp);
                        S = rheome.group.spin(R.t', w, Uk, hk, Index=Sx.index);
                        sig = R.p_fwer < opts.Alpha;  in = w >= 0.5;  out = w < 0.05;
                        D{end+1, 1} = table(ds, ms, kind, e, d, n, size(Y, 2), s0, any(sig), any(sig(out')), any(sig(in')), ...
                             mean(sig(in')), mean(R.p < opts.Alpha), S.r, S.p, 'VariableNames', {'dataset','measure','kind', ...
                             'effect','draw','n','n_tiles','seed_tile','any_fwer','fp_outside','detected','sensitivity', ...
                             'unc_rate','spin_r','spin_p'}); %#ok<AGROW>
                    end
                end
            end
        end
    end
    D = vertcat(D{:});  if isempty(D), D = table(); end
    C = i_summary(D, opts.Alpha);
    writetable(D, fullfile(outDir, 'draws.csv'));  writetable(C, fullfile(outDir, 'calibration.csv'));
    writetable(Tl, fullfile(outDir, 'tiles.csv'));  writetable(Pp, fullfile(outDir, 'participants.csv'));
    if ~isempty(C), fprintf('[group.calibrate] calibrated in %d of %d conditions\n', nnz(C.calibrated), height(C)); end
end

function C = i_summary(D, a)
    C = table();  if isempty(D), return, end
    [gi, C] = findgroups(D(:, {'dataset','measure','kind','effect'}));
    C.n_draws = splitapply(@numel, D.draw, gi);
    for f = ["any_fwer" "fp_outside" "detected" "spin"]
        if f == "spin", x = D.spin_p < a; else, x = D.(f); end
        k = splitapply(@sum, double(x), gi);
        [ph, ci] = binofit(k, C.n_draws);
        C.(f + "_rate") = ph;  C.(f + "_lo") = ci(:, 1);  C.(f + "_hi") = ci(:, 2);
    end
    C.sensitivity = splitapply(@mean, D.sensitivity, gi);  C.unc_rate = splitapply(@mean, D.unc_rate, gi);
    C.spin_r = splitapply(@mean, D.spin_r, gi);
    nul = C.effect == 0;
    C.calibrated = C.fp_outside_lo <= a & (~nul | (C.any_fwer_lo <= a & C.spin_lo <= a));
end

function X = i_read(inDir)
    f = dir(fullfile(inDir, '**', '*gaugetensor.csv'));  X = table();
    for k = 1:numel(f)
        if f(k).name == "gaugetensor.csv", [~, sub] = fileparts(f(k).folder); pre = ''; else
            sub = extractBefore(f(k).name, "_gaugetensor.csv"); pre = [sub '_']; end
        t = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string');
        t = t(:, {'hemi','depth','node_id','anisotropy','normal_share','trace'});  t.hemi = string(t.hemi);
        ds = "";  m = fullfile(f(k).folder, [pre 'metrics.csv']);
        if isfile(m)       % only the cohort is read from metrics: no outcome
            mm = readtable(m, 'TextType', 'string', 'Delimiter', ',');
            if ismember('dataset', mm.Properties.VariableNames) && height(mm), ds = string(mm.dataset(1)); end
        end
        t.subject = repmat(string(sub), height(t), 1);  t.dataset = repmat(ds, height(t), 1);
        X = [X; t]; %#ok<AGROW>
    end
    if isempty(X), error('group:calibrate:empty', 'no gaugetensor.csv under %s', inDir); end
    X.dataset(ismissing(X.dataset)) = "";
    [~, ia] = unique(X(:, {'subject','hemi','depth','node_id'}), 'last');  X = X(sort(ia), :);   % a re-run counts once
end

function [Y, sub, keep] = i_matrix(Xd, tk, ms, mx)
    sub = unique(Xd.subject);  [~, si] = ismember(Xd.subject, sub);  [~, ti] = ismember(Xd(:, {'hemi','node_id'}), tk);
    Y = nan(numel(sub), height(tk));  Y(sub2ind(size(Y), si, ti)) = Xd.(ms);
    Y(~isfinite(Y)) = NaN;
    keep = mean(isnan(Y), 1) <= mx;  Y = Y(:, keep);  r = all(isfinite(Y), 2);  Y = Y(r, :);  sub = sub(r);
end

function T = i_centres(opts, dep)
    if ~isempty(opts.Centres), T = opts.Centres;  T.hemi = string(T.hemi);  return, end
    Tp = rheome.scale.template(dep, Template=opts.Template);  T = table();
    for h = ["L" "R"]
        if ~isfield(Tp, h), continue, end
        th = Tp.(h);  c = floor(th.code / 2^(th.dmax - dep));  [ids, ~, tt] = unique(c);
        u = [accumarray(tt, th.u(:,1)) accumarray(tt, th.u(:,2)) accumarray(tt, th.u(:,3))];  u = u ./ vecnorm(u, 2, 2);
        T = [T; table(repmat(h, numel(ids), 1), ids(:), u(:,1), u(:,2), u(:,3), 'VariableNames', {'hemi','node_id','u1','u2','u3'})]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
