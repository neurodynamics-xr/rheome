function [G, P] = groupslowosc(inDir, outDir, opts)
% SCALE.GROUPSLOWOSC  The MS1 positive control across participants: one mean direction each, V-test for A->P, proportion and speed (G15, G16).
%
%   [G, P] = rheome.scale.groupslowosc(inDir, outDir)
%   [G, P] = rheome.scale.groupslowosc(inDir, outDir, SpeedRange=[1 7], Boot=2000, Seed=5)
%
% inDir holds one folder per participant (rheome.scale.run's OutDir) with slowosc.csv (rheome.scale.measure_slowosc),
% or nsp cf-slowosc's flat <sub>_slowosc.csv files (any depth).
% Per participant x estimator x param x condition (hemispheres pooled): the circular mean of the event
% directions (0 = anterior -> posterior), its resultant, and the median event speed -> participants.csv (P).
% Across participants (plan section 1.2 and G15) -> group.csv (G):
%   n             participants with at least MinEvents events
%   vtest_u/p     V-test of the participants' mean directions towards A->P (0), one-sided
%   frac_ap       proportion of participants whose mean direction is A->P (|angle| < 90 deg), with its 95%
%                 Clopper-Pearson CI (ap_lo, ap_hi)
%   speed_median  median of the participants' median speeds, 95% percentile bootstrap CI over participants
%                 (Boot resamples): speed_lo, speed_hi; and frac_speed_in_range, the proportion of participants
%                 whose median speed lies in SpeedRange (1-7 m/s, Massimini et al. 2004, as the plan quotes it)
% ⭐ The control passes for an estimator if, on the ORIGINAL events, the V-test rejects at alpha (Holm across
%   estimators x params is the reader's: the p values are reported raw) AND on the REVERSED events it does not
%   (the direction must flip), AND on the SURROGATES it does not. Nothing is pooled across cohorts here
%   (AnphySleep is the only cohort of the primary control).
% ⭐ THE LAGGED-SOURCE NULL (gate N) is not part of passes: two lagged generators may well read A->P. On the
%   original rows, lagged_excess is the median over participants of (A->P fraction of the events) minus (that of
%   their lagged nulls), lagged_p its one-sided Wilcoxon signed-rank p (real > lagged), lagged_n the pairs, and
%   beats_lagged = lagged_p < Alpha. A direction is a travel claim only where beats_lagged holds.
%
% See also: rheome.scale.measure_slowosc, rheome.scale.vtest, rheome.scale.reduce
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char = ''
        opts.SpeedRange (1,2) double = [1 7]
        opts.Boot (1,1) double = 2000
        opts.Seed (1,1) double = 5
        opts.MinEvents (1,1) double = 10
        opts.Alpha (1,1) double = 0.05
    end
    f = dir(fullfile(inDir, '**', '*slowosc.csv'));          % <sub>/slowosc.csv (run) or <sub>_slowosc.csv (nsp)
    P = table();
    for k = 1:numel(f)
        X = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string');
        if isempty(X), continue, end
        if f(k).name == "slowosc.csv", [~, sub] = fileparts(f(k).folder); else, sub = extractBefore(f(k).name, "_slowosc.csv"); end
        [gk, key] = findgroups(X(:, {'estimator','param','condition'}));
        th = deg2rad(X.angleDeg);
        z = splitapply(@(a) mean(exp(1i * a(isfinite(a)))), th, gk);
        n = splitapply(@(a) nnz(isfinite(a)), th, gk);
        sp = splitapply(@(s) median(s, 'omitnan'), X.speedMS, gk);
        key.fracAP = splitapply(@(a) mean(cos(a(isfinite(a))) > 0), th, gk);
        key.subject = repmat(string(sub), height(key), 1);
        key.nEvents = n;  key.meanAngleDeg = rad2deg(angle(z));  key.resultant = abs(z);  key.medianSpeedMS = sp;
        P = [P; key]; %#ok<AGROW>
    end
    G = table();
    if isempty(P), i_write(outDir, G, P); return, end
    P = P(P.nEvents >= opts.MinEvents, :);
    [gk, G] = findgroups(P(:, {'estimator','param','condition'}));
    rng(opts.Seed);
    for i = 1:height(G)
        m = gk == i;  th = deg2rad(P.meanAngleDeg(m));  sp = P.medianSpeedMS(m);  n = numel(th);
        [G.vtest_u(i), G.vtest_p(i)] = rheome.scale.vtest(th);
        x = nnz(cos(th) > 0);  G.n(i) = n;  G.frac_ap(i) = x / n;
        [G.ap_lo(i), G.ap_hi(i)] = i_clopper(x, n);
        z = mean(exp(1i * th));  G.mean_angle_deg(i) = rad2deg(angle(z));  G.resultant(i) = abs(z);
        sp = sp(isfinite(sp));  G.speed_median(i) = median(sp);
        if numel(sp) > 1
            bm = median(sp(randi(numel(sp), numel(sp), opts.Boot)), 1);  ci = prctile(bm, [2.5 97.5]);
        else, ci = [NaN NaN]; end
        G.speed_lo(i) = ci(1);  G.speed_hi(i) = ci(2);
        G.frac_speed_in_range(i) = mean(sp >= opts.SpeedRange(1) & sp <= opts.SpeedRange(2));
    end
    G.rejects = G.vtest_p < opts.Alpha;
    G.passes = false(height(G), 1);
    for i = find(G.condition == "original")'
        same = G.estimator == G.estimator(i) & (G.param == G.param(i) | (isnan(G.param) & isnan(G.param(i))));
        nul = same & G.condition ~= "original" & G.condition ~= "lagged";
        G.passes(i) = G.rejects(i) && ~any(G.rejects(nul));
    end
    G.lagged_n = zeros(height(G), 1);  G.lagged_excess = NaN(height(G), 1);  G.lagged_p = NaN(height(G), 1);
    for i = find(G.condition == "original")'
        same = P.estimator == G.estimator(i) & (P.param == G.param(i) | (isnan(P.param) & isnan(G.param(i))));
        a = P(same & P.condition == "original", {'subject','fracAP'});  b = P(same & P.condition == "lagged", {'subject','fracAP'});
        b.Properties.VariableNames{2} = 'fracL';  d = innerjoin(a, b, 'Keys', 'subject');
        d = d.fracAP - d.fracL;  d = d(isfinite(d));
        if isempty(d), continue, end
        G.lagged_n(i) = numel(d);  G.lagged_excess(i) = median(d);  G.lagged_p(i) = signrank(d, 0, 'tail', 'right');
    end
    G.beats_lagged = G.lagged_p < opts.Alpha;
    i_write(outDir, G, P);
end

function [lo, hi] = i_clopper(x, n)
    if n == 0, lo = NaN; hi = NaN; return, end
    if x == 0, lo = 0; else, lo = betainv(0.025, x, n - x + 1); end
    if x == n, hi = 1; else, hi = betainv(0.975, x + 1, n - x); end
end

function i_write(od, G, P)
    if isempty(od), return, end
    if ~exist(od, 'dir'), mkdir(od); end
    writetable(G, fullfile(od, 'group.csv'));  writetable(P, fullfile(od, 'participants.csv'));
end

% Author: Diellor Basha, 2026
