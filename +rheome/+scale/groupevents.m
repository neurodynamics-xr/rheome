function [G, A, TF] = groupevents(inDir, outDir, opts)
% SCALE.GROUPEVENTS  MS1 event statistics across participants, per cohort: medians, CIs, split-half ICC, the boundary acceptance and the group TF maps.
%
%   [G, A, TF] = rheome.scale.groupevents(inDir, outDir)
%   [G, A, TF] = rheome.scale.groupevents(inDir, outDir, Boot=2000, Seed=7, Accept=0.8, AcceptSNR=10)
%
% inDir holds rheome.scale.run's per-participant folders (<sub>/metrics.csv, <sub>/eventtf.mat) or nsp's flat
% files (<sub>_metrics.csv, <sub>_eventtf.mat), at any depth. Rows of analyses "events" and "eventplants".
% The cohort is the metrics' `dataset` column (`cohort` when empty); no number mixes cohorts (SPEC1 §5).
% Writes to outDir:
%   group.csv          per cohort x analysis x metric x band: n, median, its 95% percentile bootstrap CI over
%                      participants (lo, hi) and, where the measure has split halves ("<band>_h1"/"_h2"), the
%                      split-half ICC(2,1) (Shrout & Fleiss 1979) and its n
%   acceptance.csv     per cohort x "<snr>dB/tau<bound>/<model>": participants, the share passing (eventplants
%                      participant_pass) and accepted = share >= Accept. ⭐ SPEC1 §8 A's rule reads the row
%                      "<AcceptSNR>dB/tau1/<Model>" (measure_events' defaults): if not accepted, MS1 reports
%                      the tracker's window boundaries
%   group_eventtf.mat  per cohort (TF.<cohort>): f, t, n and, per alignment (on_fit off_fit on_window
%                      off_window), the median dB map over participants and its pointwise bootstrap CI
%                      (<a>, <a>_lo, <a>_hi)
%   participants.csv   the participant rows read
%
% See also: rheome.scale.measure_events, rheome.scale.measure_eventplants, rheome.scale.groupslowosc
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char = ''
        opts.Boot (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.BootTF (1,1) double {mustBeInteger, mustBePositive} = 500
        opts.Seed (1,1) double = 7
        opts.Accept (1,1) double = 0.8
        opts.AcceptSNR (1,1) double = 10
        opts.Model (1,1) string = "posterior"
    end
    rng(opts.Seed);
    f = dir(fullfile(inDir, '**', '*metrics.csv'));  P = table();
    for k = 1:numel(f)
        M = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string', 'Delimiter', ',');
        if isempty(M) || ~ismember('analysis', M.Properties.VariableNames), continue, end
        M = M(ismember(M.analysis, ["events" "eventplants"]), :);
        if isempty(M), continue, end
        c = M.dataset;  c(ismissing(c) | c == "") = M.cohort(ismissing(c) | c == "");
        M.cohortKey = c;  P = [P; M(:, {'subject','cohortKey','analysis','metric','band','value','unit'})]; %#ok<AGROW>
    end
    G = table();  A = table();  TF = struct();
    if isempty(P), i_write(outDir, G, A, TF, P); return, end
    P.subject = string(P.subject);
    [~, ia] = unique(P(:, {'subject','analysis','metric','band'}), 'last');  P = P(sort(ia), :);   % a re-run counts once
    isH = endsWith(P.band, ["_h1" "_h2"]);
    [gk, G] = findgroups(P(~isH, {'cohortKey','analysis','metric','band','unit'}));
    Q = P(~isH, :);
    for i = 1:height(G)
        v = Q.value(gk == i);  v = v(isfinite(v));  G.n(i) = numel(v);  G.median(i) = median(v);
        if numel(v) > 1
            ci = prctile(median(v(randi(numel(v), numel(v), opts.Boot)), 1), [2.5 97.5]);
        else, ci = [NaN NaN]; end
        G.lo(i) = ci(1);  G.hi(i) = ci(2);
        a = P(P.cohortKey == G.cohortKey(i) & P.analysis == G.analysis(i) & P.metric == G.metric(i), :);
        h1 = a(a.band == G.band(i) + "_h1", {'subject','value'});  h2 = a(a.band == G.band(i) + "_h2", {'subject','value'});
        G.icc(i) = NaN;  G.icc_n(i) = 0;
        if ~isempty(h1) && ~isempty(h2)
            j = innerjoin(h1, h2, 'Keys', 'subject');  Y = [j.value_h1 j.value_h2];  Y = Y(all(isfinite(Y), 2), :);
            G.icc(i) = i_icc21(Y);  G.icc_n(i) = size(Y, 1);
        end
    end
    a = G(G.analysis == "eventplants" & G.metric == "participant_pass", :);
    A = table(a.cohortKey, a.band, a.n, a.median, 'VariableNames', {'cohort','band','n','x'});
    for i = 1:height(A)
        v = P.value(P.cohortKey == A.cohort(i) & P.analysis == "eventplants" & P.metric == "participant_pass" & P.band == A.band(i));
        A.share_pass(i) = mean(v == 1);  A.accepted(i) = A.share_pass(i) >= opts.Accept;
    end
    A.x = [];  A.declared_rule = A.band == sprintf('%gdB/tau1/%s', opts.AcceptSNR, opts.Model);

    % group TF maps: subject -> cohort from the metrics
    sc = unique(P(:, {'subject','cohortKey'}));
    ft = dir(fullfile(inDir, '**', '*eventtf.mat'));  al = ["on_fit" "off_fit" "on_window" "off_window"];  Z = struct();
    for k = 1:numel(ft)
        nm = string(ft(k).name);
        if nm == "eventtf.mat", [~, sub] = fileparts(ft(k).folder); else, sub = extractBefore(nm, "_eventtf.mat"); end
        c = sc.cohortKey(sc.subject == sub);  L = load(fullfile(ft(k).folder, ft(k).name));
        if isfield(L, 'TF'), L = L.TF; end
        if isempty(c) || ~isfield(L, 'on_window'), continue, end
        c = matlab.lang.makeValidName(c(1));
        if ~isfield(Z, c), Z.(c) = struct('f', L.f, 't', L.t, 'subjects', strings(0, 1)); for x = al, Z.(c).(x) = []; end, end
        Z.(c).subjects(end + 1, 1) = sub;
        for x = al, Z.(c).(x) = cat(3, Z.(c).(x), L.(x)); end
    end
    for c = string(fieldnames(Z))'
        TF.(c) = struct('f', Z.(c).f, 't', Z.(c).t, 'n', numel(Z.(c).subjects), 'subjects', Z.(c).subjects);
        for x = al
            M3 = Z.(c).(x);  n = size(M3, 3);  TF.(c).(x) = median(M3, 3, 'omitnan');
            if n > 1
                B = zeros([size(M3, 1:2) opts.BootTF]);
                for q = 1:opts.BootTF, B(:, :, q) = median(M3(:, :, randi(n, n, 1)), 3, 'omitnan'); end
                TF.(c).(x + "_lo") = prctile(B, 2.5, 3);  TF.(c).(x + "_hi") = prctile(B, 97.5, 3);
            end
        end
    end
    i_write(outDir, G, A, TF, P);
end

function i_write(od, G, A, TF, P)
    if isempty(od), return, end
    if ~exist(od, 'dir'), mkdir(od); end
    writetable(G, fullfile(od, 'group.csv'));  writetable(A, fullfile(od, 'acceptance.csv'));
    writetable(P, fullfile(od, 'participants.csv'));  save(fullfile(od, 'group_eventtf.mat'), 'TF', '-v7');
end

% ICC(2,1) of an [n x k] table (subjects x raters), as rheome.scale.groupupstream
function r = i_icc21(Y)
    [n, k] = size(Y);  r = NaN;  if n < 2, return; end
    m = mean(Y, 'all');
    msr = k * sum((mean(Y, 2) - m).^2) / (n - 1);
    msc = n * sum((mean(Y, 1) - m).^2) / (k - 1);
    E = Y - mean(Y, 2) - mean(Y, 1) + m;  mse = sum(E.^2, 'all') / ((n - 1) * (k - 1));
    r = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n);
end

% Author: Diellor Basha, 2026
