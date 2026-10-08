function G = reduce(inDir, outDir, opts)
% SCALE.REDUCE  Per-subject tables -> group distributions, the anchor's place, cohort contrasts.
%
%   G = rheome.scale.reduce(inDir, outDir)
%   G = rheome.scale.reduce(inDir, outDir, Anchor="sub01", AnchorDir=..., Reference="norm", TwinSide="a_")
%
% inDir holds one folder per subject (rheome.scale.run's OutDir, one per subject), each
% with metrics.csv, timing.csv and optionally job.tsv (twin id). Writes to outDir:
%   metrics_long.csv     every subject's rows stacked
%   status.csv           per analysis: subjects ok / error / not_ported, median seconds and peak GB
%   distribution.csv     per (analysis, metric, band, cohort): n, median, q25, q75, min, max
%   anchor.csv           the anchor's value, its percentile in the Reference cohort, and that
%                        cohort's median and IQR
%   cohorts.csv          per cohort vs Reference: median difference and Cliff's delta (EXPLORATORY,
%                        descriptive; no test, no correction)
%   twins.csv            each subject vs its twin (the same recording, imported by two protocols):
%                        n pairs, Spearman rho, median absolute difference -- a pipeline
%                        reproducibility check, not a cohort comparison
%
% ⚠ TWINS ARE THE SAME PEOPLE. Two cohorts that share recordings are reported as their own cohorts
% and never pooled with the Reference; the twin table is the only place they meet. A twin is read
% from job.tsv ("twin<TAB><id>"); only subjects whose folder name starts with TwinSide list theirs,
% so each pair is counted once (TwinSide "" takes every subject). Anchor "" places no anchor.
%
% See also: rheome.scale.run, rheome.scale.analyses
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char
        opts.Anchor string = ""
        opts.TwinSide string = ""
        opts.AnchorDir char = ''
        opts.Reference string = "norm"
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    dd = dir(inDir);  dd = dd([dd.isdir] & ~startsWith({dd.name}, '.'));
    M = table();  Tm = table();  tw = strings(0,2);
    for k = 1:numel(dd)
        p = fullfile(inDir, dd(k).name);
        if isfile(fullfile(p, 'metrics.csv')), M = [M; i_read(fullfile(p,'metrics.csv'))]; end %#ok<AGROW>
        if isfile(fullfile(p, 'timing.csv'))
            t = readtable(fullfile(p,'timing.csv'), 'TextType','string');
            t.subject = repmat(string(dd(k).name), height(t), 1);  Tm = [Tm; t]; %#ok<AGROW>
        end
        if isfile(fullfile(p, 'job.tsv'))
            j = readlines(fullfile(p, 'job.tsv'));  j = j(startsWith(j, "twin" + char(9)));
            if ~isempty(j), v = extractAfter(j(1), char(9)); if strlength(v) > 0 && startsWith(dd(k).name, opts.TwinSide), tw(end+1,:) = [string(dd(k).name) v]; end, end %#ok<AGROW>
        end
    end
    if ~isempty(opts.AnchorDir) && isfile(fullfile(opts.AnchorDir, 'metrics.csv'))
        A = i_read(fullfile(opts.AnchorDir,'metrics.csv'));
        A.cohort(:) = "anchor";  M = [M; A];
    end
    M.band(ismissing(M.band)) = "";
    M.cohort(M.subject == opts.Anchor) = "anchor";      % placed IN the reference, never counted in it
    writetable(M, fullfile(outDir, 'metrics_long.csv'));
    G.metrics = M;

    if ~isempty(Tm)
        Tm.message(ismissing(Tm.message)) = "";
        [g, an] = findgroups(Tm.analysis);
        G.status = table(an, splitapply(@(s) sum(s=="ok"), Tm.status, g), splitapply(@(s) sum(s=="error"), Tm.status, g), ...
            splitapply(@(s) sum(s=="not_ported"), Tm.status, g), splitapply(@median, Tm.seconds, g), splitapply(@max, Tm.peak_rss_GB, g), ...
            'VariableNames', {'analysis','ok','error','not_ported','median_s','max_peak_GB'});
        writetable(G.status, fullfile(outDir, 'status.csv'));
    end

    key = M.analysis + "|" + M.metric + "|" + M.band;
    [g, k1, c1] = findgroups(key, M.cohort);
    q = @(x,p) prctile(x(~isnan(x)), p);
    D = table(k1, c1, splitapply(@(x) sum(~isnan(x)), M.value, g), splitapply(@(x) q(x,50), M.value, g), ...
              splitapply(@(x) q(x,25), M.value, g), splitapply(@(x) q(x,75), M.value, g), ...
              splitapply(@(x) min(x), M.value, g), splitapply(@(x) max(x), M.value, g), ...
              'VariableNames', {'key','cohort','n','median','q25','q75','min','max'});
    D = [i_split(D.key) D(:,2:end)];
    writetable(D, fullfile(outDir, 'distribution.csv'));  G.distribution = D;

    isA = M.subject == opts.Anchor;  isR = M.cohort == opts.Reference & ~isA;
    keys = unique(key(isA));  An = table();
    for kk = keys(:)'
        a = M.value(isA & key == kk);  r = M.value(isR & key == kk);  r = r(~isnan(r));
        if isempty(a) || isempty(r), continue; end
        An = [An; {kk, a(1), numel(r), 100*mean(r < a(1)) + 50*mean(r == a(1)), median(r), q(r,25), q(r,75)}]; %#ok<AGROW>
    end
    if ~isempty(An)
        An.Properties.VariableNames = {'key','anchor','n_ref','percentile','ref_median','ref_q25','ref_q75'};
        An = [i_split(An.key) An(:,2:end)];  writetable(An, fullfile(outDir, 'anchor.csv'));
    end
    G.anchor = An;

    C = table();  cohorts = setdiff(unique(M.cohort), [opts.Reference "anchor"]);
    for kk = unique(key)'
        r = M.value(isR & key == kk);  r = r(~isnan(r));
        for c = cohorts(:)'
            x = M.value(M.cohort == c & key == kk);  x = x(~isnan(x));
            if numel(x) < 3 || numel(r) < 3, continue; end
            C = [C; {kk, c, numel(x), numel(r), median(x) - median(r), mean(sign(x(:) - r(:)'), 'all')}]; %#ok<AGROW>
        end
    end
    if ~isempty(C)
        C.Properties.VariableNames = {'key','cohort','n','n_ref','median_diff','cliffs_delta'};
        C = [i_split(C.key) C(:,2:end)];  writetable(C, fullfile(outDir, 'cohorts.csv'));
    end
    G.cohorts = C;

    W = table();
    if ~isempty(tw)
        for kk = unique(key)'
            a = nan(size(tw,1),1);  b = a;
            for i = 1:size(tw,1)
                x = M.value(M.subject == tw(i,1) & key == kk);
                y = M.value(endsWith(M.subject, erase(tw(i,2), "-")) & key == kk);
                if ~isempty(x) && ~isempty(y), a(i) = x(1); b(i) = y(1); end
            end
            ok = ~isnan(a) & ~isnan(b);
            if nnz(ok) >= 3
                W = [W; {kk, nnz(ok), corr(a(ok), b(ok), 'Type','Spearman'), median(abs(a(ok) - b(ok)))}]; %#ok<AGROW>
            end
        end
        if ~isempty(W)
            W.Properties.VariableNames = {'key','n_pairs','spearman','median_absdiff'};
            W = [i_split(W.key) W(:,2:end)];  writetable(W, fullfile(outDir, 'twins.csv'));
        end
    end
    G.twins = W;
end

function T = i_read(f)
% ⚠ an all-empty column (band, when a subject has only band-free rows) reads back as double NaN,
% which turns every key <missing> and breaks the vertcat across subjects: force text to string
    o = detectImportOptions(f, 'TextType', 'string', 'Delimiter', ',');
    o = setvartype(o, intersect({'subject','cohort','dataset','analysis','metric','band','unit'}, o.VariableNames), 'string');
    T = readtable(f, o);
    T.band(ismissing(T.band)) = "";  T.cohort(ismissing(T.cohort)) = "";
end

function T = i_split(key)
    p = split(key(:), "|");  if isscalar(key), p = reshape(p, 1, []); end   % one key splits to a column
    T = table(p(:,1), p(:,2), p(:,3), 'VariableNames', {'analysis','metric','band'});
end

% Author: Diellor Basha, 2026
