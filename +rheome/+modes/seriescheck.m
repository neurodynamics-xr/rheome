function [T, S] = seriescheck(fitFiles, fieldsFiles, outDir)
% MODES.SERIESCHECK  Are the series a back-fit recomputed the same as the series the discovery read?
%
%   [T, S] = rheome.modes.seriescheck("fits.txt", "fields.txt", outDir)
%
% rheome.modes.fit recomputes each participant's fields (the discovery's primary variant and seed) and writes
% their SHA-256 per series to modesseries.csv; stage 1 (rheome.modes.fields) wrote the same table for the
% series the discovery clustered. Equal hashes mean the fit worked on numerically identical data, whatever
% commit it ran on. fitFiles / fieldsFiles: the participants' modesseries.csv (paths, or a text file of
% paths), named <subject>_modesseries.csv. A <subject>_provenance.json next to a file gives its commit.
% Series are matched on (subject, source, band, variant); the fit holds the primary variant only, so the
% sensitivity variants of stage 1 are not compared.
%
% Writes outDir/seriescheck.csv (one row per matched or unmatched series: sha256 and n_valid of each side,
% equal), seriescheck_summary.csv (per fit commit x source x band: n equal, n differ, n missing) and
% seriescheck.json (n_subjects, n_series, n_equal, n_differ, n_missing, all_equal, the subjects that differ).
%
% See also: rheome.modes.fit, rheome.modes.fields, rheome.modes.reduce
%
% Author: Diellor Basha, 2026

    arguments
        fitFiles string
        fieldsFiles string
        outDir (1,:) char
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    A = i_read(i_list(fitFiles));  B = i_read(i_list(fieldsFiles));
    B = B(ismember(B.variant, unique(A.variant)), :);                % stage 1's sensitivity variants: not fitted
    key = {'subject', 'source', 'band', 'variant'};
    T = outerjoin(A, B, 'Keys', key, 'MergeKeys', true);
    T = renamevars(T, {'sha256_A', 'sha256_B', 'n_valid_A', 'n_valid_B', 'commit_A', 'commit_B'}, ...
                   {'sha256_fit', 'sha256_fields', 'n_valid_fit', 'n_valid_fields', 'commit_fit', 'commit_fields'});
    T.missing = ismissing(T.sha256_fit) | ismissing(T.sha256_fields);
    T.equal = ~T.missing & T.sha256_fit == T.sha256_fields;
    T.commit_fit = fillmissing(T.commit_fit, 'constant', "unknown");
    writetable(T, fullfile(outDir, 'seriescheck.csv'));
    [g, cm, src, bn] = findgroups(T.commit_fit, T.source, T.band);
    S = table(cm, src, bn, splitapply(@nnz, T.equal, g), splitapply(@(e, m) nnz(~e & ~m), T.equal, T.missing, g), ...
              splitapply(@nnz, T.missing, g), 'VariableNames', {'commit_fit', 'source', 'band', 'n_equal', 'n_differ', 'n_missing'});
    writetable(S, fullfile(outDir, 'seriescheck_summary.csv'));
    j = struct('n_subjects', numel(unique(T.subject)), 'n_series', height(T), 'n_equal', nnz(T.equal), ...
               'n_differ', nnz(~T.equal & ~T.missing), 'n_missing', nnz(T.missing), 'all_equal', all(T.equal), ...
               'subjects_differ', {cellstr(unique(T.subject(~T.equal)))'});
    fid = fopen(fullfile(outDir, 'seriescheck.json'), 'w');  fprintf(fid, '%s', jsonencode(j, PrettyPrint=true));  fclose(fid);
    fprintf('[modes.seriescheck] %d subjects, %d series: %d equal, %d differ, %d missing\n', j.n_subjects, j.n_series, ...
            j.n_equal, j.n_differ, j.n_missing);
end

function f = i_list(f)
    if isscalar(f) && endsWith(f, ".txt"), f = string(strip(splitlines(strtrim(fileread(f))))); f(f == "") = []; end
end

function D = i_read(files)
% every file's series rows with the subject (from the file name) and the commit of its run (provenance.json)
    D = table();
    for f = files(:)'
        [d, n] = fileparts(f);  sub = extractBefore(n + "_modesseries", "_modesseries");
        t = readtable(f, 'TextType', 'string', 'Delimiter', ',');
        t = t(:, {'source', 'band', 'variant', 'n_valid', 'sha256'});
        cm = "unknown";  pj = fullfile(d, sub + "_provenance.json");
        if isfile(pj), p = jsondecode(fileread(pj));  if isfield(p, 'commit'), cm = string(p.commit); end, end   % run.m writes 'commit'
        t.subject = repmat(sub, height(t), 1);  t.commit = repmat(cm, height(t), 1);
        D = [D; t]; %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
