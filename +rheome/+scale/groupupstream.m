function G = groupupstream(inDir)
% SCALE.GROUPUPSTREAM  Gate G of AGE-1 (TP2) and WAVE-OR-LEAKAGE (TP3): reliability of U* and its spread, and nothing else.
%
%   G = rheome.scale.groupupstream('/path/to/rheome-upstream')    % <sub>_metrics.csv (nsp) or <sub>/metrics.csv (run)
%
% Why: both test plans require the split-half reliability of U* = U - U_surrogate before any outcome is read,
% computed blind (no age, no mean). This reads the "upstream" rows of every recording's metrics and returns
%   .icc          ICC(2,1) (two-way random, absolute agreement, single measure; Shrout & Fleiss 1979) of
%                 U_star_h1 against U_star_h2 over the recordings read
%   .sd           SD of U* on the whole span (the TP3 power check: TOST underpowered at n = 411 if > 0.12)
%   .n_read .n_not_read   recordings read / not read (span or kept weight; their reasons in upstream.csv)
%   .n_aperiodic_failed   recordings read whose aperiodic fit failed on some channel (gate P)
%   .passed       icc >= 0.5
% ⚠ IT RETURNS NO MEAN, BY DESIGN: not of U*, U, U_total or Delta_ap. Gate G runs only after TP2 and TP3 are
%   approved; the outcome tests read the same files afterwards.
%
% See also: rheome.scale.measure_upstream
%
% Author: Diellor Basha, 2026

    f = dir(fullfile(inDir, '**', '*metrics.csv'));
    if isempty(f), error('scale:groupupstream:empty', 'no metrics.csv under %s', inDir); end
    v = nan(numel(f), 5);                                   % read, U*_h1, U*_h2, U*, aperiodic_ok
    names = ["read" "U_star_h1" "U_star_h2" "U_star" "aperiodic_ok"];
    for k = 1:numel(f)
        M = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string');
        b = string(M.band);  b(ismissing(b)) = "";             % an all-empty column reads as NaN
        M = M(M.analysis == "upstream" & b == "", :);
        for j = 1:numel(names)
            x = M.value(M.metric == names(j));  if isscalar(x), v(k, j) = x; end
        end
    end
    has = isfinite(v(:, 1));  rd = has & v(:, 1) == 1;
    Y = v(rd, 2:3);  Y = Y(all(isfinite(Y), 2), :);
    G = struct('icc', i_icc21(Y), 'sd', std(v(rd & isfinite(v(:, 4)), 4)), 'n_read', nnz(rd), ...
               'n_not_read', nnz(has & v(:, 1) == 0), 'n_aperiodic_failed', nnz(rd & v(:, 5) ~= 1), 'n_icc', size(Y, 1));
    G.passed = G.icc >= 0.5;
end

function r = i_icc21(Y)
% ICC(2,1) of an [n x k] table (subjects x raters)
    [n, k] = size(Y);  r = NaN;  if n < 2, return; end
    m = mean(Y, 'all');
    msr = k * sum((mean(Y, 2) - m).^2) / (n - 1);
    msc = n * sum((mean(Y, 1) - m).^2) / (k - 1);
    E = Y - mean(Y, 2) - mean(Y, 1) + m;  mse = sum(E.^2, 'all') / ((n - 1) * (k - 1));
    r = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n);
end

% Author: Diellor Basha, 2026
