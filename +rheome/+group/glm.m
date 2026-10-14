function R = glm(Y, X, j, opts)
% GROUP.GLM  Tile-wise GLM: one effect's t per tile, its permutation p and the family-wise (max |t|) p over tiles.
%
%   R = rheome.group.glm(Y, X, j)
%   R = rheome.group.glm(Y, X, j, NPerm=5000, Seed=1, Subject=sid)
%
%   Y  [n x m] one value per row (recording) and tile, all finite
%   X  [n x p] the design, the intercept included (e.g. [ones(n,1) group covariates])
%   j  the column of X tested (two-sided): a group indicator or a covariate
%
% The t of column j is read at every tile. Its null is Freedman-Lane (Freedman & Lane 1983; Winkler et al. 2014):
% the residuals of the reduced model (X without column j) are permuted, its fit added back and the full model
% refitted, NPerm times (the identity counted once). The tile's p is (1 + #{|t*| >= |t|}) / (1 + NPerm); the
% family-wise p compares |t| with the max |t*| over all tiles of each permutation, so it holds FWER at any
% spatial dependence between tiles (Nichols & Holmes 2002).
%
% Repeated recordings: Subject (n x 1) averages each subject's rows first and tests the means. For an effect and
% nuisance constant within subject this is the random-intercept mixed model by summary statistics (Holmes &
% Friston 1998), exact when every subject has the same number of rows; X must be constant within subject.
%
% Returns R with .t .beta [1 x m], .p .p_fwer [1 x m], .df, .n (rows tested), .maxnull [NPerm x 1] (max |t*|),
% .untestable [1 x m]: tiles whose residuals are zero to round-off (e.g. a constant tile). Their t is NaN, their
% p and p_fwer are 1, and they are left out of the max (a round-off t would otherwise be of arbitrary size).
%
% ⚠ Rows are exchangeable under the null: one cohort per call (no exchangeability blocks across sites).
%
% See also: rheome.group.spin, rheome.group.calibrate
%
% Author: Diellor Basha, 2026

    arguments
        Y double
        X double
        j (1,1) double {mustBeInteger, mustBePositive}
        opts.NPerm (1,1) double {mustBeInteger, mustBeNonnegative} = 5000
        opts.Seed (1,1) double = 1
        opts.Subject = []
    end
    if size(X, 1) ~= size(Y, 1), error('group:glm:size', 'Y has %d rows, X %d', size(Y, 1), size(X, 1)); end
    if ~all(isfinite(Y), 'all') || ~all(isfinite(X), 'all'), error('group:glm:finite', 'Y and X must be finite'); end
    if ~isempty(opts.Subject)
        [g, ~] = findgroups(opts.Subject(:));
        Xm = splitapply(@(x) mean(x, 1), X, g);
        if any(abs(X - Xm(g, :)) > 1e-12 * max(1, abs(X)), 'all')
            error('group:glm:within', 'X varies within a subject: the summary-statistics model needs X constant per subject');
        end
        Y = splitapply(@(y) mean(y, 1), Y, g);  X = Xm;
    end
    n = size(Y, 1);  df = n - rank(X);
    if df < 1, error('group:glm:df', 'no residual degrees of freedom (n %d, rank %d)', n, rank(X)); end
    Z = X(:, setdiff(1:size(X, 2), j));
    Pi = pinv(X);  cjj = Pi(j, :) * Pi(j, :)';                 % (X'X)^-1 (j,j)
    tf = @(Yk) i_t(Yk, X, Pi, j, cjj, df);
    [R.t, R.beta] = tf(Y);
    rss = sum((Y - X * (Pi * Y)).^2, 1);
    R.untestable = sqrt(rss / n) <= 1e3 * eps * max(abs(Y), [], 1);
    R.t(R.untestable) = NaN;
    if isempty(Z), Hy = zeros(size(Y)); else, Hy = Z * (pinv(Z) * Y); end
    Rz = Y - Hy;  rs = RandStream('mt19937ar', 'Seed', opts.Seed);
    cnt = ones(1, size(Y, 2));  R.maxnull = zeros(opts.NPerm, 1);  at = abs(R.t);  at(R.untestable) = Inf;
    for k = 1:opts.NPerm
        ts = abs(tf(Rz(randperm(rs, n), :) + Hy));  ts(R.untestable) = 0;
        cnt = cnt + (ts >= at - 1e-12);  R.maxnull(k) = max(ts);
    end
    R.p = cnt / (1 + opts.NPerm);
    R.p_fwer = (1 + sum(R.maxnull >= at - 1e-12, 1)) / (1 + opts.NPerm);
    R.p(R.untestable) = 1;  R.p_fwer(R.untestable) = 1;
    R.df = df;  R.n = n;
end

function [t, b] = i_t(Y, X, Pi, j, cjj, df)
    B = Pi * Y;  e = Y - X * B;  b = B(j, :);
    t = b ./ sqrt(sum(e.^2, 1) / df * cjj);
end

% Author: Diellor Basha, 2026
