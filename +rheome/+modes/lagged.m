function L = lagged(m, R, opts)
% MODES.LAGGED  N3: the best non-travelling two-generator model of a complex mode, as the instrument shows it.
%
%   L = rheome.modes.lagged(m, R)              % m [d x 1] complex mode, R [d x d] tile resolution (operator.R)
%   L = rheome.modes.lagged(m, R, StepDeg=10)
%
% Two extended generators on DISJOINT sets of tiles, each a nonnegative pattern with ONE phase, any amplitude
% ratio and any lag: the source has exactly two phases, m ~ e^{i theta} R (x1 + e^{i delta} x2), x1, x2 >= 0,
% x1 .* x2 = 0. For each theta, delta on a StepDeg grid: nonnegative least squares on [Re; Im] with both
% generators free, each tile then given to the generator holding more of it, and refitted on that split
% (twice); the best |r| is kept. A generator may be bilateral or span many tiles (the frontal and posterior
% generators of rheome.forward.laggedpair, here on tiles); a travelling or rotating mode needs a continuum of
% source phases that two phases blurred by R cannot make. ⚠ Without the disjoint split the model is empty:
% overlapping nonnegative generators at a lag near 180 deg write any mode whose phases span < 180 deg, and
% every planted pattern (standing, travelling, rotating, lagged) read |r| = 1.000.
%
%   L.r_n3        |r| of m to its best two-generator fit (1 = indistinguishable from a lagged pair)
%   L.fit         [d x 1] the fit   L.lag_rad   delta   L.x1 L.x2 [d x 1] the generators' tile weights
%
% Why not one tile per generator: a planted frontal-occipital pair of midline blobs (each over both
% hemispheres) fitted by the best pair of single tiles or adjacent-tile patches reached only |r| = 0.63,
% and the reduction called it travelling.
%
% Author: Diellor Basha, 2026

    arguments
        m (:,1)
        R double
        opts.StepDeg (1,1) double = 10
    end
    m = double(m);  d = numel(m);  Rr = [real(R); imag(R)];  Ri = [-imag(R); real(R)];
    best = -Inf;  op = optimset('TolX', 1e-10);
    L = struct('r_n3', NaN, 'fit', nan(d, 1), 'lag_rad', NaN, 'x1', nan(d, 1), 'x2', nan(d, 1));   % a non-finite mode has no fit
    for th = deg2rad(0:opts.StepDeg:360-opts.StepDeg)
        y = m * exp(-1i*th);  y = [real(y); imag(y)];
        for de = deg2rad(0:opts.StepDeg:180)
            A = [Rr, cos(de) * Rr + sin(de) * Ri];          % R x1 + e^{i delta} R x2, stacked real/imag
            x = lsqnonneg(A, y, op);
            for pass = 1:2                                  % each tile to one generator, refit on the split
                s1 = x(1:d) >= x(d+1:end);
                x = zeros(2*d, 1);  k = [find(s1); d + find(~s1)];  x(k) = lsqnonneg(A(:, k), y, op);
            end
            f = exp(1i*th) * (R * x(1:d) + exp(1i*de) * (R * x(d+1:end)));
            r = abs(f' * m) / max(norm(f) * norm(m), eps);
            if r > best, best = r;  L = struct('r_n3', r, 'fit', f, 'lag_rad', de, 'x1', x(1:d), 'x2', x(d+1:end)); end
        end
    end
end

% Author: Diellor Basha, 2026
